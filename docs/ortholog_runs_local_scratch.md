# Why the ortholog jobs compute on node-local disk

Each ortholog pair is computed in a temporary directory on the compute node's own disk, and only
the small results are copied to `ortholog_runs/ortholog_distributions/`. The BLAST table, which is
almost everything the step writes, is discarded with that directory.

This is not a performance tweak. It is there because writing it to `the shared group directory` silently loses
data.


## The problem

`the shared group directory` is NFSv3, mounted on the compute nodes with

    vers=3, soft, proto=tcp, timeo=600, retrans=2, rsize=32768, wsize=32768

`soft` is the whole story. When the storage server does not answer, the client retries the write
twice over about three minutes and then **abandons the data**. A `hard` mount would retry
indefinitely, so the job stalls until the server recovers and nothing is lost.

Two further things turn that into silent corruption rather than a visible error:

- **The error reaches a `write()` nobody is checking.** blastp writes through buffered I/O, so the
  bytes go to the page cache and the RPC happens later during writeback, long after `write()`
  returned success. The error is attached to the file and reported to whichever `close()` or
  `fsync()` picks it up first — once, to one file descriptor. blastp, like most tools, does not
  check the result of `close()`.
- **The file keeps its correct length.** Writes at higher offsets succeed, so the size is right.
  The lost range simply has no data on the server and reads back as zeros.

The result is a BLAST table of exactly the expected size with one contiguous run of NUL bytes
through the middle of it.

### How it surfaces

Sometimes loudly. If the NULs land inside a sequence identifier, pandas reads the whole NUL block
as a single field — it contains no tab or newline — and that value becomes a filename:

    File ".../wgd_ksrates/blast_mcl.py", line 108, in get_one_v_one_orthologs_rbh
      with open(os.path.join(output_dir, '{}.ovo.tsv'.format(comb)), ...
    ValueError: embedded null byte

That is the lucky case. If the NULs land in a score or an alignment length instead, only those
lines are unparseable: the run finishes normally and the Ks distribution is quietly built on less
data, with nothing anywhere to say so.

A useful diagnostic detail from the 6 November incident: the NUL run in each ksrates log was
*exactly* the same length as the NUL run in its BLAST table (385794, 348779, 318951, 212301
bytes), because the corrupted block was read back as one field and interpolated into the
`Writing one-to-one orthologs to ...` message. The logs were never separately corrupted — there
is one corruption event per pair, in the BLAST table only.

### Observed rate

| date | concurrent blastp | corrupted | rate |
|---|---|---|---|
| 2026-11-04 | ~30 | 3 of 45 tables | 6.7% |
| 2026-11-06 | 8 (`maxForks`) | 4 of 158 tables | 2.5% |

Lowering concurrency helps and does not solve it. Extrapolated to the full run of 9101 ortholog
pairs, 2.5% is roughly 230 corrupted pairs — some crashing visibly, the rest degrading results
silently. On 6–7 November the four failures also took three `runKsratesPipeline` jobs with them:
the missing pairs were absent from the Ks database, so each pipeline recomputed them itself, which
cost nearly two extra hours and pushed all three past their four-hour wall (exit 140).


## Why moving the working directory is enough

ksrates writes and reads `ortholog_distributions/` **relative to its working directory**, in both
commands:

- `ksrates/wgd_orthologs.py`: `ortholog_dists_dir = os.path.join("ortholog_distributions", "")`
- `ksrates/compute_peaks.py` and `ksrates/fc_kde_bootstrap.py`:
  `os.path.join("ortholog_distributions", f"wgd_{sp1}_{sp2}", f"{sp1}_{sp2}.ks.tsv")`

No absolute path is involved, so changing where the job runs moves the whole ortholog analysis
without touching ksrates. Nothing in the wrapper reads these files at all — `grep` for
`ortholog_distributions`, `blast.tsv` or `ks.tsv` across `batch_ksrates/` and `run_ksrates_raw/`
returns nothing.

Everything that has to persist goes to the **Ks database server**, by an absolute path from the
config (`ks_database_server_address_path`, with `use_ks_database_server = yes`). In
`fc_kde_bootstrap.estimate_peak`, the TSV fallback files are only written in the `else` branch,
when no server address is given. So the peak and the Ks list leave the node regardless of what
happens to the working directory.


## What is kept and what is thrown away

Measured over the 116 pairs finished on 6–7 November:

| file | total | per pair | kept? |
|---|---|---|---|
| `.blast.tsv` | 21.15 GB | 182 MB average, 541 MB max | **no** |
| `.ks.tsv` | 0.12 GB | ~1 MB | yes — the peak comes from it |
| `.orthologs.tsv` | small | ~200 KB | yes |

The BLAST table is 99.4% of the bytes and is the file that corrupts. Over 9101 pairs, discarding
it keeps about 1.6 TB off the share.

Its only use would be letting a rerun of the same pair skip BLAST: `fc_wgd.ks_orthologs` re-runs
everything when `.blast.tsv` is absent. That never happens in practice, because a pair is only
queued when its peak is **missing from the Ks database**, and a successful job puts it there. To
keep them anyway, add `ortholog_distributions/wgd_*/*.blast.tsv` to the copy-back loop in
`runOrthologJob`.

The copy-back is also the error check the job was missing. `cp` reports a failed `close()`, which
is exactly what blastp does not do, so with `set -e` a lost write during the copy fails the job
instead of leaving a holed file behind.


## Choosing the scratch directory

**"Node-local" cannot be inferred from a variable name.** This was learned the expensive way. The
first version of this change wrote to `${TMPDIR:-/tmp}`, on the assumption that a scheduler's
`$TMPDIR` is node-local. On this cluster SLURM sets

    TMPDIR=/scratch/tmp/<user>
    /scratch/tmp -> <a second file server>:/mnt/datastore/tmp   nfs, soft

so the BLAST tables moved from one soft-mounted NFS share to another, and kept corrupting. The
rate barely moved - 2.5% of pairs before the change, 1.9% after - and it went unnoticed for a
week, because the check being run was "did a table appear under `ortholog_runs`", which was false
either way. Hard-coding `/tmp` instead is no better: on at least one login node here `/tmp` is itself
NFS, while on the compute nodes it is local XFS with 70 GB.

So `runOrthologJob` asks the filesystem rather than trusting a name. At job start it walks
`$SCRATCHDIR`, `$TMPDIR`, `/tmp`, `/var/tmp` and takes the first that is

1. writable,
2. **not** on a network filesystem - `stat -f -c %T` reports the type as a word (`xfs`, `nfs`,
   `tmpfs`, `lustre`), and `nfs*`, `cifs`, `smb*`, `lustre`, `gpfs`, `fuse*`, `afs`, `9p` and
   `ceph` are rejected,
3. and has at least `ortholog_scratch_min_gb` free - which also rejects a small tmpfs `/tmp`,
   where a 500 MB table would be spent out of the job's memory.

**Every job logs what it picked**, as its first line:

    scratch directory: /tmp (xfs, 56 GB free)

That line is the point. The original mistake was invisible; this makes it the first thing in
every one of the logs.

If nothing qualifies, the job warns loudly and falls back to `${TMPDIR:-/tmp}` rather than
failing - clusters with no node-local disk exist, and a network scratch mounted `hard` is safe.
The warning names the filesystem type so the reader can judge.

## Settings

`ortholog_scratch_dir` in `nextflow.config`, default `null`:

- `null` runs the probe above.
- Any other value is **used as given**, never overridden - but the job still reports its
  filesystem type and warns if it is networked or short of space. Silently ignoring an explicit
  setting would repeat the original bug in a different place.
- A path other than `/tmp` or `$HOME` must also appear in `container_mounts`. Singularity binds
  those two by itself and nothing else.

`ortholog_scratch_min_gb`, default 10: the free space a candidate must have. Budget about
**1.5 GB per concurrent pair** on the node - the table plus the `blast_tmp` working directory it
is assembled from - and remember that several pairs can land on the same node.

A `trap ... EXIT` removes the directory whether the job succeeds, fails or is killed, so a 500 MB
table is never left behind on a compute node; a sweep at job start clears anything older than 48
hours that a SIGKILL or a node crash left behind. Measured over the 524-pair fern shard: no
scratch directory survived.


## This does not fix the mount

Local scratch removes the concurrent writes that trigger the loss; it does not make `the shared group directory`
safe. Anything else writing large files from many nodes at once — including the copy-back itself,
in principle — is still exposed. `check_blast_tables.sh` remains worth running as a gate between
the ortholog stage and the pipeline stage.

The real fix is a `hard` mount, which is a sysadmin decision affecting every user of the share:
with `hard`, a storage outage makes processes block in uninterruptible sleep instead of failing,
which can leave a node needing a reboot. That is a worse failure mode for interactive users and a
much better one for this analysis, so it is not a change to ask for unilaterally. Lower-impact
alternatives worth raising instead: a larger `retrans`, or a second `hard` mount of the same
export at a different path (`nosharecache`) for this project only.
