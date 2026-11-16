#!/bin/bash
# Find Ks pipeline output tables that the NFS share corrupted while they were being written.
#
# The shared group directory is NFSv3 mounted `soft`: a write RPC that exceeds timeo x retrans is abandoned rather than retried, the
# application is not reliably told (the error surfaces during writeback, long after write() returned), and the file
# keeps its correct size with a run of NUL bytes where the lost data should be. Measured rates: 3 of 45 tables at
# ~30 concurrent blastp, still 4 of 158 after halving the concurrency. Full background, and why this is not fixed by
# running fewer jobs: code/1kp_wgms/docs/ortholog_runs_local_scratch.md
#
#   check_nfs_corruption.sh [directory ...]      (default: datasets and ortholog_runs under BASE)
#
# Prints one line per corrupted table and exits 1 if there are any. It never deletes anything.
#
# WHAT THIS IS FOR NOW. It used to be the only thing standing between a corrupted blast table and a quietly wrong Ks
# distribution, and was meant to gate the pipeline stage. Since runOrthologJob computes each pair on node-local disk
# (main.nf), the blast table never crosses NFS at all and the copy-back is checked by cp, which - unlike blastp -
# reports a failed close(). So this is a safety net to run after a batch, not a gate, and two findings mean
# different things:
#
#   ortholog_runs/ortholog_distributions/**.ks.tsv
#       The archived copy of the Ks values. A corrupted one CANNOT have poisoned the Ks database: the peak is
#       computed on the compute node from the local file, before anything is copied back. Only the archive is
#       damaged. Recompute the pair if you want the file back, otherwise it is cosmetic.
#
#   datasets/**/*.blast.tsv  or  datasets/**/*.ks.tsv
#       Written by a NESTED ksrates pipeline, which only runs its own wgdOrthologs when a pair is MISSING from the
#       Ks database. Finding anything here is a signal in itself: the ortholog stage left a gap. These feed the peak
#       directly, so a corrupted one can reach the database. Recompute the pair.
#
# A corrupted pair is redone by deleting its wgd_<a>_<b> directory and running the wrapper again (without -resume
# after a change to a process script, or every pair re-runs - see launch_batch_ksrates.sh).

set -uo pipefail

PYTHON=/software/shared/apps/x86_64/python/3.11.4/bin/python3
BASE=/group/esb/cesen/1kp/ks_analysis/1kp

# The tables that cross NFS and whose content is read back as data. Extend if a stage starts writing another one.
PATTERNS=(-name '*.blast.tsv' -o -name '*.ks.tsv' -o -name '*.ks_anchors.tsv')

DIRS=("$@")
if [ ${#DIRS[@]} -eq 0 ]; then
    DIRS=("$BASE/datasets" "$BASE/ortholog_runs")
fi

existing=()
for dir in "${DIRS[@]}"; do
    [ -d "$dir" ] && existing+=("$dir")
done
if [ ${#existing[@]} -eq 0 ]; then
    echo "none of the given directories exist: ${DIRS[*]}" >&2
    exit 2
fi

# One python process for the whole scan: a run can hold thousands of tables, and a process per file is most of the
# wall time. Reads in chunks, so a 541 MB blast table never lands in memory.
find "${existing[@]}" \( "${PATTERNS[@]}" \) -type f -print0 2>/dev/null | "$PYTHON" -I -c '
import sys

CHUNK = 1 << 22
checked = bad = 0
findings = []

for path in sys.stdin.buffer.read().split(b"\0"):
    if not path:
        continue
    path = path.decode("utf8", "surrogateescape")
    try:
        with open(path, "rb") as fh:
            total = 0
            first = None      # offset of the first NUL, to show the hole
            offset = 0
            while True:
                chunk = fh.read(CHUNK)
                if not chunk:
                    break
                n = chunk.count(0)
                if n:
                    total += n
                    if first is None:
                        first = offset + chunk.index(0)
                offset += len(chunk)
    except OSError as err:
        print(f"UNREADABLE  {path}  ({err})")
        bad += 1
        continue
    checked += 1
    if total:
        bad += 1
        findings.append((path, total, first, offset))

for path, total, first, size in findings:
    print(f"CORRUPT  {path}")
    print(f"         {total} NUL bytes from offset {first}, file is {size} bytes")

print()
print(f"checked {checked} tables, {bad} corrupted")
sys.exit(1 if bad else 0)
'
status=$?

if [ "$status" -ne 0 ]; then
    echo
    echo "To redo a corrupted pair, remove its wgd_<a>_<b> directory and run the wrapper again:"
    echo "  rm -rf <the directory holding the table printed above>"
    echo "Deleting is left to you on purpose - check the path first."
fi
exit "$status"
