# `make_branch_subtrees.py`

This script creates the tables used to set up branch-centered `ksrates` analyses:
`focal_branch_to_dataset.tsv` (one row per branch) and `subtree_datasets.tsv` (one row per distinct dataset).
See `dataset_generation_summary.md` for a one-page description of how the datasets are built. It reads the rooted 1KP species tree, visits every non-root
branch, and writes a small local subtree around that branch.

## How the Local Subtree Is Defined

For each branch, the node below the branch is called `A`.
The script then walks from `A` toward the root and records progressively older
sister clades:

```text
(((A, B), C), D)
```

The retained species are:

- `A_species`: up to four descendants of the test branch.
- `B_species`: up to two species from the sister clade of `A`.
- `C_species`: up to two species from the next older sister clade.
- `D_species`: up to two species from the next older sister clade, if present.

For `A`, the rule is based on the two child subclades below the test node:

- leaf branch: keep the single leaf;
- two descendant leaves: keep both;
- larger clades: keep up to two species from each child subclade.

When `A` is an internal node, its two child subclades can be thought of as
`A'` and `A''`. If `A_species` contains three or four leaves, those leaves are
not just the best three or four species anywhere below `A`; they are sampled
from both `A'` and `A''` so that the descendant side of the test branch remains
phylogenetically balanced.

When there are more candidate leaves than slots, representatives are chosen by
phylogenetic balance first, then by transcriptome quality.

## Two selection strategies

The script implements both, chosen with `--selection`:

| | `paranome` (default) | `quality` |
|---|---|---|
| ranks species by | median bin height of the species' own paranome Ks distribution | `0.8 x BUSCO percentile + 0.2 x TransRate percentile` |
| measures | whether the Ks distribution can carry the analysis | assembly completeness |
| extra constraint on B/C/D | expected-Ks budget from the concatenation-tree path distance | none |
| reports reliability | yes, per clade and per branch | no (the columns are left empty) |
| needs | `others/calibration/calib_common.py` and `paralog_ks_bin_heights.tsv` | `--score-table`, `--genome-busco-table`, `--busco-weight`, `--transrate-weight` |

`paranome` is the current strategy and is described in full in
[`ortholog_ks_path_distance_calibration.md`](ortholog_ks_path_distance_calibration.md). `quality`
is the earlier one, kept because the decision to drop it has not been taken; the section below
documents it and the code is unchanged.

Options specific to `paranome`: `--calibration-dir` (where `calib_common.py` lives, default
`../others/calibration`), `--bin-heights` (default: the table in the calibration data
directory), `--budget-floor` (default 2.0) and `--budget-delta` (default 0.5).

## Quality Score

Used by `--selection quality`. The script reads the 1KP species table:

```text
/path/to/1kp/source_data/1.species_dataset/1kp_paper_2019_suptab1_species.tsv
```

and scores transcriptomes as:

```text
quality_score = 0.8 * BUSCO_percentile + 0.2 * TransRate_percentile
```

Higher scores are preferred when choosing among leaves within the same selected
subclade.

## Genome BUSCO scores (methods)

The 31 genome-based species are not in the 1KP score table, and TransRate cannot be computed for them (it needs sequencing reads). To give them a quality score comparable with the transcriptomes, BUSCO completeness was measured for them as close as possible to the 1KP protocol.

**What 1KP did** (Leebens-Mack et al. 2019, *Nature* 574:679-685, Methods and Extended Data Fig. 3; Carpenter et al. 2019, *GigaScience* 8:giz126): BUSCO with default settings, using the eukaryote (`eukaryota_odb9`) and embryophyte (`embryophyta_odb9`) gene sets, on the proteins predicted from each assembly. BUSCO completeness is reported as "complete + fragmented" (= 100% - missing), and samples with more than 42.5% missing eukaryota BUSCOs were removed. The BUSCO version is not stated unambiguously in the sources we could read. The BUSCO column used here (`% BUSCOs (complete+fragmented)`) is taken to be the `eukaryota_odb9` set, because its values (30-98%, mean 87%) are in the range of that 303-gene set and not of the 1440-gene embryophyta set; the paper does not say so explicitly.

**What was done for the genomes**

1. Input: the CDS FASTA of the primary transcripts of each genome (`genome_species_manifest.tsv`, column `cds_filepath`), after the FASTA/GFF post-processing described in the genome log.
2. Translation of each CDS to protein with the standard genetic code (trailing stop removed, internal stops and ambiguous codons written as `X`), one protein file per species. Script: `others/calculate_busco_for_genomes/translate_cds_to_protein.py` (output in `genome_busco/translated_fasta/`).
3. BUSCO 3.0.2b (cluster module `busco/x86_64/3.0.2b`, protein mode `-m proteins`, default settings, `-c 4`) against `eukaryota_odb9` (303 BUSCOs). BUSCO 3 is the latest version that reads the odb9 gene sets (BUSCO 4 and later require odb10), which are needed to match 1KP. The odb9 data are no longer distributed by the BUSCO authors and were taken from the archived copy at https://zenodo.org/records/20987183. Script: `others/calculate_busco_for_genomes/run_busco_genomes.sh` (one SLURM array task per species).
4. Collection: `others/calculate_busco_for_genomes/collect_busco_scores.py` reads the BUSCO short summaries and writes `genome_busco_scores.tsv` with `busco_complete_fragmented = 100 - missing`, plus the complete, fragmented and missing percentages and the number of BUSCOs (303 for all species).

**How the score is used.** The genome BUSCO values are ranked together with the 1KP transcriptome values, giving each genome a BUSCO percentile. The quality score of a genome is this percentile alone, since it has no TransRate (for transcriptomes: 0.8 x BUSCO percentile + 0.2 x TransRate percentile). TransRate measures read-mapping errors of an assembly and is not defined for an annotated genome, so it cannot be used. This is a convention, not a claim that a genome would score like a transcriptome on TransRate: the BUSCO weight is rescaled from 0.8 to 1 so that the only available measure carries the whole score. A genome and a transcriptome are therefore compared on BUSCO alone for the genome and on BUSCO plus TransRate for the transcriptome, and a transcriptome with a high BUSCO but a low TransRate percentile scores lower than a genome with the same BUSCO percentile.

**Limitations to state when publishing**

- BUSCO 3.0.2b is not necessarily the version 1KP used, and BUSCO was applied to annotated gene models of genomes, not to proteins predicted from assemblies; the values are comparable in kind, not identical in method.
- Most genomes score above the highest transcriptome (98.3%), so they are preferred as representatives whenever a clade offers one. Pinta (62.7%) and Chlva (75.2%) rank low.
- The score only decides which species represent a clade in a dataset; it does not enter the Ks analysis itself.

## Output Columns: `focal_branch_to_dataset.tsv`

One row per branch.

- `dataset_id`: first 8 characters of the SHA-1 of the sorted selected species. It depends only on the species set, so different branches that select the same species get the same id.
- `traversal_order`: position of the branch in the preorder walk of the tree (unique per branch).
- `newick_tree`: pruned Newick tree containing only the selected species.
- `focal_species`: species that should be treated as focal species by downstream `ksrates` setup.
  This is `A_species + B_species + C_species`.
- `target_species`: all selected species in the subtree, including `D_species`.
- `A_species`, `B_species`, `C_species`, `D_species`: selected representatives by local clade.
- `n_selected_species`: total number of species retained in the row.
- `B_path`, `C_path`, `D_path`: median concatenation-tree path distance from the A species to the
  closest representative chosen for that clade. Empty under `--selection quality`.
- `B_reliability`, `C_reliability`, `D_reliability`: `reliable` (path <= 1.0), `borderline` (<= 1.5) or
  `unreliable`, for that same closest representative - the pair whose ortholog Ks peak anchors the
  rate adjustment.
- `rate_adjustment_reliability`: the best of the three, since one usable outgroup is enough to
  anchor an adjustment; `none` if no clade yielded a candidate.

A flagged branch is still analysed in full. The flag belongs to the *rate adjustment*, not to the
WGM call: clade A's own paranome Ks distribution is unaffected by how distant the outgroup is, so
the caveat is "branch assignment uncertain", not "WGM call invalid". Over the 1KP tree the chosen
representative is `reliable` for 87.3% of the 7054 B/C/D slots, `borderline` for 6.3% and `unreliable`
for 6.4%.

## Output Columns: `subtree_datasets.tsv`

One row per distinct `dataset_id`; this is the input of `prepare_branch_ksrates_configs.py`. Branches with the same `dataset_id` have the same tree and target species, so their ksrates runs are identical and they are generated once. For example the two sister branches below a node can select the same species set, with A and B swapped.

- `dataset_id`, `newick_tree`, `target_species`, `n_selected_species`: as above.
- `focal_species`: union of the focal species of the branches that share the dataset.
- `branches`: the `traversal_order` values of those branches.

To interpret a branch, look up its `dataset_id` in `focal_branch_to_dataset.tsv`, then use the output of the focal species of its `A_species` under `datasets/<dataset_id>/`. The script stops with an error if two branches with the same `dataset_id` have different trees or species.

The largest expected subtree has ten species: four from `A` and two each from
`B`, `C`, and `D`.

## Focal Species Logic

Only `A`, `B`, and `C` representatives are focal species by default:

- `A_species`: expected to share a putative WGM on the tested branch.
- `B_species`: closest species expected not to share that WGM.
- `C_species`: older comparison that gives more depth to the Ks plot, to see how the signal relates to the putative WGM on the tested branch.

`D_species` is kept in the subtree as older context for correction/outgroup
choice, but is not treated as focal by default.

## Run

Current strategy (the default; no BUSCO/TransRate options are needed):

```bash
python3 code/1kp_wgms/run_ksrates_raw/make_branch_subtrees.py \
  --tree TREE.nwk \
  --max-a 4 --max-older 2 --older-clades 3 --min-species 3 \
  --out focal_branch_to_dataset.tsv --datasets-out subtree_datasets.tsv
```

Earlier strategy, unchanged:

```bash
python3 code/1kp_wgms/run_ksrates_raw/make_branch_subtrees.py \
  --tree TREE.nwk --selection quality \
  --score-table SPECIES_TABLE.tsv --genome-busco-table GENOME_BUSCO.tsv \
  --max-a 4 --max-older 2 --older-clades 3 --min-species 3 --busco-weight 0.8 --transrate-weight 0.2 \
  --out focal_branch_to_dataset.tsv --datasets-out subtree_datasets.tsv
```

The subtree-shape options (`--max-a`, `--max-older`, `--older-clades`, `--min-species`) have no default: they are parameters of the wrapper, in `batch_ksrates/nextflow.config`, and `main.nf` passes them. `--selection` and the budget options do have defaults, so an existing wrapper invocation that passes the BUSCO options keeps working only if it also passes `--selection quality`; otherwise it now runs the paranome strategy and ignores them. The same holds for `generate_species_metadata.py` and `make_focal_colorings.py`. The example above shows how to run the script by hand.

The script requires `ete3` in the active Python environment.

## The Ks-based scoring system (`--selection paranome`, the default)

The quality score above (`0.8 x BUSCO percentile + 0.2 x TransRate percentile`) is a property of the
*assembly*. What the dataset needs is a representative whose *Ks distribution* can carry the analysis. The
scoring system described in this section measures that directly, from the paranome Ks data already in the Ks
database and from the branch lengths of the concatenation tree.

**This is now implemented and is the default.** The calibration it rests on - the power law, the
two reliability thresholds and the budget - is maintained in
`others/calibration/calib_common.py` and documented in
[`ortholog_ks_path_distance_calibration.md`](ortholog_ks_path_distance_calibration.md), which
supersedes the numbers quoted below where they differ: the conversion is the power law
`Ks = (path/0.4531)^1.266`, not the single factor 1.85, and the hard Ks 2 cut-off has become the
adaptive budget plus a reliability flag.

### Rationale

Two different things are asked of a representative, so two criteria are used:

- **Clade A** supplies the species whose paranome Ks distribution must show (or fail to show) the putative WGM
  on the tested branch. What limits that is the number of duplication events available in the Ks range where a
  peak would sit, i.e. statistical power.
- **Clades B, C and D** supply the species against which the ortholog Ks divergences are estimated. B and C are
  also focal species in their own right, and D anchors the rate adjustment, so all of them need a readable
  distribution too; in addition, the divergence between A and the representative must not be so deep that the
  ortholog Ks peak is saturated and therefore unreliable.

### Indicator 1: median bin height (all clades)

In the ksrates paranome analysis every retained paralog pair carries a gene family and a duplication node. Several
pairs can descend from the same node, so each pair is weighted by `1 / (number of pairs at its node)`, which makes
every duplication node contribute a total weight of 1. This is the same weighting ksrates uses to draw the paranome
histogram, so the sum of weights in a Ks interval is the height of the corresponding bar in the plot.

**Definition.** Weighted duplication nodes are binned in intervals of 0.1 Ks (the `bin_width_paralogs` of the
ksrates configs). The indicator is the **median of the twenty bar heights spanning Ks 0.5 to 2.5**. It is
computed once per species from the paralog export of `ksrates inspect-ks-db`, and is independent of any branch.

**Why the median of a wide window, and not the peak height.** A WGM peak occupies three to five of the twenty
bins and therefore cannot move the median. This matters because most branches of the tree do not host a WGM:
an indicator sensitive to the peak would preferentially select species that already show the signal the analysis
is testing for. Recomputing the median after removing, for every species, the five-bin block centred on its own
tallest bin changes the ranking of the 1178 species by a Spearman rho of 0.986 (typical change in the value
itself, 7%). Species with and without recent WGMs are not separated by this residual (Vitis, no lineage-specific
WGM, 0.78; Glycine soja and Populus, both with recent WGMs, 0.98 and 0.99). Recent WGM peaks lie below Ks 0.5
and fall outside the window by construction.

**It measures data, not biology.** Amborella trichopoda is in the dataset twice, as the 1KP transcriptome URDJ
and as the genome Ambtr_v1.0.27: same species and same WGM history, median bin height 38 and 80 respectively.

**Thresholds.** Below 20 the distribution was judged unreadable, 20 to 30 is a grey zone, above 30 it is usable
(see the validation below). The threshold is only used to flag a branch whose best available representative falls
below it; it never prevents a representative from being chosen.

### Indicator 2: expected ortholog Ks (clades B, C, D)

The divergence between the A representatives and a candidate has to be estimated before any ortholog Ks has been
computed, because the subtrees are defined first and the ortholog pairs to run are derived from them. It is
predicted from the tree:

```text
expected Ks = 1.85 x median over the A representatives of the path distance(A species, candidate)
```

The path distance is the sum of the branch lengths between the two tips. **Those branch lengths must be in
substitutions per site.** The ASTRAL species tree used to cut the subtrees cannot supply them: ASTRAL is given
gene tree topologies and never sees the alignment, so its internal branches are in coalescent units and its
terminal branches carry a placeholder value (1 in the 1KP tree). The branch lengths are therefore taken from the
maximum-likelihood tree of the concatenated amino-acid supermatrix (`tree concat` of
`SupFigs1-3-FigTree.tre`), whose 1178 tip labels match those of the ASTRAL tree exactly. The distance is a
property of the sequence data of the two species, not of the topology the analysis runs on, so the two trees do
not have to agree.

Within one subclade all candidates share the same divergence node with A, so the path distance varies only with
the candidate's own substitution rate; the indicator therefore selects against fast-evolving representatives, and
does so more strongly the deeper the branch is.

**Calibration of the factor 1.85.** The ortholog Ks peaks stored in the two Ks databases of the project (the 1KP
database and the master_1_2026 database, exported with `ksrates inspect-ks-db --type ortholog`) give 110 unique
species pairs whose peak position can be paired with their path distance in the concatenation tree. Restricted to
the 76 pairs with a peak below Ks 1, where the relation is linear, the ratio Ks / path distance has a median of
1.82 and a least-squares slope through the origin of 1.87; 1.85 is used. Including pairs up to Ks 2 raises the
slope to 2.19, so the fit is deliberately restricted to the unsaturated range. The ratio of individual pairs
spans 1.22 to 2.78, so the prediction is adequate for separating saturated from unsaturated divergences and not
for placing a peak.

**The saturation cut-off.** A candidate is treated as saturated once its expected Ks passes a threshold beyond
which the ortholog peak, which positions every divergence of the dataset, can no longer be estimated reliably. That
threshold is set at **Ks 2**, and it is the value the selection rule below uses as its default budget. The reason is
the bootstrap error of the ortholog peak, which over the same 110 pairs is flat up to Ks 2 and degrades beyond it:

```text
Ks band      pairs   median peak SD   SD / Ks
0 - 0.5         54        0.006        1.7 %
0.5 - 1.0       22        0.007        1.1 %
1.0 - 1.5       14        0.012        0.9 %
1.5 - 2.0        5        0.019        1.0 %
2.0 - 3.0        2        0.073        2.9 %
above 3.0       13        0.320        6.6 %
```

Only two of the 110 pairs fall between Ks 2 and 3, so the position of the transition is poorly observed; Ks 2 is
the highest value the data support and is used for that reason.

### Selection rule

Phylogenetic balance allocates the slots exactly as described above, and within each subclade the two indicators
are combined into a single criterion with no tiers:

```text
budget = max(2.0, lowest expected Ks achievable in this subclade + 0.5)
chosen = the candidate with the highest median bin height whose expected Ks <= budget
```

For clade A, where no ortholog divergence is involved, the budget does not apply and the candidate with the
highest median bin height is taken.

The margin of 0.5 Ks expresses "do not pay more than 0.5 Ks of extra saturation in order to gain data". Where no
candidate is saturated, the budget is simply Ks 2 and the median bin height decides; where every candidate is
saturated, the budget floats up to just above the best achievable value, so the best-populated of the
near-minimum candidates is taken instead of either extreme. A hard filter was tried first and rejected: ranking
the saturated case by expected Ks alone selected a hornwort with a median bin height of 8 in exchange for 0.05 Ks,
and ranking it by median bin height alone selected a representative at Ks 3.5.

A branch whose chosen representative lies above the Ks 2 cut-off, or whose A representative falls below a median
bin height of 20, is flagged rather than dropped: the dataset is still built, and the flag records that the result
rests on a saturated anchor or on a distribution too thin to read.

### Validation

The median bin height was validated against visual labelling of the unadjusted paranome plots
(`mixed_<code>_unadjusted_paranome.pdf`). The sheets were blind: they carried only the species code, name and
plot path, with the indicator values kept in a separate key file, and the species order was shuffled.

- **Round 1** (67 species, sampled at the extremes of the candidate indicators) proved uninformative for ranking:
  28 of the 31 species labelled unusable were algae or non-plant outgroups, so the sample recovered taxonomy
  rather than data quality. Duplication-node counts separate the labels perfectly only because every red alga,
  glaucophyte and chlorophyte of the dataset falls below 7000 nodes against 20.5 % of angiosperms. No absolute
  threshold on a raw count can therefore be used, and candidates are compared only among the species available in
  the same clade.
- **Round 2** (36 angiosperm transcriptomes, stratified uniformly on median bin height, labelled on two axes,
  `usable` and `peak`) returned 30 usable, 6 borderline and no unusable species, with the borderline rate falling
  monotonically across the six strata (50, 33, 17, 0, 0, 0 %). AUC against the labels 0.90.
- **Round 3** (16 bryophyte and 12 lycophyte transcriptomes, same design) produced the only unusable labels of the
  three rounds, both below a median bin height of 18, with every species above 30 labelled usable. AUC 0.93.

Splitting the label into `usable` and `peak` was necessary: Hedera helix has 16049 duplication nodes and a median
bin height of 95 but no interpretable peak, which is a biological outcome and not a defect of the data. Under a
single label it would have counted as a failure of the indicator.

Indicators tested on the same labels and rejected, as AUC (angiosperms / bryophytes+lycophytes):

```text
median bin height                0.90 / 0.93    retained
peak bin height                  0.86 / 0.92    rejected: depends on the WGM being present
duplication nodes                0.89 / 0.85    rejected: confounded with lineage and assembly size
BUSCO                            0.73 / 0.83    retained as a tie-breaker only
weighted share of Ks below 0.1    0.66 / 0.51   rejected, and higher values go with usable plots
median alignment length          0.43 /   -     rejected
weighted top family share        0.51 / 0.22    rejected
median alignment coverage        0.41 /   -     rejected
TransRate                        0.24 /   -     rejected: predicts the opposite of usability
```

BUSCO is retained only to break ties. It fails on cases the median bin height resolves (Amborella trichopoda at
93.7 and Annona muricata at 94.1 were both labelled borderline) and it saturates above 91 for every genome of the
dataset, so it cannot rank genome species at all (Klebsormidium flaccidum has a BUSCO of 100.0 and an unusable
paranome).

### Effect on the dataset table and open points

Applied to the 2354 branches with the A representatives held as they are, the rule changes the B, C and D
representatives of 1326 branches (56 %), and 372 branches (16 %) end up with at least one representative above the
Ks 2 cut-off. The expected A|B divergence of the branches is below Ks 1 for 82.3 %, between 1 and 2 for 9.1 %,
between 2 and 3 for 5.1 % and above 3 for 3.5 %, so the saturation term is inactive for most of the tree and
decisive for a minority of deep branches.

To state when publishing:

- The cut-off of Ks 2 and the margin of 0.5 Ks are reasoned choices, not fitted values.
- The conversion factor is calibrated on 110 ortholog pairs from this project only, and its per-pair spread
  (1.22 to 2.78) is large compared with the decision it informs near the cut-off.
- It has not been demonstrated that a less saturated representative yields a better rate adjustment. The
  experiment designed to test it (`ks_analysis/1kp/rate_test/`, six seed-plant sisters of equal age matched on
  BUSCO and median bin height, rate-adjusted against the same focal and outgroup species) was inconclusive: all
  six comparisons fell at Ks 3.5 to 4.8, where the estimates are saturated, and the slow and fast groups were
  confounded with angiosperms and conifers respectively. The adjusted modes of the two groups differed by 0.67
  (3.62 against 4.29) with bootstrap errors of 0.35 to 0.46, and `Ks_Focal`, which should be identical in all six
  trios, varied from 1.78 to 2.23.
- The median bin height describes the paranome. For clade D, whose species is never a focal species, it is used
  as a proxy for the gene content available to form ortholog pairs rather than as a direct measure of the
  ortholog distribution.
