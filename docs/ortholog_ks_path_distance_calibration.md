# Calibrating ortholog Ks against tree path distance, and selecting representative species

Methods record for the branch-centred WGD/WGM survey of the 1KP phylogeny (1178 species,
2354 branches). It documents three things: how the expected ortholog Ks of an unobserved
species pair is predicted from the concatenation tree, what the distribution of those
predictions looks like inside a sister clade, and how those two results constrain the rule
that picks representative species for the older clades B, C and D.

It **supersedes** the calibration described under *Indicator 2: expected ortholog Ks* in
[`make_branch_subtrees.md`](make_branch_subtrees.md), which used a single linear conversion
factor of 1.85 and a hard Ks 2 cut-off. The definition of the A-species indicator
(*Indicator 1: median bin height*) in that document still stands; section 2 below adds one
measurement to it. All numbers in sections 3 and 4 were recomputed on 2026-10-30, after 86 new
ortholog distributions were added in two runs (`calibration_run/`, `calibration_run_deep/`),
which took the collection from 118 to **204 pairs** and changed several of the earlier
conclusions. Where a number has moved, the superseded value is named.

The largest of those changes is that **not every one of those pairs carries a measurement.**
Past a path distance of about 1.0 the ortholog Ks distribution has no peak left to find, so the
calibration is now fitted on the **124 pairs that do** and the rest are reported as bounds. The
deep run described in section 6.1 is what established this; section 3.5 is the result.

Data, figures and scripts: [`ks_analysis/1kp/calibration/`](../../../ks_analysis/1kp/calibration/).
The follow-up run that will extend the calibration:
[`ks_analysis/1kp/calibration_run/`](../../../ks_analysis/1kp/calibration_run/).

---

## 1. Why a calibration is needed

Each branch of the phylogeny is analysed in a small subtree: clade A holds up to four
descendants of the test branch, and clades B, C and D are progressively older sister clades
contributing up to two species each. Choosing the B, C and D representatives requires knowing,
*before any ortholog Ks is computed*, how saturated each candidate's divergence from the A
species will be. Running every candidate pair to find out is not affordable: the 2354 branches
draw on 7054 B/C/D slots.

The predictor used is the **path distance on the concatenation ML tree** - the sum of branch
lengths along the tip-to-tip route through the MRCA, in amino-acid substitutions per site,
taken from `tree concat` in
`source_data/3.phylogenetic_tree/1kp_trees/SupFigs1-3-FigTree.tre`.

The ASTRAL species tree used elsewhere in the pipeline cannot supply this. Its internal branch
lengths are in coalescent units and its terminal branches are placeholder `:1` values, because
ASTRAL never sees the alignment. The two trees carry identical tip labels for all 1178 species,
so the concatenation tree is used purely as a per-species lookup joined to the analysis tree.

---

## 2. The A-species indicator: median bin height

Fully described in [`make_branch_subtrees.md`](make_branch_subtrees.md). One question raised in
review is answered here: **why the median of the twenty 0.1-Ks bars over Ks 0.5-2.5, rather
than the total number of Ks values in that window?**

Measured answer: **the two are equivalent on every test available.**

| statistic | AUC, angiosperm set | AUC, bryophyte+lycophyte set | AUC pooled |
|---|---|---|---|
| median bin height | 0.894 | 0.925 | 0.897 |
| window total (= 20 x mean) | 0.906 | 0.919 | **0.900** |
| mean bar height | 0.906 | 0.919 | 0.900 |
| tallest bar | 0.878 | 0.912 | 0.876 |
| lowest bar | 0.883 | 0.912 | 0.894 |

AUC is the Mann-Whitney statistic against blind visual labels of the paranome plots
(`usable = yes` versus everything else: no, doubtful, borderline). The labelled set is
**64 species** - 36 angiosperms (30 usable) and 28 bryophytes and lycophytes (20 usable) -
scored blind from the plots alone, in `ks_analysis/1kp/ks_quality_review_sheet_*.tsv`.
Reproduce with `indicator_auc.py`. Across all 1178 species the two statistics correlate
at Spearman 0.93, and the window total is marginally *better* at predicting the probability of
a thin bar somewhere in the window (Spearman with the fraction of bars below 30: -0.70 for the
total, -0.68 for the median).

They differ only in the direction of their failure. The total flatters a species carrying one
very tall recent WGD peak (*Hakea drupacea* SIIK: 93rd percentile on the total, 40th on the
median, tallest bar 556 against a median of 56). The median flatters a flat, deep distribution
with no dominant peak (the ferns *Athyrium spicant* VITX, *Polypodium* spp.: 88th percentile on
the median, 62nd on the total).

**The median is retained for interpretability, not for performance.** It is a bar height, in
the same units as what is read off the plot and as the thresholds that were measured on bar
heights (below 20 unreadable, 20-30 grey zone, above 30 usable). The window total is twenty
times a bar height and has no plot-readable meaning. Replacing it with the total would be
defensible provided the thresholds are rescaled accordingly (roughly below 400 / above 600).

---

## 3. The calibration

### 3.1 Data

193 unique species pairs with an ortholog Ks peak on record, exported with
`ksrates inspect-ks-db <address file> --type ortholog` from the two Ks databases of the
project: the 1KP database (120 pairs) and the `master_1_2026` database (73 pairs). 75 of the
1KP pairs were run specifically to populate the deep end of this calibration (section 6). For
each pair, the stored value is the **mode** of the ortholog Ks distribution with its bootstrap
standard deviation (200 replicates), and the predictor is the concatenation-tree path distance
between the same two species.

**124 of the 193 are used for fitting**: those with a path distance of at most 1.0, which are
the ones whose Ks distribution actually has a mode. The other 69 are retained in the table and
plotted, but their recorded peaks are artefacts of the `max_ks_orthologs` cutoff rather than
measurements (section 3.5). The 11 pairs of the deep run (section 6.1) are not in this count:
they lie between path 2.38 and 3.42 and none of them yielded a peak at all.

Table: [`calibration_pairs.tsv`](../../../ks_analysis/1kp/calibration/calibration_pairs.tsv).
Figure: `ortholog_ks_vs_path_distance.png`.

### 3.2 The relationship is not linear

The local conversion factor `Ks / path` rises steadily with divergence:

| ortholog Ks band | n pairs | median Ks/path | IQR | CV | min - max |
|---|---|---|---|---|---|
| 0 - 1 | 81 | 1.82 | 1.69 - 1.97 | 13.3% | 1.22 - 2.78 |
| 1 - 2 | 35 | 2.41 | 2.30 - 2.65 | 9.4% | 1.89 - 2.81 |
| 2 - 3 | 12 | 2.43 | 2.36 - 2.58 | 5.7% | 2.29 - 2.71 |
| > 3 | 65 | 2.70 | 2.45 - 3.07 | 16.2% | 1.96 - 3.94 |

The drift is already present inside the densely sampled Ks 0-2 range (1.82 to 2.41), so it is
not produced by the deep pairs alone. **Two corrections to the 118-pair version of this table.**
The factor climbs to about 2.4 by Ks 1-2 and then roughly plateaus, rather than continuing to
3.01: that value came from a deep arm made up entirely of *Isoetes* and *Physcomitrella* pairs.
And the deep band is far noisier than it appeared - CV 16.2% against the 7.6% measured on those
13 pairs - so the apparent precision at depth was an artefact of a single lineage combination. The direction is that **Ks grows faster than protein
distance** at large divergence. Two non-exclusive causes are consistent with this: the
multiple-hit correction applied to Ks becomes unstable as synonymous sites saturate and tends
to inflate the estimate; and the concatenation supermatrix is assembled from conserved genes,
whose amino-acid distances accumulate comparatively slowly.

### 3.3 Fitted models

**Which pairs the models are fitted on.** Only the 124 pairs with a path distance of at most
1.0. Beyond that the ortholog Ks distribution no longer has a peak to measure: it flattens into
a smear against the `max_ks_orthologs` cutoff, and the mode the kernel density estimate returns
is a property of the cutoff rather than of the species pair. Section 3.5 gives the measurement.
An earlier version of this calibration was fitted on all 193 pairs and therefore fitted partly
to those artefacts. The 69 excluded pairs are still drawn in the figure, unfilled, because
showing where the measurement fails is part of the result.

With `path` as the response and ortholog `Ks` as the predictor, over those 124 pairs. That is the
direction in which the pairs were *measured* - the bootstrap error belongs to the Ks peak - but it
is the opposite of the direction in which the calibration is *used*, so the table gives both forms
and the figure is drawn in the prediction direction (path on x, Ks on y). The same curve read the
two ways looks different: `path(Ks)` has an exponent below 1 and flattens, while its inverse
`Ks(path)` has an exponent above 1 and accelerates.

Fitting the regression the other way round instead - `log(Ks)` on `log(path)`, which is formally
the right one for prediction - gives `Ks = 2.59 * path^1.222`. Its overall median error is slightly
lower (7.0% against 7.5%) but it is worse in every individual Ks band (6.9 / 9.0 / 7.0 against
6.6 / 8.5 / 5.2), so the overall figure is a mixing artefact and the fit below is retained.

| model | fit | inverse, used for prediction |
|---|---|---|
| linear through origin, fitted on Ks <= 1 (n = 81) | `path = Ks / 1.92` | `Ks = 1.92 * path` |
| linear through origin, fitted on Ks <= 2 (n = 116) | `path = Ks / 2.28` | `Ks = 2.28 * path` |
| **power law (adopted)** | `path = 0.4531 * Ks^0.790` | `Ks = (path / 0.4531)^1.266` |
| saturating in Ks | `path = 1.352 * (1 - e^(-0.426 Ks))` | `Ks = -ln(1 - path / 1.352) / 0.426` |
| saturating in path | degenerate, see below | - |

Error of the predicted Ks, as median absolute error / standard deviation, by observed-Ks band:

| model | Ks 0-1 (n=81) | Ks 1-2 (n=35) | Ks 2-3.3 (n=8) | all |
|---|---|---|---|---|
| linear 1.92 | 0.03 / 0.07 | 0.29 / 0.11 | 0.58 / 0.33 | 0.07 / 0.23 |
| linear 2.28 | 0.09 / 0.07 | 0.10 / 0.12 | 0.28 / 0.33 | 0.10 / 0.17 |
| **power law** | **0.03 / 0.08** | **0.12 / 0.18** | **0.11 / 0.32** | **0.05 / 0.14** |
| saturating in Ks | 0.03 / 0.08 | 0.14 / 0.21 | 0.25 / 0.33 | 0.04 / 0.16 |
| saturating in path | 0.11 / 0.08 | 0.07 / 0.13 | 0.23 / 0.33 | 0.10 / 0.17 |

**A constant factor is the wrong shape, and it errs in the dangerous direction.** Because the
relationship is concave, any factor fitted where the data are dense is too small further out,
so it systematically *underestimates* the divergence of the deeper candidates - precisely the
error that lets an unmeasurable species through a budget. Signed median error (negative =
underestimate), and the share of pairs underestimated:

| model | signed err, Ks 0-1 | Ks 1-2 | Ks 2-3.3 | % under, Ks 0-1 | Ks 1-2 | Ks 2-3.3 |
|---|---|---|---|---|---|---|
| linear 1.92 | +0.02 | -0.29 | -0.58 | 33% | 97% | **100%** |
| linear 2.28 | +0.09 | -0.08 | -0.28 | 1% | 77% | **100%** |
| **power law** | **-0.00** | **-0.01** | **+0.04** | **51%** | **51%** | **38%** |
| saturating in Ks | +0.01 | -0.06 | +0.19 | 32% | 66% | 38% |
| saturating in path | +0.11 | -0.05 | -0.23 | 1% | 63% | 88% |

**The adopted calibration is the power law**

```
Ks = (path / 0.4531) ^ 1.266
```

It is monotone, has the right curvature, is unbiased in every band, and is the best of the
candidates **exactly where the survey works**. Almost all WGM inference happens at low Ks, both
because 88% of the B/C/D slots sit below a path of 1.0 and because those are the only Ks values
that are measurable at all. Over the 116 pairs with Ks <= 2:

| model | signed median error | median absolute error | as % of Ks | SD of log(obs/pred) |
|---|---|---|---|---|
| **power law (adopted)** | **-0.001** | **0.044** | **7.7%** | **0.142** |
| saturating in Ks | +0.012 | 0.039 | 7.9% | 0.145 |
| linear 1.92 | -0.006 | 0.052 | 11.9% | 0.181 |
| linear 2.28 | +0.072 | 0.095 | 18.0% | 0.181 |
| saturating in path | +0.091 | 0.103 | 20.0% | 0.181 |

Restricting further to the 81 pairs with Ks <= 1: 6.6% for the power law, against 7.6%, 8.9%,
25.0% and 28.3%. The saturating-in-Ks form is its equal on this range; the power law is preferred
because it extrapolates monotonically and has one fewer awkward property - its two-parameter
rival has an asymptote in `path` at 1.352, barely above the fitted range, so a path distance
only slightly larger maps to infinite Ks.

**Saturating forms are not identifiable on these data, and that is expected.** `Ks = Kmax (1 -
e^(-k path))` is the biologically motivated shape, since synonymous sites really do saturate.
Fitted to the trusted pairs it degenerates to a straight line (`Kmax = 2.2e4, k = 4e-5`), which
is the correct answer: saturation is not visible in a range where Ks has not saturated. Fitting
it to the full 193 pairs did return a bend, but that bend was being driven by the artefactual
deep peaks, and its ceiling was unidentified in any case (`Kmax = 17.1 +/- 4.0`). **There is
therefore no ceiling parameter in this calibration, and the earlier `KS_CEILING = 5.1` has been
removed** - see 3.5 for why that number was never a biological ceiling.

**A conservative variant for the flag.** The power law is centred, so about half of all
candidates are still underestimated. Over the fitted pairs the ratio of observed to predicted Ks
has p90 = 1.19, so multiplying the prediction by **1.19** leaves only ~10% underestimated
(93% / 83% / 88% covered, by band). This margin is for the flag that marks a candidate as
likely unmeasurable; the reported expected Ks stays unmultiplied.

**Limit of validity.** The prediction is quoted as a number only up to a path of 1.0, which the
power law maps to Ks 2.72. Past that `report_ks()` returns a bound and a reliability level
rather than a value, while `path_to_ks()` keeps returning the raw number because the selection
rule needs an *ordering* and path distance stays monotone in divergence even where the Ks peak
has saturated. The thresholds are stated on path, not on Ks, because what degrades is the
measurability of the peak rather than its value.

### 3.4 Is there a better predictor than protein path distance?

Tested, on the raw codeml output of all 118 calibration pairs (every `*.ks.tsv` produced by
the ortholog and paranome runs carries `Ka`, `Ks` and `Omega` per gene pair). **No cheaper
ruler available in these data beats the concatenation path, and two plausible corrections make
the calibration worse.**

**A Ka ruler from the same codeml run is no better.** For each pair the median Ka over its
one-to-one orthologs was compared with the concatenation path as a predictor of the stored
ortholog Ks mode. Coefficient of variation of the ratio, within band:

| ortholog Ks band | n | CV of Ks / concat path | CV of Ks / Ka (same run) |
|---|---|---|---|
| 0 - 1 | 81 | 13.3% | 14.5% |
| 1 - 2 | 35 | 9.4% | 11.7% |
| 2 - 3 | 12 | 5.7% | 9.4% |
| > 3 | 65 | 16.2% | **10.0%** |

The Ka ruler shares its genes, its alignments and its codeml run with the Ks it is predicting,
so its errors are correlated with the target and it should be flattered by the comparison.
Below Ks 3 it is not: the 33%-occupancy amino-acid supermatrix tracks ortholog Ks as well as
the very same alignments do, which is a strong result for an independent gene set under a
single empirical AA model.

**Above Ks 3 the picture reverses, and this is a correction to the 118-pair version.** With 65
deep pairs instead of 13, the Ka ruler is clearly the better of the two (CV 10.0% against
16.2%), so at that depth the supermatrix path *is* the weaker instrument. This does not give a
usable predictor - Ka for a pair is only known once that pair has been run - but it does locate
where the remaining error lives, and it means the deep end of the calibration is limited by the
supermatrix rather than by biology.

**Lineage omega does not transfer.** Omega is the ratio of the nonsynonymous to the synonymous
substitution rate, dN/dS = Ka/Ks; values below 1 indicate purifying selection. Taking a
per-species omega from each species' own paranome run (median Ka/Ks over paralog pairs with Ks
between 0.1 and 2, available for all 1178 species) and dividing the path by the pair mean makes
every band worse: 13.3% -> 16.1%, 9.4% -> 19.9%, 5.7% -> 27.9%, 16.2% -> 28.8%. The paranome
omega is 2.40x the pair's ortholog omega - duplicated genes are under relaxed constraint - and
on the 193-pair set it is not merely uninformative but mildly *anti*-correlated with it
(r = -0.22; it was +0.18 on 118 pairs, so the apparent weak signal did not survive more data).

**The pair's own omega does explain part of the residual, but circularly and inconsistently.**
Ks/path correlates with 1/omega at r = +0.62 (Ks 0-1) and +0.74 (Ks > 3), and dividing the path
by the pair's own omega tightens those bands (13.3% -> 12.3%, 7.6% -> 5.2%). But in the Ks 1-2
band the correlation is r = +0.09 and the correction makes things worse (10.2% -> 15.8%); and
the omega of a pair is only known once that pair has been run, which is what the prediction is
meant to avoid. Moreover omega is largely a restatement of depth rather than an independent
lineage property: over the 193 pairs r(omega, path) = -0.86, because as synonymous sites
saturate the inflated Ks drives Ka/Ks down. A lineage-pair grouping appears to explain 63% of
the variance in omega, but that grouping is itself confounded with depth
(angiosperm/angiosperm 0.163, angiosperm/gymnosperm 0.066).

**Conclusion.** Below Ks 3 the residual 5-13% scatter is genuine variation in the dN/dS of the
gene sets being compared, and the only instrument that measures it for a given pair is the
ortholog run itself; the predictor is close to the ceiling of what these data support. Above
Ks 3 the supermatrix path is the limiting factor rather than biology, which is where any future
improvement would have to come from. A further
consideration argues for keeping a protein ruler rather than a synonymous one: protein distance
does not saturate over the range of interest (0.062 to 1.73 substitutions per site, across pairs
spanning ortholog Ks 0.11 to 5.03), so it still carries information at Ks 2-4 - exactly where
the selection decision is made, and where any synonymous-based predictor has stopped being
estimable.

### 3.5 Where the ortholog Ks peak stops being a measurement

This is the central limitation, and it was found by running pairs deeper than anything the
calibration had covered (section 6.1) and then inspecting the *distributions* rather than the
peaks.

**The deep distributions have no peak.** For all eleven pairs of that run, between 85% and 93%
of ortholog Ks values exceed 10 and fall outside the window entirely, and within 0-10 the
histogram is flat: `GUBD_Nelnu`, at a path of 3.42, puts between 5% and 9% of its mass in every
0.5-wide bin from Ks 3 to Ks 10. The peak reported for such a distribution is wherever the
kernel density estimate happened to land. Across the eleven it ranges from 4.12 to 8.75 and is
uncorrelated with path distance. The twelfth pair, `BZSH_Poptr` at path 3.51, completed later
and behaves identically (concentration 18.3%, nominal peak 8.07).

A further point follows from these twelve being added to the database. The highest ortholog Ks
on record was 5.10 when the former ceiling was adopted; it is now **8.63**, and ten pairs exceed
the old ceiling. Nothing about the data changed - the value was simply the maximum of whatever
sample happened to have been run, which is what a sample maximum does when the sample grows. Had
it been retained as a constant, the database would now hold ten values it declares impossible.

**The failure begins far earlier than that run.** Scoring every ortholog Ks distribution
available on disk by its *concentration* - the fraction of retained Ks values within +-0.5 of
the mode, which a flat smear over the 0-10 window would give as 10%:

| path | 0-0.6 | 0.6-0.8 | 0.8-1.0 | 1.0-1.2 | 1.2-1.5 | 1.5-1.8 | 1.8-2.1 | 2.1-2.4 | 2.4-3.6 |
|---|---|---|---|---|---|---|---|---|---|
| n | 7 | 6 | 7 | 8 | 14 | 9 | 14 | 11 | 10 |
| % of Ks > 10 | 2.9 | 7.8 | 17.6 | 38.0 | 45.9 | 57.0 | 69.8 | 87.1 | 86.8 |
| median mode | 1.00 | 1.73 | 2.22 | 3.56 | 3.87 | 4.97 | 5.95 | 6.84 | 5.10 |
| **concentration** | **68%** | **55%** | **41%** | **21%** | **20%** | **17%** | **16%** | **16%** | **16%** |

A real mode survives to about path 1.0 and is gone by 1.5. These two values are the thresholds
`RELIABLE_PATH` and `BORDERLINE_PATH` in `calib_common.py`.

**Consequences.**

1. **The former `KS_CEILING = 5.1` was not a biological ceiling.** It was the largest value the
   kernel density estimate happened to return from flat smears bounded by `max_ks_orthologs =
   10`. The earlier text in this document read the plateau of median observed Ks (4.63, 4.78,
   5.01, 4.86 across successive deep path bins) as evidence that synonymous saturation imposes a
   real limit near 4.9. That reading was wrong: the plateau is what a bounded uniform
   distribution produces, and it moves with the cutoff rather than with biology. The constant
   has been removed.
2. **Sixty-nine of the 193 calibration pairs were above path 1.0, and 39 above 1.5.** Their Ks
   values, up to 5.10, are not measurements, and the earlier power law was fitted partly on
   them. Refitting on the 124 trusted pairs changes the exponent from 1.204 to 1.266. Over the
   trusted pairs themselves the two fits are almost indistinguishable (median absolute error
   7.5% against 7.6%), so the contaminated arm did little damage to the shallow range; the
   difference is entirely in the extrapolation, which is what is no longer quoted.
3. **A three-parameter model with a real ceiling is not fittable**, and should not be, because
   there is no measured ceiling to fit.

**An independent check: the peak stops being estimator-independent.** The concentration
statistic and the ksrates peak both come from a kernel density estimate, so a sceptical reader
could ask whether the flatness is an artefact of that estimator. It is not. Comparing the peak
ksrates stores (bootstrapped, 200 replicates) against a plain `gaussian_kde` mode computed on
the same Ks list, the two agree closely while a mode exists and then stop agreeing:

| path | 0-0.8 | 0.8-1.0 | 1.0-1.5 | 1.5-2.4 | 2.4-4.0 |
|---|---|---|---|---|---|
| n | 13 | 7 | 22 | 34 | 11 |
| median difference (Ks) | 0.03 | 0.04 | 0.03 | 0.94 | 0.93 |
| as % of the peak | 2% | 2% | 1% | **19%** | **15%** |

A real mode is located to 1-2% by either estimator; a flat distribution is not located at all,
and which value is reported depends on the bandwidth. Note that the agreement survives to path
1.5, beyond the point at which concentration has already fallen: between 1.0 and 1.5 the mode is
broad but still locatable, which is precisely what the borderline band is meant to express - usable
for ordering candidates, not for quoting a number.

**A measured flag, not only a predicted one.** Because the diagnostic is computed from a run's
own output, it can be reported per analysis rather than inferred from path distance in advance.
`ortholog_distribution_quality.py` does this, calling a distribution *measured* at
concentration >= 40%, *degraded* at >= 20%, and *no peak* below. Over the 86 distributions
available the path-based flag and the measured verdict agree on 83 (97%), with the three
disagreements all inside the borderline band, where disagreement is what the band is for. Path
distance serves as the a priori warning; concentration is the a posteriori verdict, and it is
the one to quote in results.

### 3.6 Other limitations

1. **Per-pair scatter is substantial**: 13% CV of the conversion factor at low Ks, with
   individual pairs ranging from 1.22 to 2.78, i.e. about 0.2 Ks standard deviation on any
   single predicted value. No selection threshold finer than that is resolvable.
2. **Ks 2-3 is thin**: 8 trusted pairs, and it is the band in which the selection rule most
   often has to decide.
3. **The two databases are partly confounded with depth**: `master_1_2026` contributes only
   pairs below Ks 2.
4. **The calibration is not transferable to a different `max_ks_orthologs`** without recomputing
   the concentration table, since the uniform floor is `2 * 0.5 / max_ks_orthologs`.

---

## 4. Expected Ks inside a sister clade

For a given test branch, every candidate in clade B (or C, or D) receives an expected ortholog
Ks, computed as the calibration applied to the median, over the A representatives, of the path
distance between that A species and the candidate. The median over A is used rather than the
mean because the A representatives are themselves a small, deliberately heterogeneous sample
(up to two species from each direct child of the test branch) and one unusually fast A species
should not shift the expectation for every candidate.

Figure: `expected_ks_distributions_twelve_branches.png`, twelve branches spanning the depth of
the tree.

| test branch (A) | B candidates | expected Ks range | spread | lowest 25% keeps | min + 0.5 keeps |
|---|---|---|---|---|---|
| Cupressaceae | 8 | 0.32 - 0.37 | 0.05 | 2 | 8 |
| Solanaceae | 12 | 0.78 - 0.99 | 0.21 | 3 | 12 |
| Pinaceae | 59 | 0.76 - 1.09 | 0.34 | 15 | 59 |
| Lamiales | 57 | 0.89 - 1.30 | 0.41 | 15 | 57 |
| Caryophyllales | 219 | 1.20 - 2.74 | 1.54 | 55 | 116 |
| Ranunculales | 528 | 0.82 - 2.65 | 1.82 | 132 | 377 |
| Monocots | 577 | 1.38 - 3.41 | 2.02 | 145 | 62 |
| Liverworts | 42 | 2.04 - 3.13 | 1.09 | 11 | 6 |
| Hornworts | 64 | 2.21 - 3.31 | 1.10 | 16 | 13 |
| Gymnosperms | 692 | 1.99 - 4.30 | 2.31 | 173 | 36 |
| Lycophytes | 847 | 2.55 - 6.28 | 3.73 | 212 | 18 |
| Chlorophytes | 995 | 6.71 - 12.07 | 5.36 | 249 | 9 |

Expected Ks is the power law of section 3.3 applied to the median path distance from the A
representatives. The chlorophyte row is past the calibrated range and should be read as "deeper
than anything measured". The *ordering* of candidates does not depend on the calibration, which
is monotone in path distance, so the shapes described below are unaffected by it; only the
values on the axis are.

Replacing the old constant factor with the power law raises the number of B/C/D slots whose
selected representative sits above Ks 2 from 842 (12%) to 1081 (15%) - the correction is modest
in aggregate because most slots are shallow, but it is decisive for the deep branches, which
are exactly the ones a constant factor flattered.

### 4.1 Location behaves as expected, width does not

Recent branches sit near Ks 0.3-0.9 and ancient ones at 2-3, as anticipated. The **spread**,
however, ranges from 0.04 to 2.63 Ks - a factor of 60 between branches. Any rule phrased as a
fraction of the clade therefore means something entirely different from one branch to the next.

### 4.2 The distributions are frequently not unimodal

Across all 2539 B/C/D slots holding at least four candidates: median skew 0.05-0.08, but ~30%
of slots skew right beyond 0.5, ~25% skew left beyond -0.5, and **~24% contain a clear internal
gap** (a jump between consecutive sorted values, inside the middle 80% of the data, larger than
one third of the total range). A Gaussian is not a safe default.

### 4.3 Where a clade is bimodal, the modes are rates, not depths

The worked case is the **lycophyte stem**: A = ZZEI, PQTO, Selmo_v1.0, PYHZ; B = all
euphyllophytes, 847 candidates. Figure:
`expected_ks_distribution_lycophyte_sister.png`; data:
`lycophyte_sister_candidates.tsv`.

| group in B | n | min | median | max | median bin height |
|---|---|---|---|---|---|
| ferns | 73 | 2.55 | 3.74 | 4.43 | 69 |
| gymnosperms | 82 | 2.87 | **3.38** | 4.21 | 80 |
| angiosperms | 692 | 3.78 | 5.25 | 6.28 | 72 |

The split is not ferns against seed plants: **ferns and gymnosperms overlap almost completely
and form the left mode together, while the angiosperms alone form the right one**, with the
gymnosperms the slowest of the three. All three groups diverge from the lycophytes at the same
node, so the two modes are not two divergence times - they are two substitution rates.
Angiosperms are faster, which places their expected ortholog Ks about one unit deeper for
exactly the same divergence.

Group assignments were verified against the NCBI lineage in
`source_data/1.species_dataset/1kp_paper_2019_suptab1_species_ncbi_taxonomy.csv`: 826 of the
847 candidates appear in that table and none disagree. The remaining 21 are published genomes,
absent from the 1KP supplementary table, assigned by hand.

### 4.4 How often is there a choice at all?

Over the 2354 branches there are 7054 B/C/D slots:

| sister-clade size | slots | share |
|---|---|---|
| 1 species | 2832 | 40.1% |
| 2 - 3 | 1683 | 23.9% |
| 4 - 10 | 1277 | 18.1% |
| 11 - 50 | 856 | 12.1% |
| > 50 | 406 | 5.8% |

For 64% of slots there are at most three candidates and for 40% exactly one, so no selection
rule changes the outcome. The rule only does work in the ~18% of slots with more than ten
candidates - which are, by construction, the deep and taxonomically heterogeneous clades.

---

## 5. The selection rule for clades B, C and D

### 5.1 The rule

```
budget  = max(2.0, min(expected Ks over the candidates of the clade) + delta)
pick    = the candidate with the highest median bin height among those with
          expected Ks <= budget
```

with `delta` currently 0.5 Ks. The rule chooses among the species available in the B clade (or
C, or D) - one representative per direct child subclade, as before. Both criteria are applied
to D as well as to B and C: D anchors the rate adjustment and therefore needs a well-populated
distribution for the same reason the others do.

The `min + delta` arm is what makes the rule adaptive: when every candidate is saturated, the
budget opens up around the least-saturated one and the best available representative is still
returned rather than none.

**No candidate is excluded for being too deep, and no branch is skipped.** It would be possible
to refuse any candidate beyond `RELIABLE_PATH`, since an ortholog Ks peak that cannot be
measured cannot anchor a rate adjustment. That is rejected here for two reasons. The first is
practical: this survey exists to re-examine the WGM calls of the 1KP paper, which include calls
on the deepest branches of the tree, and a method that declines to analyse those branches cannot
check them. The second is that it would be a misdirected refusal - path distance is monotone in
divergence whether or not Ks has saturated, so the *ordering* of candidates, which is all the
rule needs, remains valid throughout. What becomes invalid is the quoted number, not the choice.

Every branch is therefore analysed, and the output carries a reliability level
(`calib_common.path_reliability`): **reliable** at path <= 1.0, **borderline** to 1.5, **unreliable**
beyond. After the run the level is confirmed or corrected against the distribution actually
obtained (3.5). The flag belongs to the *rate adjustment*, not to the WGM call: the paranome Ks
distribution of the A species is unaffected by how distant the outgroup is. What degrades is the
positioning of the ortholog divergence line that assigns a duplication to one side of the test
branch or the other, so the correct caveat on a flagged branch is **"branch assignment
uncertain"**, not "WGM call invalid".

**What the limit would have cost, had it been applied.** Over all 2354 branches and 7054 B/C/D
slots, the fraction keeping at least one candidate within a given path limit:

| path limit | ~Ks | slots keeping a candidate | branches keeping one anywhere |
|---|---|---|---|
| 0.8 | 1.92 | 84.4% | 87.9% |
| 1.0 | 2.72 | 87.6% | 90.7% |
| 1.5 | 4.55 | 93.9% | 96.0% |

The cost is small because the closest available outgroup has a median path of 0.18 (p90 0.91):
the 1KP sampling is dense. It is also very unevenly distributed. Branches whose A species are
angiosperms never lose a candidate at a limit of 1.0, and gymnosperms, ferns, lycophytes and
bryophytes lose 0.6-4.5%; of the 218 branches that would be lost, **210 are algal**. The
land-plant branches that would be caught are identifiable by the size of the clade opposite
them - the bryophyte/vascular split (path 1.89), the lycophyte/euphyllophyte split (1.37), the
fern/seed-plant split (1.48), the gymnosperm/angiosperm split (1.25) - that is, the deep
backbone nodes carrying the 1KP paper's headline WGD claims. These are precisely the branches
that must be analysed and flagged rather than dropped.

### 5.2 A quantile budget was evaluated and rejected

The alternative considered was to replace the absolute Ks 2 term with a quantile - build the
distribution of expected Ks within the clade and keep the lowest 25%. It fails in both
directions, for two different reasons.

**In a wide clade it is too permissive.** For the lycophyte stem the lowest quarter of the 847
candidates still reaches Ks 2.99:

| rule | window | candidates kept | representative selected |
|---|---|---|---|
| lowest 25% of the clade | 1.88 - 2.99 | 212 | *Nelumbo nucifera*, expected Ks 2.89, bin height 131 |
| `max(2.0, min + 0.5)` | 1.88 - 2.38 | 61 | *Pinus ponderosa* JBND, expected Ks 2.25, bin height 121 |

The quantile buys 10 points of bin height for 0.64 Ks of additional saturation.

**In a narrow clade it discards candidates for a difference it can barely resolve.** The
calibration's error is multiplicative: in Ks units its standard deviation runs from 0.08 in the
Ks 0-1 band to 0.93 above Ks 3, a twelvefold range, while the standard deviation of
log(observed/predicted) stays near 0.16 in every band. Of that log-variance about 40% is a
property of one species and is therefore shared by all of its pairs; within a branch the A side
is fixed, so that share cancels and what remains to scramble a ranking of candidates is
**0.125 log units, about 13% of the predicted value**.

Against that, the candidates of one clade are tightly packed: over the 2539 B/C/D slots with at
least four candidates, the median spread between the least and most diverged candidate is
0.33 Ks, which is 33% of the slot's own level - a signal-to-noise ratio of about 2.6. In 16% of
slots the spread is smaller than the error outright.

Simulating that (the true expected Ks of the candidates, multiplicative noise of 0.125 log
units, 20000 replicates, selecting the candidate with the lowest predicted value):

| clade | n | level | spread | P(pick is truly the lowest) | mean true percentile of the pick |
|---|---|---|---|---|---|
| typical B/C/D slot | 10 | 0.91 | 0.30 | 0.39 | 16% |
| lycophyte clade after `min+0.5` | 18 | 2.86 | 0.46 | 0.15 | 26% |
| whole lycophyte clade | 847 | 5.06 | 3.73 | 0.14 | 1% |

(random choice would put the pick at the 50th percentile)

**This is a correction to an earlier version of this document**, which used a single additive
error of 0.17 Ks - a figure taken from the shallow bands of the superseded 118-pair calibration -
and concluded that a within-clade ranking was "barely better than chance". It is not. In a
typical slot the selected candidate lands at the 16th percentile of the true values, which is
clearly informative, and the deeper the clade the better the ordering separates its extremes
(1st percentile for the lycophyte clade as a whole). What degrades is the ability to identify
the single lowest candidate, and that matters less than the ordering.

So the case against a quantile budget does **not** rest on the ranking being noise. It rests on
the two structural points above: a quantile is a relative criterion with no knowledge of which
absolute Ks is dangerous, so in a wide clade its lowest quarter still reaches deep into
saturation, and in a narrow clade it discards three quarters of the candidates to gain a
fraction of a Ks that is of the same order as the calibration's own error.

**Third failure mode, visible in the twelve panels.** When the whole clade is unsaturated
(Cupressaceae through Ranunculales in the table of section 4) `min + 0.5` keeps every
candidate - which is the correct behaviour, since the choice should then rest purely on data
availability - while the quantile discards three quarters of them for no gain.

### 5.3 Recommended form

Keep an **absolute tolerance**, `budget = min_expKs + delta`, with `delta` set by the
calibration's own error rather than by a quantile. Candidates within `delta` of the best are
indistinguishable in saturation and should be separated on median bin height; candidates beyond
it are genuinely more saturated. The present `delta = 0.5` is roughly two standard deviations
of the prediction error and is defensible; 0.25 would be one.

The Ks 2 term is then best demoted from a gate to a **confidence flag** on the branch output,
raised when the selected representative's expected Ks exceeds 2, rather than a threshold that
can silently return nothing. This is the same demotion applied to the path limit in 5.1, and for
the same reason: a survey that must reproduce deep WGM calls cannot afford a rule that answers
"no representative" on the branches of most interest.

---

## 6. Planned extension of the calibration

[`ks_analysis/1kp/calibration_run/`](../../../ks_analysis/1kp/calibration_run/) - 76 new
ortholog Ks distributions over 59 species. Pair table:
[`recommended_ortholog_pairs.tsv`](../../../ks_analysis/1kp/calibration/recommended_ortholog_pairs.tsv).

Design:

1. **Pool.** The tree is partitioned into 15 lineage buckets (chlorophyte and streptophyte
   algae; liverwort, moss, hornwort; lycophyte; fern; conifer and other gymnosperm; basal
   angiosperm; monocot; basal, caryophyllid, rosid and asterid eudicot). Each contributes its
   four best species ranked by **published genome first, then BUSCO completeness** - not by
   paranome depth, because an ortholog Ks estimate needs a complete gene set rather than many
   paralogs. 20 of the 59 selected species are published genomes.
2. **Sampling is uniform in path distance, not in Ks**, since path distance is the predictor.
   Four pairs in each 0.1-wide band from 0.45 to 2.35 (predicted Ks about 0.8 to 4.3).
3. **Lineage replication.** Each band draws its four pairs from two or three *different* bucket
   combinations; no bucket combination is used more than twice anywhere; no species appears in
   more than five pairs. This is what the present calibration lacks.
4. Pairs already present in `calibration_pairs.tsv` are excluded.

Coverage gained:

| path band | predicted Ks (1.85) | pairs already available | pairs added |
|---|---|---|---|
| 0.45 - 0.95 | 0.83 - 1.76 | 24 | 20 |
| 0.95 - 1.45 | 1.76 - 2.68 | 5 | 20 |
| 1.45 - 1.95 | 2.68 - 3.61 | 9 | 16 |
| 1.95 - 2.35 | 3.61 - 4.35 | 0 | 16 |

Only 14 of the 118 existing pairs had a path distance of 0.95 or more, and 8 of the 19 bands
above that value were empty. That run has since been executed: 75 of its 76 pairs completed
(the *Chlamydomonas* x *Brachiomonas* pair was abandoned, see below) and the calibration now
rests on 193 pairs.

## 6.1 A second, smaller run past the calibrated range

`ks_analysis/1kp/calibration_run_deep/` - **12 pairs at path 2.35 to 3.55**, beyond the deepest
pair the power law was fitted on.

Their purpose is not to extend the fit. Of the 7054 B/C/D slots only 164 (2.3%) select a
representative past path 2.31, and every one of those is already beyond any usable ortholog Ks,
so a better value there changes no decision. They are there to answer a methods question:
**is the ortholog peak still a measurement at that depth, or the mode of a distribution piling
up against saturation?** The plateau and the growing bootstrap error in section 3.3 are what
raise it, and the answer belongs in the publication either way. The full Ks distributions should
be inspected, not just the peaks harvested, so `x_axis_max_limit_orthologs_plots` is 10 there.

**Outcome: the peak is not a measurement at that depth.** Eleven of the twelve pairs completed.
Every one returned a flat distribution - 85-93% of Ks values above 10, and within the window a
histogram indistinguishable from uniform - with a nominal peak ranging from 4.12 to 8.75 and
uncorrelated with path distance. The run did its job: it showed that the question it was asked
to settle has no numeric answer, that the apparent ceiling near Ks 5 was an artefact of the
cutoff, and that the breakdown begins at a path of about 1.0 rather than at 2.3. Section 3.5
gives the measurement and what followed from it. No three-parameter model was fitted, because
there is no ceiling in the data to fit.

*Chlamydomonas reinhardtii* is excluded from the pool. Its proteome is 70.2% GC, which pushes
alanine, glycine, proline and arginine to 42.6% of its residues against 28.1% in its partner,
with 3.19% of residues inside runs of five or more identical amino acids and a mean protein
length of 734 against 280. With blastp's SEG filter off by default this makes the all-vs-all
search pathologically slow: one such pair ran 11.5 hours and produced 4.8 MB of output where a
comparable pair produced 657 MB, and it was abandoned. Past a path of 2.35 one partner is
necessarily an alga, so lineage replication can only vary the other side; it does, across
lycophyte, hornwort, liverwort, moss, fern, conifer and four kinds of angiosperm.

---

## 7. Files and reproduction

**Code** lives in the version-controlled repository at
[`code/1kp_wgms/others/calibration/`](../others/calibration/); `ks_analysis/` is not under
version control and holds only data and figures. Every script is run directly with
`/software/shared/apps/x86_64/python/3.11.4/bin/python3` (stdlib plus numpy, scipy,
matplotlib) and takes no required arguments.

| script | what it does |
|---|---|
| `calib_common.py` | shared module: newick parsing, the two trees, memoised tip-to-tip path distances, species metadata, taxonomic groups, bin heights, subtree construction |
| `paralog_ks_bin_heights.py` | scans all 1178 paranome runs -> `paralog_ks_bin_heights.tsv` (20 bar heights, median, total per species) |
| `indicator_auc.py` | scores each summary of the bar heights against the blind review sheets (section 2) |
| `build_calibration_pairs.py` | joins the newest ortholog peak exports to path distance -> `calibration_pairs.tsv` |
| `plot_calibration.py` | fits the models on the trusted pairs, reports their errors, draws `ortholog_ks_vs_path_distance.*` (section 3) |
| `alternative_predictors.py` | the Ka-ruler and omega tests (section 3.4) |
| `plot_lycophyte_sister.py` | the lycophyte-sister figure, its candidate table and the ranking simulation (sections 4.3, 5.2) |
| `plot_sister_clade_distributions.py` | the twelve-branch figure and the clade-size table (sections 4, 4.4) |
| `ortholog_distribution_quality.py` | scores a run's ortholog Ks distributions by concentration and reports whether each peak is a measurement (section 3.5) |
| `recommend_ortholog_pairs.py` | selects the next set of pairs worth running (section 6) |
| `../root_to_tip_distances.py` | root-to-tip distances from the concatenation tree |

**Data and figures** are written to
[`ks_analysis/1kp/calibration/`](../../../ks_analysis/1kp/calibration/):
`calibration_pairs.tsv`, `paralog_ks_bin_heights.tsv`, `lycophyte_sister_candidates.tsv`,
`recommended_ortholog_pairs.tsv`, `concat_tree.nwk`, `species_metadata.json` and the four
figures in PNG and PDF.

To reproduce everything from scratch, in order:

```bash
PY=/software/shared/apps/x86_64/python/3.11.4/bin/python3
cd /group/esb/cesen/1kp/code/1kp_wgms/others/calibration
$PY paralog_ks_bin_heights.py          # ~10 min, scans 1178 paranome runs
$PY indicator_auc.py
$PY build_calibration_pairs.py         # needs a current ksrates inspect-ks-db export
$PY plot_calibration.py
$PY alternative_predictors.py
$PY plot_lycophyte_sister.py
$PY plot_sister_clade_distributions.py
$PY recommend_ortholog_pairs.py
```

Two notes on reproducibility.

`build_calibration_pairs.py` always takes the most recent
`ks_database_server_address_ortholog_peaks_*.tsv` in each of the two database folders, so the
export has to be refreshed first:

```bash
cd /group/esb/cesen/1kp/ks_analysis/1kp/ks_database
singularity exec -B /group/esb/cesen/1kp \
    /group/esb/cesen/1kp/software/ksrates_20261022_3add8d4.sif \
    ksrates inspect-ks-db ks_database_server_address.txt --type ortholog
```

`recommend_ortholog_pairs.py` excludes pairs that are already in the calibration, so it
proposes the *next* round each time it runs; it does not reconstruct a past one. The set that
was actually run for this document is recorded in
`ks_analysis/1kp/calibration_run/ortholog_pairs.tsv`.

The root-to-tip script, and `calib_common.extract_concat_newick`, both refuse to proceed when
all terminal branch lengths are identical, which is the signature of the ASTRAL tree and the
mistake section 1 guards against.
