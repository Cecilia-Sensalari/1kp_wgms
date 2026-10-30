#!/usr/bin/env python3
"""The expected-Ks distribution of one sister clade, worked through in detail.

Test branch = the lycophyte stem, so clade B is every euphyllophyte. Panel A is the
histogram by lineage, panel B what each candidate budget rule selects. Also simulates how
reliable a within-clade ranking by expected Ks actually is, given the calibration's own error.

Writes lycophyte_sister_candidates.tsv and expected_ks_distribution_lycophyte_sister.png/pdf.

    python3 plot_lycophyte_sister.py
"""
import os, sys, argparse, statistics
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

TARGET = "Lycopodiopsida"
SHOW = {"fern": "#2a78d6", "gymnosperm": "#eb6834", "angiosperm": "#1baf7a"}
NICE = {"fern": "ferns", "gymnosperm": "gymnosperms", "angiosperm": "angiosperms"}
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb"
DELTA = 0.5
# The calibration's error is MULTIPLICATIVE, not additive: in Ks units its SD runs from 0.08
# (Ks 0-1) to 0.93 (Ks > 3), a twelvefold range, while the SD of log(observed/predicted) stays
# near 0.16 across every band. Of that log-variance about 40% is a property of one species and
# so is shared by all of its pairs; within a branch the A side is fixed and that share cancels,
# leaving sqrt(total - species) = 0.125 log units - roughly 13% of the predicted value - as the
# error that actually scrambles a ranking of candidates. See alternative_predictors.py and
# section 5.2 of the methods document.
LOG_NOISE_SD = C.LOG_NOISE_SD   # refitted on the trusted pairs; see calib_common


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default=C.DATA_DIR)
    args = ap.parse_args()

    _, dist = C.load_concat_tree()
    meta = C.species_metadata()
    median_bin = C.median_bin_heights()
    group = C.groups(meta)
    lin = C.lineages()
    code_to_tip = {info["code"]: tip for tip, info in meta.items()}
    nodes = C.branch_nodes(C.load_analysis_tree())

    members = {code_to_tip.get(c, c) for c, lineage in lin.items() if TARGET in lineage}
    node = C.clade_containing(nodes, members)
    a_species = C.pick_A(node, median_bin)
    clade = max(C.sister_clades(node, 1), key=lambda s: len(s.leaves()))
    print(f"A representatives: {a_species}")

    rows = []
    for candidate in clade.leaves():
        ks = C.expected_ks(candidate, a_species, dist)
        if ks is not None:
            rows.append((candidate, group.get(candidate), ks, median_bin.get(candidate, 0)))
    rows.sort(key=lambda r: r[2])
    out = os.path.join(args.out_dir, "lycophyte_sister_candidates.tsv")
    with open(out, "w") as fh:
        fh.write("code\tgroup\texpected_ks_model\texpected_ks_reported\t"
                 "median_bin_height\n")
        for code, grp, ks, bar in rows:
            fh.write(f"{code}\t{grp}\t{ks:.3f}\t{C.report_ks(ks)}\t{bar:.1f}\n")
    print(f"{len(rows)} candidates -> {out}")
    for name in SHOW:
        v = sorted(r[2] for r in rows if r[1] == name)
        if v:
            print(f"  {name:<12} n={len(v):>3}  min {v[0]:.2f}  median {statistics.median(v):.2f}"
                  f"  max {v[-1]:.2f}")

    values = np.array([r[2] for r in rows])
    bars_h = np.array([r[3] for r in rows])
    groups_of = np.array([r[1] for r in rows])
    quartile = np.percentile(values, 25)
    budget = max(2.0, values.min() + DELTA)
    pick_q = int(np.argmax(np.where(values <= quartile, bars_h, -1)))
    pick_b = int(np.argmax(np.where(values <= budget, bars_h, -1)))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.9), facecolor=SURFACE,
                             gridspec_kw=dict(width_ratios=[1, 1.08]))
    for ax in axes:
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.grid(True, color=GRID, lw=.6, alpha=.85)
        ax.set_axisbelow(True)
        ax.tick_params(colors=INK_2, labelsize=9)

    ax = axes[0]
    edges = np.arange(values.min() - 0.06, values.max() + 0.08, 0.06)
    bottom = np.zeros(len(edges) - 1)
    for name in ("fern", "gymnosperm", "angiosperm"):
        counts = np.histogram(values[groups_of == name], bins=edges)[0]
        ax.bar(edges[:-1], counts, width=np.diff(edges) * 0.86, bottom=bottom,
               color=SHOW[name], linewidth=0, align="edge", zorder=3)
        bottom = bottom + counts
    ax.set_ylim(0, (bottom.max() or 1) * 1.22)
    for i, name in enumerate(("fern", "gymnosperm", "angiosperm")):
        v = values[groups_of == name]
        y = ax.get_ylim()[1] * (0.96 - i * 0.068)
        ax.scatter([edges[0] + 0.01], [y], s=58, c=SHOW[name], marker="s",
                   edgecolors=SURFACE, linewidths=.8, zorder=6)
        ax.text(edges[0] + 0.07, y, f"{NICE[name]}   n={len(v)},  median {np.median(v):.2f}",
                color=INK_2, fontsize=9.3, va="center", zorder=6,
                bbox=dict(fc=SURFACE, ec="none", pad=1.8, alpha=.95))
    ax.axvline(quartile, color=INK, lw=1.6, ls=":", zorder=5)
    ax.axvline(budget, color=INK, lw=1.6, ls="--", zorder=5)
    ax.annotate(f"lowest 25% of the clade\nends at {quartile:.2f}", xy=(quartile, bottom.max() * .72),
                xytext=(quartile + 0.06, bottom.max() * .82), color=INK, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=INK, lw=.9))
    ax.annotate(f"min+{DELTA}\n= {budget:.2f}", xy=(budget, bottom.max() * .30),
                xytext=(budget - 0.52, bottom.max() * .46), color=INK, fontsize=9,
                ha="right", arrowprops=dict(arrowstyle="-", color=INK, lw=.9),
                bbox=dict(fc=SURFACE, ec="none", pad=1.8, alpha=.95))
    ax.set_xlabel("expected ortholog Ks to the lycophyte A-species",
                  color=INK, fontsize=10.5)
    ax.set_ylabel("number of candidate species", color=INK, fontsize=10.5)
    ax.set_title("A   The B clade is bimodal, and the modes are rates, not depths",
                 color=INK, fontsize=11.5, loc="left", pad=9)

    ax = axes[1]
    ax.axvspan(values.min() - 0.02, quartile, color=INK_2, alpha=.055, lw=0)
    ax.axvspan(values.min() - 0.02, budget, color=INK_2, alpha=.085, lw=0)
    for name in ("fern", "gymnosperm", "angiosperm"):
        m = groups_of == name
        ax.scatter(values[m], bars_h[m], s=26, c=SHOW[name], edgecolors=SURFACE,
                   linewidths=.7, alpha=.9, zorder=3, label=f"{NICE[name]} (n={m.sum()})")
    for index, text, dy in ((pick_q, "25%-quantile rule picks", 95), (pick_b, "min+0.5 rule picks", 140)):
        ax.scatter([values[index]], [bars_h[index]], s=190, facecolors="none",
                   edgecolors=INK, linewidths=2, zorder=6)
        ax.annotate(f"{text}\n{rows[index][0]}\nexpKs {values[index]:.2f}, "
                    f"bin height {bars_h[index]:.0f}",
                    xy=(values[index], bars_h[index]),
                    xytext=(values[index] + (0.13 if dy == 95 else -0.36), bars_h[index] + dy),
                    color=INK, fontsize=9, arrowprops=dict(arrowstyle="-", color=INK, lw=.9))
    ax.set_xlabel("expected ortholog Ks", color=INK, fontsize=10.5)
    ax.set_ylabel("median bin height  (data available)", color=INK, fontsize=10.5)
    ax.set_title("B   What each rule selects", color=INK, fontsize=11.5, loc="left", pad=9)
    ax.legend(frameon=False, fontsize=9, loc="upper right", labelcolor=INK_2)

    fig.suptitle(f"Expected ortholog Ks of the {len(rows)} candidates in the sister clade "
                 f"of the lycophytes", color=INK, fontsize=13, x=0.011, ha="left", y=.985)
    fig.tight_layout(rect=[0, 0, 1, .945])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.out_dir, "expected_ks_distribution_lycophyte_sister." + ext),
                    dpi=190, facecolor=SURFACE)
    print("wrote", os.path.join(args.out_dir, "expected_ks_distribution_lycophyte_sister.png"))

    rng = np.random.default_rng(0)
    draws = 20000
    print("\n=== how reliable is a within-clade ranking by expected Ks? ===")
    print(f"(multiplicative noise, SD {LOG_NOISE_SD} in log units = "
          f"{100 * (np.exp(LOG_NOISE_SD) - 1):.0f}% of the predicted value, "
          f"{draws} replicates)")
    print(f"{'clade':<34}{'n':>5}{'level':>8}{'spread':>9}{'P(pick is lowest)':>20}"
          f"{'true pct. of pick':>20}")

    def simulate(label, truth):
        truth = np.asarray(truth, dtype=float)
        n = len(truth)
        observed = truth * np.exp(rng.normal(0, LOG_NOISE_SD, size=(draws, n)))
        chosen = np.take_along_axis(np.tile(truth, (draws, 1)),
                                    observed.argmin(1)[:, None], 1)[:, 0]
        hit = (truth[observed.argmin(1)] == truth.min()).mean()
        percentile = 100 * ((truth < chosen[:, None]).sum(1) / (n - 1)).mean()
        print(f"{label:<34}{n:>5}{np.median(truth):>8.2f}"
              f"{truth.max() - truth.min():>9.2f}{hit:>20.2f}{percentile:>19.0f}%")

    # a typical B/C/D slot: the median of all 2539 slots with at least four candidates is
    # 10 candidates, a level of 0.91 Ks and a spread of 33% of that level
    simulate("typical B/C/D slot (median of all)",
             np.linspace(0.91 * (1 - 0.33 / 2), 0.91 * (1 + 0.33 / 2), 10))
    simulate("this clade, after min+%.1f" % DELTA, values[values <= budget])
    simulate("this clade, all candidates", values)
    print("\nRandom choice would put the pick at the 50th percentile.")


if __name__ == "__main__":
    main()
