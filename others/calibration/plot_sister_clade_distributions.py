#!/usr/bin/env python3
"""Expected ortholog Ks of the candidates in the sister clade, for twelve named branches.

One panel per test branch, from Cupressaceae to the chlorophytes, with the histogram
coloured by lineage and the two candidate budget rules drawn on it. Also reports how many
B/C/D slots across all 2354 branches actually offer a choice.

    python3 plot_sister_clade_distributions.py
"""
import os, sys, argparse, statistics, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

# validated categorical slots 1-6, assigned in this fixed order
COLOUR = dict(zip(C.GROUP_ORDER,
                  ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]))
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb"

PANELS = [("Cupressaceae", "Cupressaceae"), ("Solanaceae", "Solanaceae"),
          ("Pinaceae", "Pinaceae"), ("Lamiales", "Lamiales"),
          ("Caryophyllales", "Caryophyllales"), ("Ranunculales", "Ranunculales"),
          ("Monocots", "Liliopsida"), ("Liverworts", "Marchantiophyta"),
          ("Hornworts", "Anthocerotophyta"), ("Gymnosperms", "Acrogymnospermae"),
          ("Lycophytes", "Lycopodiopsida"), ("Chlorophytes", "Chlorophyta")]
DELTA = 0.5


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
    tree = C.load_analysis_tree()
    nodes = C.branch_nodes(tree)

    def members(key):
        return {code_to_tip.get(c, c) for c, lineage in lin.items() if key in lineage}

    data = []
    for label, key in PANELS:
        node = C.clade_containing(nodes, members(key))
        if node is None:
            print(f"{label}: no clean clade, skipped", file=sys.stderr)
            continue
        sisters = C.sister_clades(node, 1)
        if not sisters:
            continue
        a_species = C.pick_A(node, median_bin)
        rows = []
        for candidate in max(sisters, key=lambda s: len(s.leaves())).leaves():
            p = C.median_path(candidate, a_species, dist)
            if p is not None:
                rows.append((candidate, C.path_to_ks(p), median_bin.get(candidate, 0),
                             group.get(candidate), p))
        if len(rows) >= 6:
            data.append((label, rows))

    fig, axes = plt.subplots(4, 3, figsize=(14.5, 13.2), facecolor=SURFACE)
    for ax, (label, rows) in zip(axes.ravel(), data):
        ax.set_facecolor(SURFACE)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
        ax.grid(True, color=GRID, lw=.6, alpha=.85)
        ax.set_axisbelow(True)
        ax.tick_params(colors=INK_2, labelsize=8.5)

        values = np.array([r[1] for r in rows])
        groups_of = np.array([r[3] for r in rows])
        paths = np.array([r[4] for r in rows])
        low, high = values.min(), values.max()
        pad = max(0.04, (high - low) * 0.10)
        edges = np.linspace(low - pad, high + pad, 26)
        bottom = np.zeros(len(edges) - 1)
        for name in C.GROUP_ORDER:
            counts = np.histogram(values[groups_of == name], bins=edges)[0]
            if counts.sum() == 0:
                continue
            ax.bar(edges[:-1], counts, width=np.diff(edges) * 0.86, bottom=bottom,
                   color=COLOUR[name], linewidth=0, align="edge", zorder=3)
            bottom = bottom + counts

        # the two reliability thresholds, drawn in the plot's Ks units: past the first the
        # ortholog Ks peak is a broad mode, past the second it is a flat smear and the number
        # on this axis is an extrapolation rather than anything that was ever measured
        grey_from, unreliable_from = C.path_to_ks(C.RELIABLE_PATH), C.path_to_ks(C.GREY_PATH)
        for start, stop, alpha, text in (
                (grey_from, min(unreliable_from, high + pad), .07, "broad peak"),
                (unreliable_from, high + pad, .16, "no measurable\npeak")):
            if high + pad > start:
                ax.axvspan(start, stop, color=INK_2, alpha=alpha, lw=0, zorder=1)
                # only label a band wide enough to hold the words, and sit it low in the
                # panel where the bars have usually run out
                if (stop - start) > (high + pad - (low - pad)) * .17:
                    ax.text(start + (stop - start) * .05, .30, text,
                            transform=ax.get_xaxis_transform(), color=INK_2, fontsize=8,
                            va="center", zorder=7, clip_on=True,
                            bbox=dict(fc=SURFACE, ec="none", pad=1.4, alpha=.85))

        quartile = np.percentile(values, 25)
        budget = max(2.0, low + DELTA)
        ax.axvline(quartile, color=INK, lw=1.5, ls=":", zorder=5)
        if budget <= high + pad:
            ax.axvline(budget, color=INK, lw=1.5, ls="--", zorder=5)
        ax.set_ylim(0, (bottom.max() or 1) * 1.42)
        ax.set_title(f"{label}   n={len(values)} candidates", color=INK, fontsize=10.5,
                     loc="left", pad=6)
        ax.set_xlim(low - pad, high + pad)
        lo_path, hi_path = paths[values.argmin()], paths[values.argmax()]
        reliable = C.path_reliability(hi_path) == "reliable"
        shown_spread = f"{high - low:.2f}" if reliable else "-"
        for y, text in ((.955, f"expKs {C.report_ks(low, lo_path)} - "
                               f"{C.report_ks(high, hi_path)}   spread {shown_spread}"),
                        (.875, f"lowest 25% keeps {int((values <= quartile).sum())}   ·   "
                               f"min+{DELTA} keeps {int((values <= budget).sum())}")):
            # clip_on keeps a long line from enlarging the axes' tight bounding box, which
            # tight_layout would otherwise honour by shrinking every panel
            ax.text(.02, y, text, transform=ax.transAxes, color=INK_2, fontsize=8.6,
                    va="top", zorder=7, clip_on=True,
                    bbox=dict(fc=SURFACE, ec="none", pad=1.6, alpha=.92))
        ax.set_xlabel("expected ortholog Ks  (power law)", color=INK, fontsize=9)
    for ax in axes[:, 0]:
        ax.set_ylabel("candidate species", color=INK, fontsize=9)

    handles = [plt.Line2D([], [], marker="s", ls="", ms=9, color=COLOUR[g]) for g in C.GROUP_ORDER]
    fig.legend(handles, C.GROUP_ORDER, frameon=False, ncol=6, loc="upper right",
               bbox_to_anchor=(.995, .978), fontsize=10, labelcolor=INK_2,
               handletextpad=.4, columnspacing=1.6)
    fig.suptitle("Expected ortholog Ks of the candidates in the sister clade, for twelve "
                 "branches of the 1KP tree", color=INK, fontsize=14, x=.008, ha="left", y=.985)
    fig.text(.008, .962, f"dotted line = lowest 25% of the clade   ·   dashed line = min + {DELTA}",
             color=INK_2, fontsize=9.5, ha="left")
    fig.tight_layout(rect=[0, 0, 1, .952])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.out_dir, "expected_ks_distributions_twelve_branches." + ext),
                    dpi=165, facecolor=SURFACE)
    print("wrote", os.path.join(args.out_dir, "expected_ks_distributions_twelve_branches.png"))

    sizes = collections.Counter()
    total = 0
    for node in nodes:
        for clade in C.sister_clades(node):
            total += 1
            k = len(clade.leaves())
            sizes["1" if k == 1 else "2-3" if k <= 3 else "4-10" if k <= 10
                  else "11-50" if k <= 50 else ">50"] += 1
    print(f"\nB/C/D slots over all {len(nodes)} branches: {total}")
    for key in ("1", "2-3", "4-10", "11-50", ">50"):
        print(f"  clade size {key:<6} {sizes[key]:>5}  ({100 * sizes[key] / total:4.1f}%)")


if __name__ == "__main__":
    main()
