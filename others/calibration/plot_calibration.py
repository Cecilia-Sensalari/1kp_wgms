#!/usr/bin/env python3
"""Fit and plot ortholog Ks against concatenation-tree path distance.

Fits four models to ks_analysis/1kp/calibration/calibration_pairs.tsv, all with
path distance as the response and the ortholog Ks peak as the predictor:

    1. linear through the origin, fitted on pairs with Ks <= 1
    2. linear through the origin, fitted on pairs with Ks <= 2
    3. power law          path = A * Ks^c          (c < 1 == concave)
    4. saturating         path = Pmax * (1 - e^(-k Ks))

All four are fitted ONLY on the pairs with path distance <= C.RELIABLE_PATH. Beyond that the
ortholog Ks distribution flattens into a smear against the max_ks_orthologs cutoff and the
peak the KDE returns is an artefact of the cutoff rather than a measurement, so fitting to
those points would be fitting to noise (see the methods document, 3.5). The excluded pairs are
still drawn, in grey, because showing where the measurement fails is the point.

and reports, for each, the error of the Ks it predicts back from a path distance.
Panel A is the point cloud with the fits; panel B the prediction error.

    python3 plot_calibration.py
"""
import os, sys, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import calib_common as C

SERIES_1, SERIES_2 = "#2a78d6", "#eb6834"        # validated categorical slots 1 and 2
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb"
BANDS = [(0, 1), (1, 2), (2, 3), (3, 99)]


def band_label(lo, hi):
    return f"{lo}-{hi}" if hi < 99 else f"{lo}+"


def load(path):
    ks, sd, dist, source = [], [], [], []
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            f = line.rstrip("\n").split("\t")
            ks.append(float(f[2])); sd.append(float(f[3]))
            dist.append(float(f[4])); source.append(f[5])
    return np.array(ks), np.array(sd), np.array(dist), np.array(source)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default=C.CALIBRATION_PAIRS)
    ap.add_argument("--out-dir", default=C.DATA_DIR)
    args = ap.parse_args()

    ks_all, sd_all, dist_all, source_all = load(args.pairs)
    trusted = dist_all <= C.RELIABLE_PATH
    ks, sd, dist, source = ks_all[trusted], sd_all[trusted], dist_all[trusted], source_all[trusted]
    n, n_excluded = len(ks), int((~trusted).sum())

    # ---- the four models, all fitted on the trusted pairs only ------------------------
    def slope_through_origin(mask):
        return (ks[mask] * dist[mask]).sum() / (ks[mask] ** 2).sum()

    under_1, under_2 = ks <= 1.0, ks <= 2.0
    slope_1, slope_2 = slope_through_origin(under_1), slope_through_origin(under_2)
    log_slope, log_intercept = np.polyfit(np.log(ks), np.log(dist), 1)
    power_a, power_c = np.exp(log_intercept), log_slope
    (sat_max, sat_k), _ = curve_fit(lambda x, p, k: p * (1 - np.exp(-k * x)),
                                    ks, dist, p0=[2.5, 0.6], maxfev=200000)

    def from_linear(d, slope):
        return d / slope

    def from_power(d):
        return (d / power_a) ** (1.0 / power_c)

    def from_saturating(d):
        return -np.log(np.clip(1 - d / sat_max, 1e-9, None)) / sat_k

    print(f"n pairs fitted = {n}   (1kp {int((source == '1kp').sum())}, "
          f"master1 {int((source == 'master1').sum())});  "
          f"{n_excluded} excluded at path > {C.RELIABLE_PATH}")
    print(f"adopted constants in calib_common: POWER_A = {power_a:.4f}, "
          f"POWER_EXPONENT = {1 / power_c:.3f}")
    print("\n=== fits:  y = concat path distance,  x = ortholog Ks peak ===")
    print(f"linear thru origin, Ks<=1 (n={under_1.sum():3d}) : path = Ks / {1 / slope_1:.2f}")
    print(f"linear thru origin, Ks<=2 (n={under_2.sum():3d}) : path = Ks / {1 / slope_2:.2f}")
    print(f"power law        (all n={n}) : path = {power_a:.4f} * Ks^{power_c:.3f}")
    print(f"saturating       (all n={n}) : path = {sat_max:.3f}*(1-exp(-{sat_k:.3f}*Ks))")

    residuals = {f"linear 1/{1 / slope_1:.2f}": from_linear(dist, slope_1) - ks,
                 f"linear 1/{1 / slope_2:.2f}": from_linear(dist, slope_2) - ks,
                 "power law": from_power(dist) - ks,
                 "saturating": from_saturating(dist) - ks}

    print("\n=== error of PREDICTED Ks (median |err| / SD), by observed-Ks band ===")
    header = "".join(f"{'Ks ' + band_label(*b):>17}" for b in BANDS)
    print(f"{'model':<16}{header}{'all':>17}")
    for name, res in residuals.items():
        row = f"{name:<16}"
        for lo, hi in BANDS:
            m = (ks >= lo) & (ks < hi)
            row += f"{np.median(np.abs(res[m])):>8.2f}/{np.std(res[m]):>8.2f}" if m.sum() \
                else f"{'-':>17}"
        print(row + f"{np.median(np.abs(res)):>8.2f}/{np.std(res):>8.2f}")
    print("n per band:", {band_label(lo, hi): int(((ks >= lo) & (ks < hi)).sum())
                          for lo, hi in BANDS})

    print("\n=== local conversion factor  Ks/path  by band ===")
    print(f"{'band':<8}{'n':>4}{'median':>9}{'IQR':>20}{'CV%':>8}{'min-max':>18}")
    for lo, hi in BANDS:
        m = (ks >= lo) & (ks < hi)
        if not m.sum():
            continue
        factor = ks[m] / dist[m]
        q1, q3 = np.percentile(factor, [25, 75])
        print(f"{band_label(lo, hi):<8}{m.sum():>4}{np.median(factor):>9.2f}"
              f"{f'{q1:.2f}-{q3:.2f}':>20}{100 * np.std(factor) / np.mean(factor):>8.1f}"
              f"{f'{factor.min():.2f}-{factor.max():.2f}':>18}")

    # ---- figure ----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.8), facecolor=SURFACE)
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
    # Drawn in the direction the calibration is USED: path in, Ks out. The models are fitted the
    # other way round (path as the response, which is how the pairs were measured and which puts
    # the bootstrap error on the response), so each curve below is the fitted relation inverted.
    # Showing the fitted orientation instead needs both formulas in the legend to stay honest,
    # since path(Ks) flattens while its inverse Ks(path) accelerates - the same curve, read twice.
    grid = np.linspace(1e-3, dist_all.max() * 1.04, 600)
    ax.plot(grid, grid / slope_1, color=INK_2, lw=1.5, ls="--", alpha=.75,
            label=f"linear, fitted on Ks$\\leq$1   (Ks = {1 / slope_1:.2f}$\\times$path)")
    ax.plot(grid, grid / slope_2, color=INK_2, lw=1.5, ls=":", alpha=.75,
            label=f"linear, fitted on Ks$\\leq$2   (Ks = {1 / slope_2:.2f}$\\times$path)")
    # the saturating fit has an asymptote in path at sat_max, so its inverse only exists below it
    sat_grid = grid[grid < sat_max * 0.999]
    ax.plot(sat_grid, -np.log(1 - sat_grid / sat_max) / sat_k, color=INK_2, lw=1.5, ls="-.",
            alpha=.75, label="saturating (unidentifiable, not used)")
    ax.plot(grid, (grid / power_a) ** (1 / power_c), color=INK, lw=2.6,
            label=f"ADOPTED  ·  power law\n"
                  f"   Ks = (path/{power_a:.4f})$^{{{1 / power_c:.3f}}}$")

    right = dist_all.max() * 1.05
    ax.axvspan(C.RELIABLE_PATH, C.BORDERLINE_PATH, color="#c8a23a", alpha=.13, lw=0, zorder=0)
    ax.axvspan(C.BORDERLINE_PATH, right, color=INK_2, alpha=.17, lw=0, zorder=0)
    top = ks_all.max() * 1.06
    for x, text in ((C.RELIABLE_PATH, f"path {C.RELIABLE_PATH}  ·  peak broadens, still locatable"),
                    (C.BORDERLINE_PATH, f"path {C.BORDERLINE_PATH}  ·  flat smear, these Ks are not measurements")):
        ax.axvline(x, color=INK_2, lw=1.1, ls=(0, (1, 2)), zorder=4)
        ax.text(x - right * .008, 0.02, text, color=INK, fontsize=8.3, rotation=90,
                ha="right", va="bottom", transform=ax.get_xaxis_transform(), zorder=8,
                bbox=dict(fc=SURFACE, ec="none", pad=1.6, alpha=.88))

    ax.scatter(dist_all[~trusted], ks_all[~trusted], s=34, c="none", marker="o",
               edgecolors=INK_2, linewidths=1.0, alpha=.55, zorder=2,
               label=f"excluded from the fit  (n={n_excluded})")
    for tag, colour, marker, label in (("1kp", SERIES_1, "o", "1kp database"),
                                       ("master1", SERIES_2, "^", "master_1_2026 database")):
        m = source == tag
        ax.errorbar(dist[m], ks[m], yerr=sd[m], fmt="none", ecolor=colour, elinewidth=.9, alpha=.45)
        ax.scatter(dist[m], ks[m], s=40, c=colour, marker=marker, edgecolors=SURFACE,
                   linewidths=1.1, zorder=3, label=f"{label}  (n={m.sum()})")
    ax.set_xlim(-0.03, right)
    ax.set_ylim(-0.15, top)
    ax.set_xlabel("concat-tree path distance (subst./site)   —   the input", color=INK, fontsize=10.5)
    ax.set_ylabel("ortholog Ks peak (mode; bars = bootstrap SD)", color=INK, fontsize=10.5)
    ax.set_title("A   The point cloud, in the direction the calibration is used",
                 color=INK, fontsize=12, loc="left", pad=9)
    ax.legend(frameon=True, framealpha=.88, edgecolor="none", facecolor=SURFACE,
              fontsize=8.5, loc="upper left", labelcolor=INK_2)

    ax = axes[1]
    ax.axhspan(-0.25, 0.25, color=INK_2, alpha=.07, lw=0)
    for name, colour, marker in ((f"linear 1/{1 / slope_2:.2f}", SERIES_2, "x"),
                                 ("power law", SERIES_1, "o")):
        ax.scatter(ks, residuals[name], s=36, c=colour, marker=marker,
                   edgecolors=SURFACE if marker == "o" else None, linewidths=1.0, alpha=.9,
                   zorder=3,
                   label=f"{name}   (SD {np.std(residuals[name]):.2f} Ks)"
                         + ("   <- adopted" if name == "power law" else ""))
    ax.axhline(0, color=INK, lw=1.2)
    ax.set_xlabel("observed ortholog Ks peak", color=INK, fontsize=10.5)
    ax.set_ylabel("predicted $-$ observed Ks", color=INK, fontsize=10.5)
    ax.set_title("B   Error when predicting Ks from path distance", color=INK, fontsize=12,
                 loc="left", pad=9)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left", labelcolor=INK_2)

    fig.suptitle(f"Ortholog Ks peak vs concat-tree path distance  ·  {n} pairs fitted, "
                 f"{n_excluded} too deep to measure", color=INK, fontsize=13, x=0.011,
                 ha="left", y=.985)
    fig.tight_layout(rect=[0, 0, 1, .945])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(args.out_dir, "ortholog_ks_vs_path_distance." + ext),
                    dpi=190, facecolor=SURFACE)
    print("\nwrote", os.path.join(args.out_dir, "ortholog_ks_vs_path_distance.png"))


if __name__ == "__main__":
    main()
