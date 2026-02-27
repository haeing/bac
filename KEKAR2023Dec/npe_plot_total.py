#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Make heatmaps for:
  - Experiment Np.e. map (from kek_result.root)
  - Simulation Np.e. map (Gaussian mean of nhMppc from Geant4 ROOT files)
And make:
  - Relative-difference map: (Data - Sim) / Data
  - (Optional) Absolute-difference map: Data - Sim

Notes
-----
- Simulation files are assumed like: .../kek_aerogel{2|3}_x{m?dd}_y{m?dd}.root
- Experiment selection uses HV=58 and threshold rule: {2:60, 3:75}
- Comparison excludes y = ±54 (outside aerogel), by default.
"""

import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

# ---------- EXP (PyROOT) ----------
import ROOT
ROOT.gROOT.SetBatch(True)

# ---------- SIM (uproot + scipy) ----------
import uproot
from scipy.optimize import curve_fit

# =========================
# CONFIG
# =========================
# Experiment ROOT
EXP_FILE = "kek_result.root"
EXP_TREE = "tree"
BR_X     = "x_pos"
BR_Y     = "y_pos"
BR_THICK = "bac_thick"
BR_HV    = "bac_HV"
BR_THRE  = "bac_thre"
BR_NPE   = "bac_npe"

HV_TARGET = 58
THRE_BY_THICK = {2: 60, 3: 75}

# Simulation ROOT file pattern (capture aerogel2/3 automatically)
SIM_PATTERN = "../../../data/KEKAR2023Dec/g4_root/kek_aerogel*_x*_y*.root"
SIM_TREE = "tree"
SIM_BR   = "nhMppc"

# Fixed grid (mm)
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=int)
ys = np.array([-54, -36, -18, 0, 18, 36, 54], dtype=int)

# Comparison mask (exclude aerogel-outside rows)
EXCLUDE_Y = set([54, -54])   # set() to compare all

# Gaussian-fit settings for simulation
NBINS = 120
FIT_WINDOW_SIGMA = 2.0

# Output directory
OUTDIR = "plots_compare"
os.makedirs(OUTDIR, exist_ok=True)


# =========================
# Utilities
# =========================
def parse_xy_aerogel(fname: str):
    """
    Example:
      kek_aerogel3_xm16_ym54.root -> aerogel=3, x=-16, y=-54
      kek_aerogel2_x16_y54.root   -> aerogel=2, x=16,  y=54
    """
    m = re.search(r"aerogel(?P<a>\d+)_x(?P<x>m?\d+)_y(?P<y>m?\d+)\.root$", os.path.basename(fname))
    if not m:
        raise ValueError(f"Cannot parse aerogel/x/y from: {fname}")

    def tok2int(tok):
        return -int(tok[1:]) if tok.startswith("m") else int(tok)

    a = int(m.group("a"))
    x = tok2int(m.group("x"))
    y = tok2int(m.group("y"))
    return a, x, y


def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fit_gaussian_mean(data):
    """
    data: 1D numpy array (nhMppc)
    return: (mu, mu_err)
    """
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if len(data) < 50:
        return np.nan, np.nan

    xmin, xmax = np.percentile(data, [0.5, 99.5])
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmin >= xmax:
        return np.nan, np.nan

    counts, edges = np.histogram(data, bins=NBINS, range=(xmin, xmax))
    centers = 0.5 * (edges[:-1] + edges[1:])

    ipeak = int(np.argmax(counts))
    mu0 = float(centers[ipeak])

    sigma0 = float(np.std(data))
    if not np.isfinite(sigma0) or sigma0 <= 0:
        sigma0 = max((xmax - xmin) / 20.0, 1.0)

    lo = mu0 - FIT_WINDOW_SIGMA * sigma0
    hi = mu0 + FIT_WINDOW_SIGMA * sigma0
    mask = (centers >= lo) & (centers <= hi) & (counts > 0)

    xfit = centers[mask]
    yfit = counts[mask]
    if len(xfit) < 8:
        return np.nan, np.nan

    yerr = np.sqrt(yfit)
    yerr[yerr == 0] = 1.0
    A0 = float(np.max(yfit))

    try:
        popt, pcov = curve_fit(
            gauss, xfit, yfit,
            p0=[A0, mu0, sigma0],
            sigma=yerr,
            absolute_sigma=True,
            maxfev=20000
        )
        mu = float(popt[1])
        mu_err = float(np.sqrt(pcov[1, 1])) if (pcov is not None and pcov.shape == (3, 3)) else np.nan
        return mu, mu_err
    except Exception:
        return np.nan, np.nan


def plot_heatmap(m, title, outname, cbar_label="Np.e.", fmt="{:.2f}", vmin=None, vmax=None):
    fig, ax = plt.subplots(figsize=(7.2, 6.2))

    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]

    im = ax.imshow(
        m, origin="lower", aspect="auto",
        extent=[xs[0]-dx/2, xs[-1]+dx/2, ys[0]-dy/2, ys[-1]+dy/2],
        vmin=vmin, vmax=vmax
    )

    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title(title)
    ax.set_xticks(xs)
    ax.set_yticks(ys)

    # annotate cell values
    for iy, yy in enumerate(ys):
        for ix, xx in enumerate(xs):
            if np.isfinite(m[iy, ix]):
                ax.text(xx, yy, fmt.format(m[iy, ix]),
                        ha="center", va="center", fontsize=9, color="black")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)

    fig.tight_layout()
    fig.savefig(outname, dpi=200)
    plt.close(fig)
    print(f"Saved: {outname}")


def apply_compare_mask(m):
    """Return a masked copy excluding y in EXCLUDE_Y."""
    m2 = np.array(m, copy=True)
    for iy, yy in enumerate(ys):
        if yy in EXCLUDE_Y:
            m2[iy, :] = np.nan
    return m2


# =========================
# Experiment map
# =========================
def load_experiment_arrays():
    f = ROOT.TFile.Open(EXP_FILE)
    t = f.Get(EXP_TREE)
    if not t:
        raise RuntimeError(f"Cannot find TTree '{EXP_TREE}' in {EXP_FILE}")

    x, y, thick, hv, thre, npe = [], [], [], [], [], []
    for ev in t:
        x.append(int(getattr(ev, BR_X)))
        y.append(int(getattr(ev, BR_Y)))
        thick.append(int(getattr(ev, BR_THICK)))
        hv.append(int(getattr(ev, BR_HV)))
        thre.append(int(getattr(ev, BR_THRE)))
        npe.append(float(getattr(ev, BR_NPE)))

    return (np.array(x, dtype=int),
            np.array(y, dtype=int),
            np.array(thick, dtype=int),
            np.array(hv, dtype=int),
            np.array(thre, dtype=int),
            np.array(npe, dtype=float))


def make_exp_map(th, x, y, thick, hv, thre, npe):
    if th not in THRE_BY_THICK:
        raise ValueError(f"No threshold rule for thickness={th}")
    thre_target = THRE_BY_THICK[th]

    sel = (hv == HV_TARGET) & (thick == th) & (thre == thre_target)
    print(f"[EXP thick={th}] select HV=={HV_TARGET} & thre=={thre_target} -> {int(np.sum(sel))} entries")

    m = np.full((len(ys), len(xs)), np.nan, dtype=float)
    for iy, yy in enumerate(ys):
        for ix, xx in enumerate(xs):
            mask = sel & (x == xx) & (y == yy)
            if np.any(mask):
                m[iy, ix] = float(np.mean(npe[mask]))
    return m


# =========================
# Simulation map (aerogel 2/3)
# =========================
def make_sim_map(aerogel_id):
    files = sorted(glob.glob(SIM_PATTERN))
    if not files:
        raise RuntimeError(f"No files matched: {SIM_PATTERN}")

    mean_map = np.full((len(ys), len(xs)), np.nan, dtype=float)
    err_map  = np.full((len(ys), len(xs)), np.nan, dtype=float)

    nused = 0
    for fn in files:
        try:
            a, x, y = parse_xy_aerogel(fn)
        except ValueError:
            continue
        if a != aerogel_id:
            continue
        if x not in xs or y not in ys:
            continue

        with uproot.open(fn) as f:
            t = f[SIM_TREE]
            data = t[SIM_BR].array(library="np")

        mu, mu_err = fit_gaussian_mean(data)

        iy = int(np.where(ys == y)[0][0])
        ix = int(np.where(xs == x)[0][0])
        mean_map[iy, ix] = mu
        err_map[iy, ix]  = mu_err

        nused += 1

    print(f"[SIM aerogel={aerogel_id}] used files: {nused}")
    return mean_map, err_map


# =========================
# Compare + Plot
# =========================
def make_relative_diff(data_map, sim_map):
    """
    rel = (data - sim) / data
    """
    rel = np.full_like(data_map, np.nan, dtype=float)
    good = np.isfinite(data_map) & np.isfinite(sim_map) & (data_map != 0.0)
    rel[good] = (data_map[good] - sim_map[good]) / data_map[good]
    return rel


def make_absolute_diff(data_map, sim_map):
    diff = np.full_like(data_map, np.nan, dtype=float)
    good = np.isfinite(data_map) & np.isfinite(sim_map)
    diff[good] = data_map[good] - sim_map[good]
    return diff


def main():
    # Load experiment arrays once
    x, y, thick, hv, thre, npe = load_experiment_arrays()

    # Do both thickness/aerogel = 2 and 3
    for th in [2, 3]:
        # EXP map
        exp_map = make_exp_map(th, x, y, thick, hv, thre, npe)
        exp_map_cmp = apply_compare_mask(exp_map)

        # SIM map
        sim_mean_map, sim_err_map = make_sim_map(th)
        sim_mean_cmp = apply_compare_mask(sim_mean_map)

        # Plots: exp/sim
        plot_heatmap(
            exp_map,
            title=f"Data Np.e. map (HV={HV_TARGET}, thick={th}, thre={THRE_BY_THICK[th]})",
            outname=os.path.join(OUTDIR, f"data_npe_map_thick{th}.png"),
            cbar_label="Np.e.", fmt="{:.2f}"
        )
        plot_heatmap(
            sim_mean_map,
            title=f"Simulation Np.e. map (aerogel={th}) Gaussian mean",
            outname=os.path.join(OUTDIR, f"sim_npe_mean_map_aerogel{th}.png"),
            cbar_label="Np.e.", fmt="{:.2f}"
        )
        plot_heatmap(
            sim_err_map,
            title=f"Simulation Np.e. map (aerogel={th}) Gaussian mean error",
            outname=os.path.join(OUTDIR, f"sim_npe_meanerr_map_aerogel{th}.png"),
            cbar_label="Np.e. error", fmt="{:.3f}"
        )

        # Compare maps (exclude y=±54 by default)
        rel = make_relative_diff(exp_map_cmp, sim_mean_cmp)
        diff = make_absolute_diff(exp_map_cmp, sim_mean_cmp)

        # relative difference heatmap (percent)
        plot_heatmap(
            100.0 * rel,
            title=f"Relative difference (Data−Sim)/Data [%] (exclude y=±54), aerogel={th}",
            outname=os.path.join(OUTDIR, f"reldiff_percent_aerogel{th}_excludeY54.png"),
            cbar_label="[%]", fmt="{:.1f}",
            vmin=-30, vmax=30  # adjust as you like
        )

        # absolute difference heatmap
        plot_heatmap(
            diff,
            title=f"Absolute difference Data−Sim [p.e.] (exclude y=±54), aerogel={th}",
            outname=os.path.join(OUTDIR, f"absdiff_pe_aerogel{th}_excludeY54.png"),
            cbar_label="[p.e.]", fmt="{:.2f}",
            vmin=-10, vmax=10  # adjust as you like
        )

    print("DONE.")


if __name__ == "__main__":
    main()
