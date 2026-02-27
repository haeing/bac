#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

# EXP
import ROOT
ROOT.gROOT.SetBatch(True)

# SIM
import uproot
from scipy.optimize import curve_fit

# =========================
# CONFIG
# =========================
OUTDIR = "plots_compare"
os.makedirs(OUTDIR, exist_ok=True)

# Experiment ROOT
EXP_FILE = "kek_result.root"
EXP_TREE = "tree"
BR_X       = "x_pos"
BR_Y       = "y_pos"
BR_THICK   = "bac_thick"
BR_HV      = "bac_HV"
BR_THRE    = "bac_thre"
BR_NPE     = "bac_npe"
BR_NPE_ERR = "bac_npe_err"   # <-- 추가

HV_TARGET = 58
THRE_BY_THICK = {2: 60, 3: 75}

# Simulation ROOT file pattern (aerogel2/3)
SIM_PATTERN = "../../../data/KEKAR2023Dec/g4_root/kek_aerogel*_x*_y*.root"
SIM_TREE = "tree"
SIM_BR   = "nhMppc"

# Fixed grid (mm)
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=int)
ys = np.array([-54, -36, -18, 0, 18, 36, 54], dtype=int)

# Exclude outside-aerogel rows
EXCLUDE_Y = set([54, -54])

# Gaussian-fit settings for simulation
NBINS = 120
FIT_WINDOW_SIGMA = 2.0

# Histogram settings
REL_BINS = 30
REL_RANGE = (-0.30, 0.30)   # (Data-Sim)/Data
PULL_BINS = 30
PULL_RANGE = (-5.0, 5.0)

# =========================
# Utilities
# =========================
def parse_xy_aerogel(fname: str):
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

def apply_exclude_y_mask(m):
    m2 = np.array(m, copy=True)
    for iy, yy in enumerate(ys):
        if yy in EXCLUDE_Y:
            m2[iy, :] = np.nan
    return m2

# =========================
# Experiment maps (mean + error)
# =========================
def load_experiment_arrays():
    f = ROOT.TFile.Open(EXP_FILE)
    t = f.Get(EXP_TREE)
    if not t:
        raise RuntimeError(f"Cannot find TTree '{EXP_TREE}' in {EXP_FILE}")

    x, y, thick, hv, thre, npe, npe_err = [], [], [], [], [], [], []
    for ev in t:
        x.append(int(getattr(ev, BR_X)))
        y.append(int(getattr(ev, BR_Y)))
        thick.append(int(getattr(ev, BR_THICK)))
        hv.append(int(getattr(ev, BR_HV)))
        thre.append(int(getattr(ev, BR_THRE)))
        npe.append(float(getattr(ev, BR_NPE)))
        npe_err.append(float(getattr(ev, BR_NPE_ERR)))

    return (np.array(x, dtype=int),
            np.array(y, dtype=int),
            np.array(thick, dtype=int),
            np.array(hv, dtype=int),
            np.array(thre, dtype=int),
            np.array(npe, dtype=float),
            np.array(npe_err, dtype=float))

def make_exp_map(th, x, y, thick, hv, thre, npe, npe_err):
    thre_target = THRE_BY_THICK[th]
    sel = (hv == HV_TARGET) & (thick == th) & (thre == thre_target)

    m = np.full((len(ys), len(xs)), np.nan, dtype=float)
    e = np.full((len(ys), len(xs)), np.nan, dtype=float)

    for iy, yy in enumerate(ys):
        for ix, xx in enumerate(xs):
            mask = sel & (x == xx) & (y == yy)
            if np.any(mask):
                vals = npe[mask]
                errs = npe_err[mask]

                m[iy, ix] = float(np.mean(vals))

                # 평균의 오차: sqrt(sum(err_i^2))/N (독립 오차 가정)
                N = float(np.sum(mask))
                e[iy, ix] = float(np.sqrt(np.sum(errs**2)) / N) if N > 0 else np.nan

    return m, e

# =========================
# Simulation maps (mean + fit error)
# =========================
def make_sim_map(aerogel_id):
    files = sorted(glob.glob(SIM_PATTERN))
    if not files:
        raise RuntimeError(f"No files matched: {SIM_PATTERN}")

    mean_map = np.full((len(ys), len(xs)), np.nan, dtype=float)
    err_map  = np.full((len(ys), len(xs)), np.nan, dtype=float)

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

    return mean_map, err_map

# =========================
# Vectors: relative diff / pull
# =========================
def rel_diff_vector(data_map, sim_map):
    dm = apply_exclude_y_mask(data_map)
    sm = apply_exclude_y_mask(sim_map)
    good = np.isfinite(dm) & np.isfinite(sm) & (dm != 0.0)
    return (dm[good] - sm[good]) / dm[good]

def pull_vector(data_map, data_err_map, sim_map, sim_err_map):
    dm = apply_exclude_y_mask(data_map)
    de = apply_exclude_y_mask(data_err_map)
    sm = apply_exclude_y_mask(sim_map)
    se = apply_exclude_y_mask(sim_err_map)

    good = np.isfinite(dm) & np.isfinite(sm) & np.isfinite(de) & np.isfinite(se)
    denom = np.sqrt(de[good]**2 + se[good]**2)
    good2 = denom > 0
    return (dm[good][good2] - sm[good][good2]) / denom[good2]

def summarize(v):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {}
    return {
        "N": int(v.size),
        "mean": float(np.mean(v)),
        "std": float(np.std(v, ddof=1)) if v.size > 1 else np.nan,
        "median": float(np.median(v)),
        "p16": float(np.percentile(v, 16)),
        "p84": float(np.percentile(v, 84)),
    }

# =========================
# Plot helpers
# =========================
def plot_hist(v, title, xlabel, outname, bins=30, vrange=None, note=None):
    s = summarize(v)
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.hist(v, bins=bins, range=vrange, histtype="stepfilled", alpha=0.7)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")

    if s:
        txt = (
            f"N={s['N']}\n"
            f"mean={s['mean']:.3f}\n"
            f"median={s['median']:.3f}\n"
            f"std={s['std']:.3f}\n"
            f"16–84%: [{s['p16']:.3f}, {s['p84']:.3f}]"
        )
        if note:
            txt = note + "\n" + txt
        ax.text(0.98, 0.98, txt, transform=ax.transAxes,
                ha="right", va="top", fontsize=10,
                bbox=dict(boxstyle="round", alpha=0.15))

    fig.tight_layout()
    fig.savefig(outname, dpi=200)
    plt.close(fig)
    print(f"Saved: {outname}")

def plot_hist_overlay(v2, v3, title, xlabel, outname, bins=30, vrange=None, labels=("aerogel=2","aerogel=3")):
    s2 = summarize(v2)
    s3 = summarize(v3)

    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ax.hist(v2, bins=bins, range=vrange, histtype="step", linewidth=2, label=labels[0])
    ax.hist(v3, bins=bins, range=vrange, histtype="step", linewidth=2, label=labels[1])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")
    ax.legend()



    fig.tight_layout()
    fig.savefig(outname, dpi=200)
    plt.close(fig)
    print(f"Saved: {outname}")

# =========================
# Main
# =========================
def main():
    x, y, thick, hv, thre, npe, npe_err = load_experiment_arrays()

    # aerogel=2
    exp2, exp2e = make_exp_map(2, x, y, thick, hv, thre, npe, npe_err)
    sim2, sim2e = make_sim_map(2)
    rel2  = rel_diff_vector(exp2, sim2)
    pull2 = pull_vector(exp2, exp2e, sim2, sim2e)

    # aerogel=3
    exp3, exp3e = make_exp_map(3, x, y, thick, hv, thre, npe, npe_err)
    sim3, sim3e = make_sim_map(3)
    rel3  = rel_diff_vector(exp3, sim3)
    pull3 = pull_vector(exp3, exp3e, sim3, sim3e)

    # --- Relative diff hists
    plot_hist(rel2, "Relative difference (exclude y=±54)", r"$(\mathrm{Data}-\mathrm{Sim})/\mathrm{Data}$",
              os.path.join(OUTDIR, "reldiff_hist_aerogel2.png"),
              bins=REL_BINS, vrange=REL_RANGE, note="aerogel=2")

    plot_hist(rel3, "Relative difference (exclude y=±54)", r"$(\mathrm{Data}-\mathrm{Sim})/\mathrm{Data}$",
              os.path.join(OUTDIR, "reldiff_hist_aerogel3.png"),
              bins=REL_BINS, vrange=REL_RANGE, note="aerogel=3")

    plot_hist_overlay(rel2, rel3, "",
                      r"Relative difference in $N_{\mathrm{pe}}$",
                      os.path.join(OUTDIR, "reldiff_hist_overlay_a2_a3.png"),
                      bins=REL_BINS, vrange=REL_RANGE)

    # --- Pull hists
    plot_hist(pull2, "Pull distribution (exclude y=±54)", r"$(\mathrm{Data}-\mathrm{Sim})/\sqrt{\sigma_{\rm data}^2+\sigma_{\rm sim}^2}$",
              os.path.join(OUTDIR, "pull_hist_aerogel2.png"),
              bins=PULL_BINS, vrange=PULL_RANGE, note="aerogel=2")

    plot_hist(pull3, "Pull distribution (exclude y=±54)", r"$(\mathrm{Data}-\mathrm{Sim})/\sqrt{\sigma_{\rm data}^2+\sigma_{\rm sim}^2}$",
              os.path.join(OUTDIR, "pull_hist_aerogel3.png"),
              bins=PULL_BINS, vrange=PULL_RANGE, note="aerogel=3")

    plot_hist_overlay(pull2, pull3, "Pull distribution (exclude y=±54)",
                      r"$(\mathrm{Data}-\mathrm{Sim})/\sqrt{\sigma_{\rm data}^2+\sigma_{\rm sim}^2}$",
                      os.path.join(OUTDIR, "pull_hist_overlay_a2_a3.png"),
                      bins=PULL_BINS, vrange=PULL_RANGE)

    print("\n=== Summary (exclude y=±54) ===")
    print("rel aerogel=2:", summarize(rel2))
    print("rel aerogel=3:", summarize(rel3))
    print("pull aerogel=2:", summarize(pull2))
    print("pull aerogel=3:", summarize(pull3))
    print("DONE.")

if __name__ == "__main__":
    main()
