import ROOT
import glob, re
import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

ROOT.gROOT.SetBatch(True)

# =========================
# Common grid (mm)
# =========================
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=int)
ys = np.array([-54, -36, -18, 0, 18, 36, 54], dtype=int)

# =========================
# Experiment ROOT settings
# =========================
EXP_FILE = "kek_result.root"
EXP_TREE = "tree"
BR_X     = "x_pos"
BR_Y     = "y_pos"
BR_THICK = "bac_thick"
BR_HV    = "bac_HV"
BR_THRE  = "bac_thre"
BR_NPE   = "bac_npe"
BR_NPEE  = "bac_npe_err"

HV_TARGET = 58
THRE_BY_THICK = {2: 60, 3: 75}

# =========================
# Simulation ROOT settings
# =========================
SIM_PATTERN = "../../../data/KEKAR2023Dec/g4_root/kek_aerogel*_x*_y*.root"
SIM_TREE = "tree"
SIM_BR   = "nhMppc"

NBINS = 120
FIT_WINDOW_SIGMA = 2.0

# -------------------------
# Utils
# -------------------------
def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def fit_gaussian_mean(data):
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if len(data) < 30:
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
    if len(xfit) < 6:
        return np.nan, np.nan

    yerr = np.sqrt(yfit)
    yerr[yerr == 0] = 1.0

    A0 = float(np.max(yfit))

    try:
        popt, pcov = curve_fit(
            gauss, xfit, yfit,
            p0=[A0, mu0, sigma0],
            sigma=yerr, absolute_sigma=True,
            maxfev=10000
        )
        mu = float(popt[1])
        mu_err = float(np.sqrt(pcov[1, 1])) if pcov is not None and pcov.shape == (3, 3) else np.nan
        return mu, mu_err
    except Exception:
        return np.nan, np.nan

def parse_sim_file(fn):
    """
    kek_aerogel3_xm16_ym54.root -> (thick=3, x=-16, y=-54)
    kek_aerogel2_x16_y54.root   -> (thick=2, x=16, y=54)
    """
    m = re.search(r"kek_aerogel(?P<th>\d+)_x(?P<x>m?\d+)_y(?P<y>m?\d+)\.root$", fn)
    if not m:
        raise ValueError(f"Cannot parse: {fn}")

    th = int(m.group("th"))

    def tok2int(tok):
        return -int(tok[1:]) if tok.startswith("m") else int(tok)

    x = tok2int(m.group("x"))
    y = tok2int(m.group("y"))
    return th, x, y

# -------------------------
# Load experiment as dict: exp[(th, x, y)] = (mean, err)
# -------------------------
def load_experiment():
    f = ROOT.TFile.Open(EXP_FILE)
    t = f.Get(EXP_TREE)
    if not t:
        raise RuntimeError(f"Cannot find TTree '{EXP_TREE}' in {EXP_FILE}")

    exp = {}
    for th in [2, 3]:
        thre_target = THRE_BY_THICK[th]
        # collect per (x,y) -> list of (npe, npe_err)
        buf = {}
        for ev in t:
            if int(getattr(ev, BR_HV)) != HV_TARGET:   continue
            if int(getattr(ev, BR_THICK)) != th:       continue
            if int(getattr(ev, BR_THRE)) != thre_target: continue

            xx = int(getattr(ev, BR_X))
            yy = int(getattr(ev, BR_Y))
            val = float(getattr(ev, BR_NPE))
            err = float(getattr(ev, BR_NPEE)) if hasattr(ev, BR_NPEE) else np.nan

            buf.setdefault((xx, yy), []).append((val, err))

        # average per bin
        for (xx, yy), arr in buf.items():
            arr = np.asarray(arr, dtype=float)
            m = float(np.nanmean(arr[:, 0]))
            # if bac_npe_err is already "error of mean" per entry, averaging isn't strictly correct,
            # but this is a pragmatic choice: take RMS of provided errs / sqrt(N) if all finite.
            e = arr[:, 1]
            if np.all(np.isfinite(e)) and len(e) > 0:
                e_mean = float(np.sqrt(np.sum(e**2)) / len(e))  # conservative-ish
            else:
                e_mean = np.nan
            exp[(th, xx, yy)] = (m, e_mean)

    return exp

# -------------------------
# Load simulation as dict: sim[(th, x, y)] = (mu, mu_err)
# -------------------------
def load_simulation():
    sim = {}
    files = sorted(glob.glob(SIM_PATTERN))
    if not files:
        raise RuntimeError(f"No sim files matched: {SIM_PATTERN}")

    for fn in files:
        th, x, y = parse_sim_file(fn)
        if th not in (2, 3): 
            continue
        if x not in xs or y not in ys:
            continue

        with uproot.open(fn) as f:
            data = f[SIM_TREE][SIM_BR].array(library="np")

        mu, mu_err = fit_gaussian_mean(data)
        sim[(th, x, y)] = (mu, mu_err)

    return sim

# -------------------------
# Build 1D lines
# -------------------------
def line_y0(dct, th, use_err=True):
    # y=0, x varies, exclude y=±54 not relevant here anyway
    X = []
    Y = []
    E = []
    for x in xs:
        key = (th, int(x), 0)
        if key in dct:
            yv, ev = dct[key]
            X.append(x)
            Y.append(yv)
            E.append(ev if use_err else np.nan)
    return np.array(X), np.array(Y), np.array(E)

def line_x0(dct, th, use_err=True):
    # x=0, y varies, exclude ±54
    YPOS = []
    VAL = []
    ERR = []
    for y in ys:
        if y in (-54, 54):
            continue
        key = (th, 0, int(y))
        if key in dct:
            v, e = dct[key]
            YPOS.append(y)
            VAL.append(v)
            ERR.append(e if use_err else np.nan)
    return np.array(YPOS), np.array(VAL), np.array(ERR)

# -------------------------
# Plot helpers
# -------------------------
def plot_overlay_y0(exp, sim, outname="compare_y0.png"):
    fig, ax = plt.subplots(figsize=(7.2, 4.8))

    xerr = 9.0/2.0  # mm
    for th in [2, 3]:
        X, Y, E = line_y0(exp, th, use_err=True)
        ax.errorbar(X, Y, yerr=E, xerr=xerr, fmt='o', capsize=2,
                    label=f"Data aerogel{th} (y=0)")
        
        #Xs, Ys, Es = line_y0(sim, th, use_err=True)
        #ax.errorbar(Xs, Ys, yerr=Es, xerr=xerr, fmt='s', capsize=2,
        #           label=f"Sim aerogel{th} (y=0)")

    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Np.e.")
    ax.set_title("Data vs Simulation (y = 0 mm)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outname, dpi=200)
    plt.close(fig)

def plot_overlay_x0(exp, sim, outname="compare_x0.png"):
    fig, ax = plt.subplots(figsize=(7.2, 4.8))

    xerr = 13.0/2.0  # mm (trigger height in y)
    for th in [2, 3]:
        YPOS, VAL, ERR = line_x0(exp, th, use_err=True)
        ax.errorbar(YPOS, VAL, yerr=ERR, xerr=xerr, fmt='o', capsize=2,
                    label=f"Data aerogel{th} (x=0)")

        #YPOSs, VALs, ERRs = line_x0(sim, th, use_err=True)
        #ax.errorbar(YPOSs, VALs, yerr=ERRs, xerr=xerr, fmt='s', capsize=2,
        #           label=f"Sim aerogel{th} (x=0)")

    ax.set_xlabel("Y (mm)")
    ax.set_ylabel("Np.e.")
    ax.set_title("Data vs Simulation (x = 0 mm)  (excluding y = ±54 mm)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outname, dpi=200)
    plt.close(fig)

# =========================
# Main
# =========================
exp = load_experiment()
sim = load_simulation()

plot_overlay_y0(exp, sim, "compare_1D_y0.png")
plot_overlay_x0(exp, sim, "compare_1D_x0.png")

print("Saved: compare_1D_y0.png, compare_1D_x0.png")
