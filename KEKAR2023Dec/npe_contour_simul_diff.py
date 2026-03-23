import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.ndimage import gaussian_filter

ROOT.gROOT.SetBatch(True)

# =========================================================
# INPUT: measurement data ROOT
# =========================================================
DATA_FILE = "kek_result.root"
DATA_TREE = "tree"

BR_X        = "x_pos"
BR_Y        = "y_pos"
BR_THICK    = "bac_thick"
BR_HV       = "bac_HV"
BR_THRE     = "bac_thre"
BR_NPE      = "bac_npe"
BR_NPE_ERR  = "bac_npe_err"   # optional

# =========================================================
# INPUT: simulation ROOT
# thickness-dependent file
# use {th} placeholder
# =========================================================
SIM_FILE_TEMPLATE = "../../../data/KEKAR2023Dec/g4_root/kek_aerogel{th}_x0_y0_total.root"
SIM_TREE = "tree"

SIM_BR_X   = "evtposx"
SIM_BR_Y   = "evtposy"
SIM_BR_NPE = "nhMppc"

# =========================================================
# Detector / scan geometry
# =========================================================
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=float)
ys = np.array([-36, -18, 0, 18, 36], dtype=float)

TRIG_WX = 9.0
TRIG_WY = 13.0

DET_SIZE = 115.0
DET_HALF = DET_SIZE / 2.0
DET_XMIN, DET_XMAX = -DET_HALF, +DET_HALF
DET_YMIN, DET_YMAX = -DET_HALF, +DET_HALF

# =========================================================
# Selection
# =========================================================
HV_TARGET = 58
THRE_BY_THICK = {2: 60, 3: 75}

# =========================================================
# Plot configuration
# =========================================================
BIN_MM = 2.0
SMOOTH_SIGMA = 1.2        # in bin unit for simulation contour
CMAP = "viridis"
ZMIN = 0
ZMAX = 50
NLEVELS = 15

FIGSIZE = (7.8, 6.8)

SHOW_BOX_TEXT = True
BOX_TEXT_FONTSIZE = 8
BOX_EDGE_COLOR = "white"
BOX_EDGE_ALPHA = 1
USE_TEXT_BBOX = True

# =========================================================
# Helpers
# =========================================================
def rect_edges_from_center(xc, yc, wx, wy):
    return xc - wx / 2.0, xc + wx / 2.0, yc - wy / 2.0, yc + wy / 2.0


def get_sim_file(th):
    return SIM_FILE_TEMPLATE.format(th=th)


# =========================================================
# Measurement data
# =========================================================
def read_measurement_tree():
    rf = ROOT.TFile.Open(DATA_FILE)
    if not rf or rf.IsZombie():
        raise RuntimeError(f"Cannot open data ROOT file: {DATA_FILE}")

    tt = rf.Get(DATA_TREE)
    if not tt:
        raise RuntimeError(f"Cannot find TTree '{DATA_TREE}' in {DATA_FILE}")

    has_npe_err = (tt.GetBranch(BR_NPE_ERR) is not None)

    x, y, thick, hv, thre, npe = [], [], [], [], [], []
    npe_err = []

    for ev in tt:
        x.append(float(getattr(ev, BR_X)))
        y.append(float(getattr(ev, BR_Y)))
        thick.append(int(getattr(ev, BR_THICK)))
        hv.append(int(getattr(ev, BR_HV)))
        thre.append(int(getattr(ev, BR_THRE)))
        npe.append(float(getattr(ev, BR_NPE)))
        if has_npe_err:
            npe_err.append(float(getattr(ev, BR_NPE_ERR)))

    out = {
        "x": np.array(x, dtype=float),
        "y": np.array(y, dtype=float),
        "thick": np.array(thick, dtype=int),
        "hv": np.array(hv, dtype=int),
        "thre": np.array(thre, dtype=int),
        "npe": np.array(npe, dtype=float),
        "has_npe_err": has_npe_err,
    }

    if has_npe_err:
        out["npe_err"] = np.array(npe_err, dtype=float)

    print(f"[DATA] loaded entries = {len(out['x'])}")
    return out


def build_measurement_maps(data_dict, th):
    thre_target = THRE_BY_THICK[th]

    x = data_dict["x"]
    y = data_dict["y"]
    thick = data_dict["thick"]
    hv = data_dict["hv"]
    thre = data_dict["thre"]
    npe = data_dict["npe"]

    sel = (hv == HV_TARGET) & (thick == th) & (thre == thre_target)
    print(f"[DATA] thick={th}, thre={thre_target}, selected = {int(np.sum(sel))}")

    m_val = np.full((len(ys), len(xs)), np.nan, dtype=float)
    m_err = np.full((len(ys), len(xs)), np.nan, dtype=float)
    m_n   = np.zeros((len(ys), len(xs)), dtype=int)

    for iy, yy in enumerate(ys):
        for ix, xx in enumerate(xs):
            mask = sel & (x == xx) & (y == yy)
            n = int(np.sum(mask))
            if n <= 0:
                continue

            vals = npe[mask]
            m_val[iy, ix] = np.mean(vals)
            m_n[iy, ix] = n

            if data_dict["has_npe_err"]:
                ee = data_dict["npe_err"][mask]
                m_err[iy, ix] = np.sqrt(np.sum(ee * ee)) / n
            else:
                if n > 1:
                    m_err[iy, ix] = np.std(vals, ddof=1) / np.sqrt(n)
                else:
                    m_err[iy, ix] = 0.0

    return m_val, m_err, thre_target


# =========================================================
# Simulation
# =========================================================
def read_simulation_tree(sim_file):
    rf = ROOT.TFile.Open(sim_file)
    if not rf or rf.IsZombie():
        raise RuntimeError(f"Cannot open simulation ROOT file: {sim_file}")

    tt = rf.Get(SIM_TREE)
    if not tt:
        raise RuntimeError(f"Cannot find TTree '{SIM_TREE}' in {sim_file}")

    xsim, ysim, zsim = [], [], []

    for ev in tt:
        xsim.append(float(getattr(ev, SIM_BR_X)))
        ysim.append(float(getattr(ev, SIM_BR_Y)))
        zsim.append(float(getattr(ev, SIM_BR_NPE)))

    xsim = np.array(xsim, dtype=float)
    ysim = np.array(ysim, dtype=float)
    zsim = np.array(zsim, dtype=float)

    good = np.isfinite(xsim) & np.isfinite(ysim) & np.isfinite(zsim)
    xsim = xsim[good]
    ysim = ysim[good]
    zsim = zsim[good]

    print(f"[SIM] file = {sim_file}")
    print(f"[SIM] loaded events = {len(xsim)}")
    return xsim, ysim, zsim


def build_sim_mean_map(xsim, ysim, zsim, bin_mm=2.0):
    x_edges = np.arange(DET_XMIN, DET_XMAX + bin_mm, bin_mm)
    y_edges = np.arange(DET_YMIN, DET_YMAX + bin_mm, bin_mm)

    sum_map, _, _ = np.histogram2d(
        xsim, ysim,
        bins=[x_edges, y_edges],
        weights=zsim
    )
    cnt_map, _, _ = np.histogram2d(
        xsim, ysim,
        bins=[x_edges, y_edges]
    )

    sum_map = sum_map.T
    cnt_map = cnt_map.T

    mean_map = np.full_like(sum_map, np.nan, dtype=float)
    mask = cnt_map > 0
    mean_map[mask] = sum_map[mask] / cnt_map[mask]

    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    return x_centers, y_centers, mean_map


def smooth_nan_map(z, sigma=1.2):
    z0 = np.where(np.isfinite(z), z, 0.0)
    w0 = np.where(np.isfinite(z), 1.0, 0.0)

    z_s = gaussian_filter(z0, sigma=sigma, mode="nearest")
    w_s = gaussian_filter(w0, sigma=sigma, mode="nearest")

    out = np.full_like(z, np.nan, dtype=float)
    good = w_s > 1e-8
    out[good] = z_s[good] / w_s[good]
    return out


def calc_sim_in_rectangle(xsim, ysim, zsim, xc, yc, wx, wy):
    x0, x1, y0, y1 = rect_edges_from_center(xc, yc, wx, wy)
    mask = (xsim >= x0) & (xsim < x1) & (ysim >= y0) & (ysim < y1)

    n = int(np.sum(mask))
    if n <= 0:
        return np.nan, np.nan, 0

    vals = zsim[mask]
    mean = np.mean(vals)
    if n > 1:
        err = np.std(vals, ddof=1) / np.sqrt(n)
    else:
        err = 0.0
    return mean, err, n


# =========================================================
# Relative difference
# =========================================================
def calc_relative_difference_percent(d, sd, s, ss):
    """
    delta = 100 * (d - s) / d
    """
    if not np.isfinite(d) or not np.isfinite(s) or d == 0:
        return np.nan, np.nan

    delta = 100.0 * (d - s) / d
    sigma = 100.0 * np.sqrt((s * sd / (d * d)) ** 2 + (ss / d) ** 2)
    return delta, sigma


def add_measurement_boxes_with_diff(ax, m_data, m_data_err, xsim, ysim, zsim):
    for iy, yc in enumerate(ys):
        for ix, xc in enumerate(xs):
            d = m_data[iy, ix]
            sd = m_data_err[iy, ix]

            if not np.isfinite(d):
                continue

            sim_mean, sim_err, n_sim = calc_sim_in_rectangle(
                xsim, ysim, zsim, xc, yc, TRIG_WX, TRIG_WY
            )

            x0, x1, y0, y1 = rect_edges_from_center(xc, yc, TRIG_WX, TRIG_WY)

            ax.add_patch(
                Rectangle(
                    (x0, y0), x1 - x0, y1 - y0,
                    fill=False,
                    edgecolor=BOX_EDGE_COLOR,
                    linewidth=0.8,
                    alpha=BOX_EDGE_ALPHA
                )
            )

            if not SHOW_BOX_TEXT:
                continue

            if not np.isfinite(sim_mean):
                txt = "N/A"
            else:
                delta, sigma = calc_relative_difference_percent(d, sd, sim_mean, sim_err)
                if np.isfinite(delta) and np.isfinite(sigma):
                    txt = f"{delta:+.0f}\n±{sigma:.0f}%"
                else:
                    txt = "N/A"

            bbox = dict(boxstyle="round,pad=0.15", fc=(1, 1, 1, 0.45), ec="none") if USE_TEXT_BBOX else None

            ax.text(
                xc, yc, txt,
                ha="center", va="center",
                fontsize=BOX_TEXT_FONTSIZE,
                color="black",
                bbox=bbox
            )


# =========================================================
# Plot
# =========================================================
def plot_simulation_contour_with_measurement_diff(
    xc, yc, zmap, m_data, m_data_err, xsim, ysim, zsim, outname, th, thre_target
):
    X, Y = np.meshgrid(xc, yc)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    levels = np.linspace(ZMIN, ZMAX, NLEVELS)
    cf = ax.contourf(X, Y, zmap, levels=levels, cmap=CMAP, vmin=ZMIN, vmax=ZMAX)
    ax.contour(X, Y, zmap, levels=levels, colors="k", linewidths=0.22, alpha=0.18)

    ax.plot(
        [DET_XMIN, DET_XMAX, DET_XMAX, DET_XMIN, DET_XMIN],
        [DET_YMIN, DET_YMIN, DET_YMAX, DET_YMAX, DET_YMIN],
        color="black", linewidth=1.5
    )

    add_measurement_boxes_with_diff(ax, m_data, m_data_err, xsim, ysim, zsim)

    ax.set_xlim(DET_XMIN, DET_XMAX)
    ax.set_ylim(DET_YMIN, DET_YMAX)

    ax.set_xticks([-40, -20, 0, 20, 40])
    ax.set_yticks([-40, -20, 0, 20, 40])

    ax.set_xlabel("X (mm)", fontsize=18)
    ax.set_ylabel("Y (mm)", fontsize=18)
    ax.tick_params(axis="both", labelsize=14)

    ax.set_aspect("equal", adjustable="box")

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("Np.e.", fontsize=18)
    cbar.set_ticks(np.arange(ZMIN, ZMAX + 1, 10))
    cbar.ax.tick_params(labelsize=14)

    fig.tight_layout()
    fig.savefig(outname, dpi=220)
    plt.close(fig)
    print(f"Saved: {outname}")


# =========================================================
# MAIN
# =========================================================
if __name__ == "__main__":
    data_dict = read_measurement_tree()

    for th in [2, 3]:
        m_data, m_data_err, thre_target = build_measurement_maps(data_dict, th)

        sim_file = get_sim_file(th)
        xsim, ysim, zsim = read_simulation_tree(sim_file)

        xc, yc, zmean = build_sim_mean_map(xsim, ysim, zsim, bin_mm=BIN_MM)
        zsmooth = smooth_nan_map(zmean, sigma=SMOOTH_SIGMA)

        outname = f"simulation_contour_with_data_diff_HV{HV_TARGET}_thick{th}_thre{thre_target}.png"

        plot_simulation_contour_with_measurement_diff(
            xc, yc, zsmooth,
            m_data, m_data_err,
            xsim, ysim, zsim,
            outname=outname,
            th=th,
            thre_target=thre_target
        )

    print("DONE.")
