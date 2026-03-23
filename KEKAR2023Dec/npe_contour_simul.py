import ROOT
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

ROOT.gROOT.SetBatch(True)

# =========================================================
# INPUT
# =========================================================
SIM_FILE = "../../../data/KEKAR2023Dec/g4_root/kek_aerogel3_x0_y0_total.root"
SIM_TREE = "tree"              # <- 바꿔줘

BR_X   = "evtposx"
BR_Y   = "evtposy"
BR_NPE = "nhMppc"

# =========================================================
# Detector geometry
# =========================================================
DET_SIZE = 115.0
DET_HALF = DET_SIZE / 2.0
DET_XMIN, DET_XMAX = -DET_HALF, +DET_HALF
DET_YMIN, DET_YMAX = -DET_HALF, +DET_HALF

# =========================================================
# Plot / bin settings
# =========================================================
BIN_MM = 2.0            # 2 mm binning
SMOOTH_SIGMA = 1.2      # smoothing in bin unit
NLEVELS = 16
CMAP = "viridis"

# =========================================================
# Read simulation tree
# =========================================================
def read_simulation(filename, tree_name):
    rf = ROOT.TFile.Open(filename)
    if not rf or rf.IsZombie():
        raise RuntimeError(f"Cannot open ROOT file: {filename}")

    tt = rf.Get(tree_name)
    if not tt:
        raise RuntimeError(f"Cannot find TTree '{tree_name}' in {filename}")

    xs, ys, zs = [], [], []

    for ev in tt:
        xs.append(float(getattr(ev, BR_X)))
        ys.append(float(getattr(ev, BR_Y)))
        zs.append(float(getattr(ev, BR_NPE)))

    xs = np.array(xs, dtype=float)
    ys = np.array(ys, dtype=float)
    zs = np.array(zs, dtype=float)

    good = np.isfinite(xs) & np.isfinite(ys) & np.isfinite(zs)
    xs = xs[good]
    ys = ys[good]
    zs = zs[good]

    print(f"[SIM] loaded events = {len(xs)}")
    return xs, ys, zs


# =========================================================
# Build 2D mean map of nhMppc
# =========================================================
def build_mean_response_map(x, y, z, bin_mm=2.0):
    x_edges = np.arange(DET_XMIN, DET_XMAX + bin_mm, bin_mm)
    y_edges = np.arange(DET_YMIN, DET_YMAX + bin_mm, bin_mm)

    # weighted sum
    sum_map, _, _ = np.histogram2d(
        x, y,
        bins=[x_edges, y_edges],
        weights=z
    )

    # counts
    cnt_map, _, _ = np.histogram2d(
        x, y,
        bins=[x_edges, y_edges]
    )

    # histogram2d returns [xbin, ybin], transpose for plotting
    sum_map = sum_map.T
    cnt_map = cnt_map.T

    mean_map = np.full_like(sum_map, np.nan, dtype=float)
    mask = cnt_map > 0
    mean_map[mask] = sum_map[mask] / cnt_map[mask]

    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    return x_centers, y_centers, mean_map, cnt_map


# =========================================================
# Smooth map while handling empty bins
# =========================================================
def smooth_nan_map(z, sigma=1.2):
    """
    Smooth a map containing NaNs by filtering value and weight separately.
    """
    z0 = np.where(np.isfinite(z), z, 0.0)
    w0 = np.where(np.isfinite(z), 1.0, 0.0)

    z_s = gaussian_filter(z0, sigma=sigma, mode="nearest")
    w_s = gaussian_filter(w0, sigma=sigma, mode="nearest")

    out = np.full_like(z, np.nan, dtype=float)
    good = w_s > 1e-8
    out[good] = z_s[good] / w_s[good]
    return out


# =========================================================
# Plot contour
# =========================================================
def plot_simulation_contour(xc, yc, zmap, outname):
    X, Y = np.meshgrid(xc, yc)

    fig, ax = plt.subplots(figsize=(7.6, 6.6))

    finite = np.isfinite(zmap)
    if not np.any(finite):
        raise RuntimeError("No valid bins to plot.")

    #vmin = np.nanmin(zmap)
    #vmax = np.nanmax(zmap)
    vmin = 0
    vmax = 50

    if np.isclose(vmin, vmax):
        vmax = vmin + 1e-6

    levels = np.linspace(vmin, vmax, NLEVELS)

    cf = ax.contourf(X, Y, zmap, levels=levels, cmap=CMAP)
    ax.contour(X, Y, zmap, levels=levels, colors="k", linewidths=0.25, alpha=0.22)

    # detector boundary
    ax.plot(
        [DET_XMIN, DET_XMAX, DET_XMAX, DET_XMIN, DET_XMIN],
        [DET_YMIN, DET_YMIN, DET_YMAX, DET_YMAX, DET_YMIN],
        color="black", linewidth=1.5
    )

    ax.set_xlabel("X (mm)", fontsize=18)
    ax.set_ylabel("Y (mm)", fontsize=18)
    ax.tick_params(axis="both", labelsize=14)

    ax.set_xlim(DET_XMIN, DET_XMAX)
    ax.set_ylim(DET_YMIN, DET_YMAX)
    ax.set_aspect("equal", adjustable="box")


    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("Np.e.", fontsize=18)
    cbar.set_ticks(np.arange(0, 51, 10))
    cbar.ax.tick_params(labelsize=14)

    fig.tight_layout()
    fig.savefig(outname, dpi=220)
    plt.close(fig)

    print(f"Saved: {outname}")


# =========================================================
# MAIN
# =========================================================
if __name__ == "__main__":
    xsim, ysim, nh = read_simulation(SIM_FILE, SIM_TREE)

    xc, yc, zmean, zcount = build_mean_response_map(
        xsim, ysim, nh,
        bin_mm=BIN_MM
    )

    zsmooth = smooth_nan_map(zmean, sigma=SMOOTH_SIGMA)

    plot_simulation_contour(
        xc, yc, zsmooth,
        outname="simulation_mean_nhMppc_contour_3.png"
    )

    print("DONE.")
