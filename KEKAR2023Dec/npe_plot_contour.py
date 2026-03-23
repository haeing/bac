import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.ndimage import gaussian_filter

ROOT.gROOT.SetBatch(True)

# =========================
# INPUT ROOT
# =========================
FILE = "kek_result.root"
TREE = "tree"

BR_X     = "x_pos"
BR_Y     = "y_pos"
BR_THICK = "bac_thick"
BR_HV    = "bac_HV"
BR_THRE  = "bac_thre"
BR_NPE   = "bac_npe"
BR_NPE_ERR = "bac_npe_err"   # optional

# =========================
# GRID (mm)
# =========================
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=float)
ys = np.array([-36, -18, 0, 18, 36], dtype=float)

# =========================
# Selection rules
# =========================
HV_TARGET = 58
THRE_BY_THICK = {2: 60, 3: 75}

# =========================
# Trigger footprint (mm)
# =========================
TRIG_WX = 9.0
TRIG_WY = 13.0

# =========================
# Detector active area (mm)
# =========================
DET_SIZE = 115.0
DET_HALF = DET_SIZE / 2.0
DET_XMIN, DET_XMAX = -DET_HALF, +DET_HALF
DET_YMIN, DET_YMAX = -DET_HALF, +DET_HALF


# -------------------------
# Geometry helpers
# -------------------------
def rect_edges_from_center(xc, yc, wx, wy):
    return xc - wx / 2.0, xc + wx / 2.0, yc - wy / 2.0, yc + wy / 2.0


def overlap_1d(a0, a1, b0, b1):
    return max(0.0, min(a1, b1) - max(a0, b0))


# -------------------------
# Extrapolation helper
# -------------------------
def idw_fill_full_detector(z, x_centers, y_centers, xs_meas, ys_meas, m_meas,
                           power=2.0, eps=1e-6):
    """
    Fill NaN pixels using inverse-distance weighting (IDW)
    from measured center values.
    """
    z_filled = np.array(z, copy=True)

    valid_meas = np.isfinite(m_meas)
    meas_x = xs_meas[valid_meas]
    meas_y = ys_meas[valid_meas]
    meas_v = m_meas[valid_meas]

    if len(meas_v) == 0:
        return z_filled

    X, Y = np.meshgrid(x_centers, y_centers)
    miss = ~np.isfinite(z_filled)

    if not np.any(miss):
        return z_filled

    xm = X[miss]
    ym = Y[miss]

    # distance from each missing pixel to each measurement center
    dx = xm[:, None] - meas_x[None, :]
    dy = ym[:, None] - meas_y[None, :]
    d2 = dx * dx + dy * dy

    # exact hit protection
    exact = d2 < eps
    weights = 1.0 / np.maximum(d2, eps) ** (power / 2.0)

    # if exact match exists, use it directly
    vals = np.sum(weights * meas_v[None, :], axis=1) / np.sum(weights, axis=1)

    if np.any(exact):
        row_has_exact = np.any(exact, axis=1)
        exact_idx = np.argmax(exact, axis=1)
        vals[row_has_exact] = meas_v[exact_idx[row_has_exact]]

    z_filled[miss] = vals
    return z_filled


# -------------------------
# Read ROOT tree -> numpy arrays
# -------------------------
f = ROOT.TFile.Open(FILE)
t = f.Get(TREE)
if not t:
    raise RuntimeError(f"Cannot find TTree '{TREE}' in {FILE}")

x, y, thick, hv, thre, npe = [], [], [], [], [], []
npe_err = []

has_npe_err = (t.GetBranch(BR_NPE_ERR) is not None)

for ev in t:
    x.append(float(getattr(ev, BR_X)))
    y.append(float(getattr(ev, BR_Y)))
    thick.append(int(getattr(ev, BR_THICK)))
    hv.append(int(getattr(ev, BR_HV)))
    thre.append(int(getattr(ev, BR_THRE)))
    npe.append(float(getattr(ev, BR_NPE)))
    if has_npe_err:
        npe_err.append(float(getattr(ev, BR_NPE_ERR)))

x     = np.array(x, dtype=float)
y     = np.array(y, dtype=float)
thick = np.array(thick, dtype=int)
hv    = np.array(hv, dtype=int)
thre  = np.array(thre, dtype=int)
npe   = np.array(npe, dtype=float)
if has_npe_err:
    npe_err = np.array(npe_err, dtype=float)


# -------------------------
# Build map (mean Npe at each measurement center)
# -------------------------
def make_map(th):
    if th not in THRE_BY_THICK:
        raise ValueError(f"No threshold rule for thickness={th}")

    thre_target = THRE_BY_THICK[th]
    sel = (hv == HV_TARGET) & (thick == th) & (thre == thre_target)

    print(f"[thick={th}] select: HV=={HV_TARGET}, thre=={thre_target} -> {int(np.sum(sel))} entries")

    m = np.full((len(ys), len(xs)), np.nan, dtype=float)
    m_err = np.full((len(ys), len(xs)), np.nan, dtype=float)

    for iy, yy in enumerate(ys):
        for ix, xx in enumerate(xs):
            mask = sel & (x == xx) & (y == yy)
            if np.any(mask):
                m[iy, ix] = np.mean(npe[mask])

                if has_npe_err:
                    ee = npe_err[mask]
                    m_err[iy, ix] = np.sqrt(np.sum(ee * ee)) / np.sum(mask)

    return m, m_err, thre_target


# -------------------------
# Area-aware contour plot with full-detector extrapolation
# -------------------------
def plot_area_averaged_contour_map_full(
    m, th, thre_target, outname,
    pixel_mm=1.0,
    smooth_sigma_mm=3.0,
    nlevels=24,
    idw_power=2.0,
    draw_measurement_boxes=True,
    draw_centers=True,
    annotate_values=False
):
    """
    1) Treat each measured value as rectangle-average over 9x13 mm footprint
    2) Fill fine pixels by area overlap
    3) Extrapolate uncovered detector area with IDW from measurement centers
    4) Smooth and draw contour over the full 115x115 mm detector
    """

    m = np.asarray(m, dtype=float)

    # Fine grid edges
    x_edges = np.arange(DET_XMIN, DET_XMAX + pixel_mm, pixel_mm)
    y_edges = np.arange(DET_YMIN, DET_YMAX + pixel_mm, pixel_mm)

    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    nx = len(x_centers)
    ny = len(y_centers)

    sum_map = np.zeros((ny, nx), dtype=float)
    weight_map = np.zeros((ny, nx), dtype=float)

    # Fill measured rectangle regions
    for iy, yc in enumerate(ys):
        for ix, xc in enumerate(xs):
            val = m[iy, ix]
            if not np.isfinite(val):
                continue

            rx0, rx1, ry0, ry1 = rect_edges_from_center(xc, yc, TRIG_WX, TRIG_WY)

            rx0 = max(rx0, DET_XMIN)
            rx1 = min(rx1, DET_XMAX)
            ry0 = max(ry0, DET_YMIN)
            ry1 = min(ry1, DET_YMAX)

            if rx1 <= rx0 or ry1 <= ry0:
                continue

            ix0 = max(0, int(np.floor((rx0 - DET_XMIN) / pixel_mm)))
            ix1 = min(nx - 1, int(np.floor((rx1 - DET_XMIN) / pixel_mm)))
            iy0 = max(0, int(np.floor((ry0 - DET_YMIN) / pixel_mm)))
            iy1 = min(ny - 1, int(np.floor((ry1 - DET_YMIN) / pixel_mm)))

            for jy in range(iy0, iy1 + 1):
                py0 = y_edges[jy]
                py1 = y_edges[jy + 1]
                oy = overlap_1d(py0, py1, ry0, ry1)
                if oy <= 0:
                    continue

                for jx in range(ix0, ix1 + 1):
                    px0 = x_edges[jx]
                    px1 = x_edges[jx + 1]
                    ox = overlap_1d(px0, px1, rx0, rx1)
                    if ox <= 0:
                        continue

                    area = ox * oy
                    sum_map[jy, jx] += val * area
                    weight_map[jy, jx] += area

    # Average where covered by measurement rectangles
    z = np.full_like(sum_map, np.nan, dtype=float)
    covered = weight_map > 0
    z[covered] = sum_map[covered] / weight_map[covered]

    # Extrapolate uncovered detector area
    xx_meas, yy_meas = np.meshgrid(xs, ys)
    z_full = idw_fill_full_detector(
        z=z,
        x_centers=x_centers,
        y_centers=y_centers,
        xs_meas=xx_meas,
        ys_meas=yy_meas,
        m_meas=m,
        power=idw_power
    )

    # Smooth after full fill
    sigma_pix = smooth_sigma_mm / pixel_mm
    z_smooth = gaussian_filter(z_full, sigma=sigma_pix, mode="nearest")

    # Plot
    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    X, Y = np.meshgrid(x_centers, y_centers)

    vmin = np.nanmin(z_smooth)
    vmax = np.nanmax(z_smooth)
    levels = np.linspace(vmin, vmax, nlevels)

    cf = ax.contourf(X, Y, z_smooth, levels=levels, cmap="viridis")
    ax.contour(X, Y, z_smooth, levels=levels, colors="k", linewidths=0.45, alpha=0.4)

    # Detector boundary
    ax.plot(
        [DET_XMIN, DET_XMAX, DET_XMAX, DET_XMIN, DET_XMIN],
        [DET_YMIN, DET_YMIN, DET_YMAX, DET_YMAX, DET_YMIN],
        color="black", linewidth=1.6
    )

    # Original measurement rectangles
    if draw_measurement_boxes:
        for iy, yc in enumerate(ys):
            for ix, xc in enumerate(xs):
                val = m[iy, ix]
                if not np.isfinite(val):
                    continue

                rx0, rx1, ry0, ry1 = rect_edges_from_center(xc, yc, TRIG_WX, TRIG_WY)
                rx0 = max(rx0, DET_XMIN)
                rx1 = min(rx1, DET_XMAX)
                ry0 = max(ry0, DET_YMIN)
                ry1 = min(ry1, DET_YMAX)

                if rx1 > rx0 and ry1 > ry0:
                    ax.add_patch(
                        Rectangle(
                            (rx0, ry0), rx1 - rx0, ry1 - ry0,
                            fill=False, edgecolor="white", linewidth=0.6, alpha=0.35
                        )
                    )

    # Measurement centers
    if draw_centers:
        valid = np.isfinite(m)
        ax.scatter(
            xx_meas[valid], yy_meas[valid],
            s=42, facecolors="none", edgecolors="white", linewidths=1.1, zorder=3
        )

    # Optional numeric labels
    if annotate_values:
        for iy, yc in enumerate(ys):
            for ix, xc in enumerate(xs):
                val = m[iy, ix]
                if np.isfinite(val):
                    ax.text(
                        xc, yc, f"{val:.2f}",
                        ha="center", va="center",
                        fontsize=8, color="white"
                    )

    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title(
        f"Full-detector extrapolated Np.e. contour map (HV={HV_TARGET}, thick={th}, thre={thre_target})\n"
        f"footprint = {TRIG_WX:.0f} × {TRIG_WY:.0f} mm², detector = {DET_SIZE:.0f} × {DET_SIZE:.0f} mm²"
    )

    ax.set_xlim(DET_XMIN - 5, DET_XMAX + 5)
    ax.set_ylim(DET_YMIN - 5, DET_YMAX + 5)
    ax.set_aspect("equal", adjustable="box")

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("Np.e.")

    fig.tight_layout()
    fig.savefig(outname, dpi=200)
    plt.close(fig)

    print(f"Saved: {outname}")


# =========================
# Run
# =========================
for th in [2, 3]:
    m, m_err, thre_target = make_map(th)

    plot_area_averaged_contour_map_full(
        m, th, thre_target,
        outname=f"data_npe_full_extrapolated_HV{HV_TARGET}_thick{th}_thre{thre_target}.png",
        pixel_mm=1.0,
        smooth_sigma_mm=3.0,
        nlevels=24,
        idw_power=2.0,
        draw_measurement_boxes=True,
        draw_centers=False,
        annotate_values=False
    )

print("DONE.")
