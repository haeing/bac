import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
import matplotlib.cm as cm

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
# optional (if exists)
BR_NPE_ERR = "bac_npe_err"  # ok if not used

# =========================
# GRID (mm)
# =========================
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=float)
#ys = np.array([-54, -36, -18, 0, 18, 36, 54], dtype=float)
ys = np.array([-36, -18, 0, 18, 36], dtype=float)

# =========================
# Selection rules
# =========================
HV_TARGET = 58
THRE_BY_THICK = {2: 60, 3: 75}

# =========================
# Trigger footprint (mm)
# =========================
TRIG_WX = 9.0   # x-direction
TRIG_WY = 13.0  # y-direction

# =========================
# Aerogel active area (mm)
# =========================
AERO_SIZE = 115.0
AERO_HALF = AERO_SIZE / 2.0  # 57.5 mm
# assume aerogel center at (0,0) in this map coordinate
AERO_XMIN, AERO_XMAX = -AERO_HALF, +AERO_HALF
AERO_YMIN, AERO_YMAX = -AERO_HALF, +AERO_HALF


# -------------------------
# Utils: rectangle clipping to aerogel box
# -------------------------
def intersect_rect_with_box(xc, yc, wx, wy, xmin, xmax, ymin, ymax):
    """
    Returns (x0, y0, w, h) of the intersection rectangle in global coords,
    or None if no overlap.
    Rectangle is centered at (xc,yc) with width wx, height wy.
    Box is [xmin,xmax]×[ymin,ymax].
    """
    rx0 = xc - wx/2
    rx1 = xc + wx/2
    ry0 = yc - wy/2
    ry1 = yc + wy/2

    ix0 = max(rx0, xmin)
    ix1 = min(rx1, xmax)
    iy0 = max(ry0, ymin)
    iy1 = min(ry1, ymax)

    if ix1 <= ix0 or iy1 <= iy0:
        return None
    return (ix0, iy0, ix1-ix0, iy1-iy0)


# -------------------------
# Read ROOT tree → numpy arrays
# -------------------------
f = ROOT.TFile.Open(FILE)
t = f.Get(TREE)
if not t:
    raise RuntimeError(f"Cannot find TTree '{TREE}' in {FILE}")

x, y, thick, hv, thre, npe = [], [], [], [], [], []
npe_err = []

# branch existence check (optional)
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
# Build map (mean Npe) at each (x,y)
# -------------------------
def make_map(th):
    if th not in THRE_BY_THICK:
        raise ValueError(f"No threshold rule for thickness={th}")
    thre_target = THRE_BY_THICK[th]

    sel = (hv == HV_TARGET) & (thick == th) & (thre == thre_target)
    print(f"[thick={th}] select: HV=={HV_TARGET} & thre=={thre_target} -> {int(np.sum(sel))} entries")

    m = np.full((len(ys), len(xs)), np.nan, dtype=float)
    m_err = np.full((len(ys), len(xs)), np.nan, dtype=float)

    for iy, yy in enumerate(ys):
        for ix, xx in enumerate(xs):
            mask = sel & (x == xx) & (y == yy)
            if np.any(mask):
                m[iy, ix] = np.mean(npe[mask])
                if has_npe_err:
                    # 평균 에러를 어떻게 정의할지: 보통은 sqrt(sum(err^2))/N (uncorrelated 가정)
                    # 혹은 sample std/sqrt(N)도 가능. 여기선 branch err를 "measurement err"로 보고 합성.
                    ee = npe_err[mask]
                    m_err[iy, ix] = np.sqrt(np.sum(ee*ee)) / np.sum(mask)

    return m, m_err, thre_target


# -------------------------
# Plot: trigger footprint + aerogel boundary (clipped)
# -------------------------
def plot_trigger_footprint_map(m, th, thre_target, outname,
                               annotate=True, show_outside_centers=True):
    """
    - Each grid point is drawn as a 9×13 mm² rectangle.
    - Only the overlap with aerogel 115×115 mm² is filled (clipped).
    - Aerogel boundary is drawn as a black square.
    """
    m = np.asarray(m, dtype=float)

    finite = np.isfinite(m)
    if not np.any(finite):
        raise RuntimeError("Map has no finite values.")

    vmin = np.nanmin(m)
    vmax = np.nanmax(m)
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap("viridis")

    fig, ax = plt.subplots(figsize=(7.2, 6.2))

    # Aerogel boundary
    aero_rect = Rectangle((AERO_XMIN, AERO_YMIN), AERO_SIZE, AERO_SIZE,
                          facecolor="none", edgecolor="black", linewidth=1.5)
    ax.add_patch(aero_rect)

    # Draw rectangles
    for iy, y0 in enumerate(ys):
        for ix, x0 in enumerate(xs):
            val = m[iy, ix]
            if not np.isfinite(val):
                continue

            clipped = intersect_rect_with_box(
                x0, y0, TRIG_WX, TRIG_WY,
                AERO_XMIN, AERO_XMAX, AERO_YMIN, AERO_YMAX
            )

            # if no overlap with aerogel: optionally still show center marker (not filled)
            if clipped is None:
                if show_outside_centers:
                    ax.plot([x0], [y0], marker="x", markersize=6)
                    if annotate:
                        ax.text(x0, y0, f"{val:.1f}", ha="center", va="center", fontsize=8)
                continue

            x_left, y_bot, w, h = clipped
            rect = Rectangle((x_left, y_bot), w, h,
                             facecolor=cmap(norm(val)),
                             edgecolor="k", linewidth=0.5)
            ax.add_patch(rect)

            if annotate:
                ax.text(x0, y0, f"{val:.2f}", ha="center", va="center",
                        fontsize=9, color="black")

    ax.set_xlabel("X (mm)", fontsize=20)
    ax.set_ylabel("Y (mm)", fontsize=20)
    ax.set_title(f"Data Np.e. footprint map (HV={HV_TARGET}, thick={th}, thre={thre_target})\n"
                 f"Trigger area = {TRIG_WX:.0f}×{TRIG_WY:.0f} mm$^2$, Aerogel = {AERO_SIZE:.0f}×{AERO_SIZE:.0f} mm$^2$")

    # limits: show a bit margin around aerogel
    margin = 10
    ax.set_xlim(AERO_XMIN - margin, AERO_XMAX + margin)
    ax.set_ylim(AERO_YMIN - margin, AERO_YMAX + margin)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xticks(xs)
    ax.set_yticks(ys)
    ax.grid(False)

    # colorbar
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Np.e.", fontsize=20)

    fig.tight_layout()
    fig.savefig(outname, dpi=200)
    plt.close(fig)
    print(f"Saved: {outname}")


# =========================
# Run for thick 2 and 3
# =========================
for th in [2, 3]:
    m, m_err, thre_target = make_map(th)
    plot_trigger_footprint_map(
        m, th, thre_target,
        outname=f"data_npe_footprint_HV{HV_TARGET}_thick{th}_thre{thre_target}.png",
        annotate=True,
        show_outside_centers=True
    )

print("DONE.")
