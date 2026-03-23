import ROOT
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse
from scipy.optimize import curve_fit

ROOT.gROOT.SetBatch(True)

# =========================================================
# INPUT ROOT (NPE)
# =========================================================
FILE = "kek_result.root"
TREE = "tree"

BR_X     = "x_pos"
BR_Y     = "y_pos"
BR_THICK = "bac_thick"
BR_HV    = "bac_HV"
BR_THRE  = "bac_thre"
BR_NPE   = "bac_npe"

# =========================================================
# INPUT ROOT (BEAM)
# =========================================================
BEAM_FILE   = "../../../data/E72_Simul/E72_Beam_Simul_735.root"
BEAM_TREE   = None      # None -> first TTree automatically
BEAM_BRANCH = "BAC"
BEAM_MAX_ENTRIES = None

# =========================================================
# GRID / SELECTION
# =========================================================
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=float)
ys = np.array([-54, -36, -18, 0, 18, 36, 54], dtype=float)

HV_TARGET = 58
THRE_BY_THICK = {2: 60, 3: 75}

# measurement rectangle size
TRIG_WX = 9.0
TRIG_WY = 13.0

# detector size
DET_SIZE = 115.0
DET_HALF = DET_SIZE / 2.0
DET_XMIN, DET_XMAX = -DET_HALF, +DET_HALF
DET_YMIN, DET_YMAX = -DET_HALF, +DET_HALF

# beam ellipse
BEAM_SIGMA_LEVELS = [1, 2]

# smoothness of reconstructed NPE map
# start with these; increase slightly if still rough
RECO_SIGMA_X = 7.0   # mm
RECO_SIGMA_Y = 9.0   # mm

# plotting density
PIXEL_MM = 1.0
NLEVELS = 14

# beam histogram setting
BEAM_BINS = 70


# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------
def rect_edges_from_center(xc, yc, wx, wy):
    return xc - wx / 2.0, xc + wx / 2.0, yc - wy / 2.0, yc + wy / 2.0


def find_first_ttree(root_file):
    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("TTree"):
            return obj.GetName()
    return None


# ---------------------------------------------------------
# Read NPE
# ---------------------------------------------------------
f = ROOT.TFile.Open(FILE)
t = f.Get(TREE)
if not t:
    raise RuntimeError(f"Cannot find TTree '{TREE}' in {FILE}")

x, y, thick, hv, thre, npe = [], [], [], [], [], []

for ev in t:
    x.append(float(getattr(ev, BR_X)))
    y.append(float(getattr(ev, BR_Y)))
    thick.append(int(getattr(ev, BR_THICK)))
    hv.append(int(getattr(ev, BR_HV)))
    thre.append(int(getattr(ev, BR_THRE)))
    npe.append(float(getattr(ev, BR_NPE)))

x     = np.array(x, dtype=float)
y     = np.array(y, dtype=float)
thick = np.array(thick, dtype=int)
hv    = np.array(hv, dtype=int)
thre  = np.array(thre, dtype=int)
npe   = np.array(npe, dtype=float)


def make_map(th):
    thre_target = THRE_BY_THICK[th]
    sel = (hv == HV_TARGET) & (thick == th) & (thre == thre_target)

    print(f"[NPE] thick={th}, thre={thre_target}, selected entries = {int(np.sum(sel))}")

    m = np.full((len(ys), len(xs)), np.nan, dtype=float)

    for iy, yy in enumerate(ys):
        for ix, xx in enumerate(xs):
            mask = sel & (x == xx) & (y == yy)
            if np.any(mask):
                m[iy, ix] = np.mean(npe[mask])

    return m, thre_target


# ---------------------------------------------------------
# Read beam positions
# BAC is vector<TParticle>, so use BAC[0]
# ---------------------------------------------------------
def read_beam_positions(filename, tree_name=None, branch_name="BAC", max_entries=None):
    rf = ROOT.TFile.Open(filename)
    if not rf or rf.IsZombie():
        raise RuntimeError(f"Cannot open beam file: {filename}")

    if tree_name is None:
        tree_name = find_first_ttree(rf)
        if tree_name is None:
            raise RuntimeError(f"No TTree found in {filename}")

    tt = rf.Get(tree_name)
    if not tt:
        raise RuntimeError(f"Cannot find TTree '{tree_name}' in {filename}")

    xv, yv = [], []
    nent = tt.GetEntries()
    if max_entries is not None:
        nent = min(nent, max_entries)

    print(f"[BEAM] reading {nent} entries from {filename}:{tree_name}")

    for i, ev in enumerate(tt):
        if max_entries is not None and i >= max_entries:
            break

        bac_vec = getattr(ev, branch_name)

        if len(bac_vec) == 0:
            continue

        # use the first BAC particle
        bac = bac_vec[0]

        xv.append(float(bac.Vx()))
        yv.append(float(bac.Vy()))

    xv = np.array(xv, dtype=float)
    yv = np.array(yv, dtype=float)

    good = np.isfinite(xv) & np.isfinite(yv)
    xv = xv[good]
    yv = yv[good]

    print(f"[BEAM] valid positions = {len(xv)}")
    return xv, yv


# ---------------------------------------------------------
# 2D Gaussian fit for beam
# ---------------------------------------------------------
def gaussian2d_rotated(coords, A, x0, y0, sx, sy, theta, C):
    xg, yg = coords
    ct = np.cos(theta)
    st = np.sin(theta)

    xp =  ct * (xg - x0) + st * (yg - y0)
    yp = -st * (xg - x0) + ct * (yg - y0)

    expo = -0.5 * ((xp / sx) ** 2 + (yp / sy) ** 2)
    return A * np.exp(expo) + C


def fit_beam_2d_gaussian(xbeam, ybeam, bins=70, xrange=(-57.5, 57.5), yrange=(-57.5, 57.5)):
    H, xedges, yedges = np.histogram2d(
        xbeam, ybeam,
        bins=[bins, bins],
        range=[xrange, yrange]
    )
    H = H.T

    xcent = 0.5 * (xedges[:-1] + xedges[1:])
    ycent = 0.5 * (yedges[:-1] + yedges[1:])
    X, Y = np.meshgrid(xcent, ycent)
    Z = H.astype(float)

    total = np.sum(Z)
    if total <= 0:
        raise RuntimeError("Beam histogram is empty")

    x0_init = np.sum(X * Z) / total
    y0_init = np.sum(Y * Z) / total

    dx = X - x0_init
    dy = Y - y0_init
    sxx = np.sum(Z * dx * dx) / total
    syy = np.sum(Z * dy * dy) / total
    sxy = np.sum(Z * dx * dy) / total

    cov = np.array([[sxx, sxy], [sxy, syy]], dtype=float)
    evals, evecs = np.linalg.eigh(cov)
    order = np.argsort(evals)[::-1]
    evals = evals[order]
    evecs = evecs[:, order]

    sx_init = np.sqrt(max(evals[0], 1.0))
    sy_init = np.sqrt(max(evals[1], 1.0))
    theta_init = np.arctan2(evecs[1, 0], evecs[0, 0])

    A_init = np.max(Z) - np.min(Z)
    C_init = np.min(Z)

    p0 = [A_init, x0_init, y0_init, sx_init, sy_init, theta_init, C_init]
    lower = [0.0, xrange[0], yrange[0], 1.0, 1.0, -np.pi/2, 0.0]
    upper = [np.inf, xrange[1], yrange[1], 100.0, 100.0, np.pi/2, np.max(Z)]

    xdata = np.vstack((X.ravel(), Y.ravel()))
    ydata = Z.ravel()

    try:
        popt, _ = curve_fit(
            gaussian2d_rotated,
            xdata, ydata,
            p0=p0,
            bounds=(lower, upper),
            maxfev=50000
        )
    except Exception as e:
        print("[BEAM] curve_fit failed, using moment estimate")
        print("[BEAM] reason:", e)
        popt = p0

    A, x0, y0, sx, sy, theta, C = popt

    if sy > sx:
        sx, sy = sy, sx
        theta += np.pi / 2.0

    while theta > np.pi / 2:
        theta -= np.pi
    while theta < -np.pi / 2:
        theta += np.pi

    print(f"[BEAM] center = ({x0:.2f}, {y0:.2f}) mm")
    print(f"[BEAM] sigma major/minor = ({sx:.2f}, {sy:.2f}) mm")
    print(f"[BEAM] angle = {np.degrees(theta):.2f} deg")

    return {
        "x0": x0,
        "y0": y0,
        "sx": sx,
        "sy": sy,
        "theta": theta
    }


# ---------------------------------------------------------
# Smooth NPE reconstruction using Gaussian kernels
# ---------------------------------------------------------
def build_gaussian_map(m, xs, ys, x_centers, y_centers, sigma_x=7.0, sigma_y=9.0):
    X, Y = np.meshgrid(x_centers, y_centers)

    z_sum = np.zeros_like(X, dtype=float)
    w_sum = np.zeros_like(X, dtype=float)

    for iy, yc in enumerate(ys):
        for ix, xc in enumerate(xs):
            val = m[iy, ix]
            if not np.isfinite(val):
                continue

            dx = X - xc
            dy = Y - yc

            kernel = np.exp(
                -0.5 * ((dx / sigma_x) ** 2 + (dy / sigma_y) ** 2)
            )

            z_sum += val * kernel
            w_sum += kernel

    z = np.full_like(z_sum, np.nan, dtype=float)
    good = w_sum > 1e-12
    z[good] = z_sum[good] / w_sum[good]
    return z


# ---------------------------------------------------------
# Beam ellipses: 1σ and 2σ only, no center mark
# ---------------------------------------------------------
def add_beam_sigma_ellipses(ax, fit, sigma_levels=(1, 2),
                            edgecolor="red", linewidth=1.8,
                            label_fontsize=15):
    x0 = fit["x0"]
    y0 = fit["y0"]
    sx = fit["sx"]
    sy = fit["sy"]
    theta = fit["theta"]
    theta_deg = np.degrees(theta)

    ct = np.cos(theta)
    st = np.sin(theta)

    for nsig in sigma_levels:
        ell = Ellipse(
            xy=(0, y0),
            width=2.0 * nsig * sx,
            height=2.0 * nsig * sy,
            angle=theta_deg,
            fill=False,
            edgecolor=edgecolor,
            linewidth=linewidth
        )
        ax.add_patch(ell)

        # small sigma label near major axis
        lx = x0 + nsig * sx * ct
        ly = y0 + nsig * sx * st
        ax.text(
            0, ly + 3.0,
            f"{nsig}$\\sigma$",
            color=edgecolor,
            fontsize=label_fontsize,
            ha="left",
            va="bottom"
        )


# ---------------------------------------------------------
# Plot NPE + beam sigma overlay
# ---------------------------------------------------------
def plot_area_averaged_contour_map_smooth(
    m, th, thre_target, outname,
    beam_fit=None,
    pixel_mm=1.0,
    sigma_x=7.0,
    sigma_y=9.0,
    nlevels=14,
):
    m = np.asarray(m, dtype=float)

    x_edges = np.arange(DET_XMIN, DET_XMAX + pixel_mm, pixel_mm)
    y_edges = np.arange(DET_YMIN, DET_YMAX + pixel_mm, pixel_mm)

    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    X, Y = np.meshgrid(x_centers, y_centers)

    # smooth reconstructed map
    z_smooth = build_gaussian_map(
        m, xs, ys,
        x_centers, y_centers,
        sigma_x=sigma_x,
        sigma_y=sigma_y
    )

    fig, ax = plt.subplots(figsize=(7.8, 6.8))

    #vmin = np.nanmin(z_smooth)
    #vmax = np.nanmax(z_smooth)
    vmin = 0
    vmax = 50
    levels = np.linspace(vmin, vmax, nlevels)

    cf = ax.contourf(X, Y, z_smooth, levels=levels, cmap="viridis")

    # make contour lines weak; comment out next line if you want even cleaner look
    ax.contour(X, Y, z_smooth, levels=levels, colors="k", linewidths=0.25, alpha=0.22)

    # detector boundary
    ax.plot(
        [DET_XMIN, DET_XMAX, DET_XMAX, DET_XMIN, DET_XMIN],
        [DET_YMIN, DET_YMIN, DET_YMAX, DET_YMAX, DET_YMIN],
        color="black", linewidth=1.6
    )

    # measurement rectangles only
    for iy, yc in enumerate(ys):
        for ix, xc in enumerate(xs):
            val = m[iy, ix]
            if not np.isfinite(val):
                continue

            rx0, rx1, ry0, ry1 = rect_edges_from_center(xc, yc, TRIG_WX, TRIG_WY)
            ax.add_patch(
                Rectangle(
                    (rx0, ry0), rx1 - rx0, ry1 - ry0,
                    fill=False, edgecolor="white", linewidth=0.6, alpha=1
                )
            )

    # beam sigma ellipses
    if beam_fit is not None:
        add_beam_sigma_ellipses(
            ax, beam_fit,
            sigma_levels=BEAM_SIGMA_LEVELS,
            edgecolor="red",
            linewidth=1.8,
            label_fontsize=12
        )

    ax.set_xlabel("X (mm)",fontsize=16)
    ax.set_ylabel("Y (mm)",fontsize=16)
    ax.set_xticks([-40, -20, 0, 20, 40])
    ax.set_yticks([-40, -20, 0, 20, 40])
    ax.set_xlim(-57.5, 57.5)
    ax.set_ylim(-57.5, 57.5)
    ax.tick_params(axis='both', labelsize=14)
    
    ax.set_aspect("equal", adjustable="box")

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_ticks(np.arange(0, 51, 10))
    cbar.set_label("Np.e.",fontsize=18)
    cbar.ax.tick_params(labelsize=14)

    fig.tight_layout()
    fig.savefig(outname, dpi=220)
    plt.close(fig)
    print(f"Saved: {outname}")


# =========================================================
# MAIN
# =========================================================
if __name__ == "__main__":
    xbeam, ybeam = read_beam_positions(
        BEAM_FILE,
        tree_name=BEAM_TREE,
        branch_name=BEAM_BRANCH,
        max_entries=BEAM_MAX_ENTRIES
    )

    beam_fit = fit_beam_2d_gaussian(
        xbeam, ybeam,
        bins=BEAM_BINS,
        xrange=(DET_XMIN, DET_XMAX),
        yrange=(DET_YMIN, DET_YMAX)
    )

    for th in [2, 3]:
        m, thre_target = make_map(th)

        plot_area_averaged_contour_map_smooth(
            m, th, thre_target,
            outname=f"data_npe_smooth_HV{HV_TARGET}_thick{th}_thre{thre_target}.png",
            beam_fit=beam_fit,
            pixel_mm=PIXEL_MM,
            sigma_x=RECO_SIGMA_X,
            sigma_y=RECO_SIGMA_Y,
            nlevels=NLEVELS,
        )

    print("DONE.")
