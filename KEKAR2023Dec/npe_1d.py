import ROOT
import numpy as np
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)

# =========================
# INPUT / BRANCH
# =========================
FILE = "kek_result.root"
TREE = "tree"

BR_X     = "x_pos"
BR_Y     = "y_pos"
BR_THICK = "bac_thick"
BR_HV    = "bac_HV"
BR_THRE  = "bac_thre"
BR_NPE   = "bac_npe"
BR_ERR   = "bac_npe_err"

# =========================
# GRID
# =========================
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=int)
ys_all = np.array([-54, -36, -18, 0, 18, 36, 54], dtype=int)

# use only inside-aerogel rows
ys = np.array([y for y in ys_all if y not in (-54, 54)], dtype=int)

# =========================
# CONDITIONS
# =========================
HV_TARGET    = 58
THICK_TARGET = 2
THRE_TARGET  = 60

xerr = 9.0 / 2.0  # trigger width ±4.5 mm

# =========================
# Load ROOT
# =========================
f = ROOT.TFile.Open(FILE)
t = f.Get(TREE)
if not t:
    raise RuntimeError("Tree not found")

x, y, thick, hv, thre, npe, npe_err = [], [], [], [], [], [], []

for ev in t:
    x.append(int(getattr(ev, BR_X)))
    y.append(int(getattr(ev, BR_Y)))
    thick.append(int(getattr(ev, BR_THICK)))
    hv.append(int(getattr(ev, BR_HV)))
    thre.append(int(getattr(ev, BR_THRE)))
    npe.append(float(getattr(ev, BR_NPE)))
    npe_err.append(float(getattr(ev, BR_ERR)))

x = np.array(x)
y = np.array(y)
thick = np.array(thick)
hv = np.array(hv)
thre = np.array(thre)
npe = np.array(npe)
npe_err = np.array(npe_err)

# global selection
sel0 = (hv == HV_TARGET) & (thick == THICK_TARGET) & (thre == THRE_TARGET)

# =========================
# helper
# =========================
def mean_and_err(vals, errs):
    vals = np.asarray(vals)
    errs = np.asarray(errs)
    m = np.mean(vals)
    em = np.sqrt(np.sum(errs**2)) / len(vals)
    return m, em

# =========================
# Plot overlay
# =========================
fig, ax = plt.subplots(figsize=(7.5, 5.5))

markers = ["o", "s", "^", "D", "v"]
colors  = ["C0", "C1", "C2", "C3", "C4"]

for iy, yy in enumerate(ys):
    xpts, ypts, yerrs = [], [], []

    for xx in xs:
        mask = sel0 & (x == xx) & (y == yy)
        if np.any(mask):
            m, em = mean_and_err(npe[mask], npe_err[mask])
            xpts.append(xx)
            ypts.append(m)
            yerrs.append(em)

    if len(xpts) == 0:
        continue

    ax.errorbar(
        xpts, ypts,
        xerr=xerr,
        yerr=yerrs,
        fmt=markers[iy],
        color=colors[iy],
        capsize=3,
        label=f"y = {yy} mm"
    )

ax.set_xlabel("X position (mm)")
#ax.set_ylim(10, 40)
ax.set_ylabel("Np.e.")
ax.set_xticks(xs)
ax.grid(True, alpha=0.3)
#ax.legend(title="Beam Y position")

fig.tight_layout()
fig.savefig("npe_vs_x_overlay_by_y_2.png", dpi=200)
plt.close(fig)

