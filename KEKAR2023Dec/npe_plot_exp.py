import ROOT
import numpy as np
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)

FILE = "kek_result.root"
TREE = "tree"

# branches
BR_X     = "x_pos"
BR_Y     = "y_pos"
BR_THICK = "bac_thick"   # or bac_thickness
BR_HV    = "bac_HV"
BR_THRE  = "bac_thre"
BR_NPE   = "bac_npe"

# =========================
# FIXED X/Y GRID (mm)
# =========================
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=int)
ys = np.array([-54, -36, -18, 0, 18, 36, 54], dtype=int)

# =========================
# Selection rules
# =========================
HV_TARGET = 58
THRE_BY_THICK = {2: 60, 3: 75}

# =========================
# Load ROOT file / tree
# =========================
f = ROOT.TFile.Open(FILE)
t = f.Get(TREE)
if not t:
    raise RuntimeError(f"Cannot find TTree '{TREE}' in {FILE}")

# =========================
# Read data into arrays
# =========================
x, y, thick, hv, thre, npe = [], [], [], [], [], []

for ev in t:
    x.append(int(getattr(ev, BR_X)))
    y.append(int(getattr(ev, BR_Y)))
    thick.append(int(getattr(ev, BR_THICK)))
    hv.append(int(getattr(ev, BR_HV)))
    thre.append(int(getattr(ev, BR_THRE)))
    npe.append(float(getattr(ev, BR_NPE)))

x     = np.array(x, dtype=int)
y     = np.array(y, dtype=int)
thick = np.array(thick, dtype=int)
hv    = np.array(hv, dtype=int)
thre  = np.array(thre, dtype=int)
npe   = np.array(npe, dtype=float)

# =========================
# Build Np.e. map for a given thickness with conditions
# =========================
def make_map(th):
    if th not in THRE_BY_THICK:
        raise ValueError(f"No threshold rule for thickness={th}")

    thre_target = THRE_BY_THICK[th]

    # global selection mask for this thickness
    sel = (hv == HV_TARGET) & (thick == th) & (thre == thre_target)

    nsel = int(np.sum(sel))
    print(f"[thick={th}] select: HV=={HV_TARGET} & thre=={thre_target} -> {nsel} entries")

    m = np.full((len(ys), len(xs)), np.nan, dtype=float)

    for iy, yy in enumerate(ys):
        for ix, xx in enumerate(xs):
            mask = sel & (x == xx) & (y == yy)
            if np.any(mask):
                m[iy, ix] = np.mean(npe[mask])  # 평균 Np.e.

    return m, thre_target

# =========================
# Plot function
# =========================
def plot_map(m, th, thre_target, outname):
    fig, ax = plt.subplots(figsize=(7, 6))

    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]

    im = ax.imshow(
        m,
        origin="lower",
        aspect="auto",
        extent=[xs[0]-dx/2, xs[-1]+dx/2, ys[0]-dy/2, ys[-1]+dy/2]
    )

    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title(f"Np.e. map (HV={HV_TARGET}, thick={th}, thre={thre_target})")

    ax.set_xticks(xs)
    ax.set_yticks(ys)

    # annotate
    for iy, yy in enumerate(ys):
        for ix, xx in enumerate(xs):
            if np.isfinite(m[iy, ix]):
                ax.text(xx, yy, f"{m[iy, ix]:.2f}",
                        ha="center", va="center", fontsize=10, color="black")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Np.e.")

    fig.tight_layout()
    fig.savefig(outname, dpi=200)
    plt.close(fig)

# =========================
# Make plots for thickness 2 and 3
# =========================
for th in [2, 3]:
    m, thre_target = make_map(th)
    plot_map(m, th, thre_target, f"npe_map_HV{HV_TARGET}_thick{th}_thre{thre_target}.png")

print("DONE.")
