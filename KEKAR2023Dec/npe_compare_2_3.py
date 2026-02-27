import ROOT
import numpy as np
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)

# =========================
# INPUT
# =========================
FILE = "kek_result.root"
TREE = "tree"

BR_X     = "x_pos"
BR_Y     = "y_pos"
BR_THICK = "bac_thick"
BR_HV    = "bac_HV"
BR_THRE  = "bac_thre"
BR_NPE   = "bac_npe"

# =========================
# FIXED GRID (mm)
# =========================
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=int)
ys = np.array([-54, -36, -18, 0, 18, 36, 54], dtype=int)

# =========================
# Selection rules (same as yours)
# =========================
HV_TARGET = 58
THRE_BY_THICK = {2: 60, 3: 75}

# If you want to exclude y=±54 (optional)
EXCLUDE_EDGE_Y = False  # True면 y=±54 제외
EDGE_Y = set([-54, 54])

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
# Helper: mean and error per (x,y) for given thickness
# =========================
def mean_err_map_for_thickness(th):
    if th not in THRE_BY_THICK:
        raise ValueError(f"No threshold rule for thickness={th}")
    thre_target = THRE_BY_THICK[th]

    sel = (hv == HV_TARGET) & (thick == th) & (thre == thre_target)

    # maps
    mean_map = np.full((len(ys), len(xs)), np.nan, dtype=float)
    err_map  = np.full((len(ys), len(xs)), np.nan, dtype=float)
    n_map    = np.zeros((len(ys), len(xs)), dtype=int)

    for iy, yy in enumerate(ys):
        if EXCLUDE_EDGE_Y and (yy in EDGE_Y):
            continue
        for ix, xx in enumerate(xs):
            mask = sel & (x == xx) & (y == yy)
            if np.any(mask):
                vals = npe[mask]
                n = int(np.sum(mask))
                mu = float(np.mean(vals))
                # error of the mean (SEM). If you prefer RMS as error: np.std(vals, ddof=1)
                sig = float(np.std(vals, ddof=1)) if n > 1 else 0.0
                sem = sig / np.sqrt(n) if n > 1 else np.nan

                mean_map[iy, ix] = mu
                err_map[iy, ix]  = sem
                n_map[iy, ix]    = n

    return mean_map, err_map, n_map, thre_target

# =========================
# Build maps for 2-tile and 3-tile
# =========================
m2, e2, n2, thr2 = mean_err_map_for_thickness(2)
m3, e3, n3, thr3 = mean_err_map_for_thickness(3)

print(f"[thick=2] thre={thr2}, filled points: {np.sum(np.isfinite(m2))}")
print(f"[thick=3] thre={thr3}, filled points: {np.sum(np.isfinite(m3))}")

# =========================
# Ratio per position: r = m3/m2 (only where both exist)
# =========================
ratios = []
ratio_errs = []
pos_list = []

for iy, yy in enumerate(ys):
    if EXCLUDE_EDGE_Y and (yy in EDGE_Y):
        continue
    for ix, xx in enumerate(xs):
        if np.isfinite(m2[iy, ix]) and np.isfinite(m3[iy, ix]) and m2[iy, ix] != 0:
            r = float(m3[iy, ix] / m2[iy, ix])

            # error propagation (optional): r * sqrt( (e3/m3)^2 + (e2/m2)^2 )
            # if SEM is nan (n==1), this becomes nan; we keep it.
            if np.isfinite(e2[iy, ix]) and np.isfinite(e3[iy, ix]) and m3[iy, ix] != 0:
                re = r * np.sqrt((e3[iy, ix] / m3[iy, ix])**2 + (e2[iy, ix] / m2[iy, ix])**2)
            else:
                re = np.nan

            ratios.append(r)
            ratio_errs.append(re)
            pos_list.append((xx, yy))

ratios = np.array(ratios, dtype=float)
ratio_errs = np.array(ratio_errs, dtype=float)

print(f"Ratio entries: {len(ratios)}")

# =========================
# Plot 1) Histogram of ratios
# =========================
plt.figure(figsize=(7, 5))
# bins: you can tune. Here: auto with a cap.
nbins = 20
finite = np.isfinite(ratios)
plt.hist(ratios[finite], bins=nbins)

plt.xlabel(r"$\langle N_{\mathrm{pe}}\rangle_{\mathrm{3\,tiles}} / \langle N_{\mathrm{pe}}\rangle_{\mathrm{2\,tiles}}$")
plt.ylabel("Counts")
plt.title(f"Ratio distribution (HV={HV_TARGET}, thre2={thr2}, thre3={thr3})")
plt.tight_layout()
plt.savefig("ratio_hist_aerogel3_over_aerogel2.png", dpi=200)
plt.close()
print("Saved: ratio_hist_aerogel3_over_aerogel2.png")

# =========================
# Plot 2) (Optional but useful) scatter: ratio vs position index
# =========================
# If you want a quick check of outliers by position
plt.figure(figsize=(8, 4))
idx = np.arange(len(ratios))
plt.errorbar(idx, ratios, yerr=np.where(np.isfinite(ratio_errs), ratio_errs, 0.0), fmt='o', capsize=2)
plt.xlabel("Position index (each valid grid point)")
plt.ylabel(r"Ratio $\langle N_{\mathrm{pe}}\rangle_{3} / \langle N_{\mathrm{pe}}\rangle_{2}$")
plt.title("Point-by-point ratio with propagated errors")
plt.tight_layout()
plt.savefig("ratio_point_by_point.png", dpi=200)
plt.close()
print("Saved: ratio_point_by_point.png")

# =========================
# Print top outliers (optional)
# =========================
if len(ratios) > 0:
    order = np.argsort(ratios)
    print("\nLowest 5 ratios:")
    for k in order[:5]:
        xx, yy = pos_list[k]
        print(f"  (x={xx}, y={yy}) r={ratios[k]:.3f}  m2={m2[np.where(ys==yy)[0][0], np.where(xs==xx)[0][0]]:.3f}  m3={m3[np.where(ys==yy)[0][0], np.where(xs==xx)[0][0]]:.3f}")

    print("\nHighest 5 ratios:")
    for k in order[-5:][::-1]:
        xx, yy = pos_list[k]
        print(f"  (x={xx}, y={yy}) r={ratios[k]:.3f}  m2={m2[np.where(ys==yy)[0][0], np.where(xs==xx)[0][0]]:.3f}  m3={m3[np.where(ys==yy)[0][0], np.where(xs==xx)[0][0]]:.3f}")
