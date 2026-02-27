import glob
import re
import numpy as np
import uproot
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# =========================
# 설정
# =========================
PATTERN = "../../../data/KEKAR2023Dec/g4_root/kek_aerogel2_*.root"   # 파일 패턴
TREE = "tree"
BR = "nhMppc"

# 고정 grid (원래 쓰던 것)
xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=int)
ys = np.array([-54, -36, -18, 0, 18, 36, 54], dtype=int)

# 히스토그램/피팅 설정
NBINS = 120
FIT_WINDOW_SIGMA = 2.0   # 피크 주변 몇 sigma 범위로 fit 할지(초기 추정 기반)

# =========================
# 유틸
# =========================
def parse_xy(fname: str):
    """
    kek_aerogel3_xm16_ym54.root  -> x=-16, y=-54
    kek_aerogel3_x16_y54.root    -> x=16,  y=54
    """
    m = re.search(r"_x(?P<x>m?\d+)_y(?P<y>m?\d+)\.root$", fname)
    if not m:
        raise ValueError(f"Cannot parse x/y from: {fname}")

    def tok2int(tok):
        return -int(tok[1:]) if tok.startswith("m") else int(tok)

    x = tok2int(m.group("x"))
    y = tok2int(m.group("y"))
    return x, y

def gauss(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

def fit_gaussian_mean(data):
    """
    data: 1D numpy array (nhMppc)
    return: (mu, mu_err)
    """
    data = np.asarray(data, dtype=float)
    data = data[np.isfinite(data)]
    if len(data) < 20:
        return np.nan, np.nan

    # 히스토그램
    xmin, xmax = np.percentile(data, [0.5, 99.5])
    if not np.isfinite(xmin) or not np.isfinite(xmax) or xmin >= xmax:
        return np.nan, np.nan

    counts, edges = np.histogram(data, bins=NBINS, range=(xmin, xmax))
    centers = 0.5 * (edges[:-1] + edges[1:])

    # 피크 위치(최대 bin)
    ipeak = int(np.argmax(counts))
    mu0 = centers[ipeak]

    # sigma 초기 추정: 피크 주변 데이터의 std로
    # (너무 작아지는 것 방지)
    sigma0 = np.std(data)
    if not np.isfinite(sigma0) or sigma0 <= 0:
        sigma0 = max((xmax - xmin) / 20.0, 1.0)

    # fit 구간: mu0 ± FIT_WINDOW_SIGMA*sigma0
    lo = mu0 - FIT_WINDOW_SIGMA * sigma0
    hi = mu0 + FIT_WINDOW_SIGMA * sigma0
    mask = (centers >= lo) & (centers <= hi) & (counts > 0)

    xfit = centers[mask]
    yfit = counts[mask]
    if len(xfit) < 6:
        return np.nan, np.nan

    # 가중치(포아송): sigma_y ~ sqrt(N)
    yerr = np.sqrt(yfit)
    yerr[yerr == 0] = 1.0

    A0 = float(np.max(yfit))

    try:
        popt, pcov = curve_fit(
            gauss, xfit, yfit,
            p0=[A0, mu0, sigma0],
            sigma=yerr,
            absolute_sigma=True,
            maxfev=10000
        )
        A, mu, sigma = popt
        mu_err = float(np.sqrt(pcov[1, 1])) if pcov is not None and pcov.shape == (3, 3) else np.nan
        return float(mu), mu_err
    except Exception:
        return np.nan, np.nan

# =========================
# 메인: 파일들 읽고 map 채우기
# =========================
files = sorted(glob.glob(PATTERN))
if not files:
    raise RuntimeError(f"No files matched: {PATTERN}")

mean_map = np.full((len(ys), len(xs)), np.nan, dtype=float)
err_map  = np.full((len(ys), len(xs)), np.nan, dtype=float)

for fn in files:
    x, y = parse_xy(fn)
    if x not in xs or y not in ys:
        # grid 밖은 무시
        continue

    with uproot.open(fn) as f:
        t = f[TREE]
        data = t[BR].array(library="np")

    mu, mu_err = fit_gaussian_mean(data)

    iy = int(np.where(ys == y)[0][0])
    ix = int(np.where(xs == x)[0][0])
    mean_map[iy, ix] = mu
    err_map[iy, ix]  = mu_err

    print(f"{fn}: (x={x}, y={y}) mean={mu:.3f} ± {mu_err:.3f}")

# =========================
# Plot heatmap (mean)
# =========================
def plot_map(m, title, outname, fmt="{:.2f}"):
    fig, ax = plt.subplots(figsize=(7, 6))

    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]

    im = ax.imshow(
        m, origin="lower", aspect="auto",
        extent=[xs[0]-dx/2, xs[-1]+dx/2, ys[0]-dy/2, ys[-1]+dy/2]
    )

    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_title(title)
    ax.set_xticks(xs)
    ax.set_yticks(ys)

    for iy, yy in enumerate(ys):
        for ix, xx in enumerate(xs):
            if np.isfinite(m[iy, ix]):
                ax.text(xx, yy, fmt.format(m[iy, ix]),
                        ha="center", va="center", fontsize=10, color="black")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Np.e.")

    fig.tight_layout()
    fig.savefig(outname, dpi=200)
    plt.close(fig)

plot_map(mean_map, "Simulation Np.e. (nhMppc) Gaussian mean", "sim_npe_mean_map_aerogel2.png")
plot_map(err_map,  "Simulation Np.e. (nhMppc) Gaussian mean error", "sim_npe_meanerr_map_aerogel2.png", fmt="{:.3f}")

print("Saved: sim_npe_mean_map.png, sim_npe_meanerr_map.png")
