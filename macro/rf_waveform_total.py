import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def load_and_interp(csv_path, dt=0.0, dv=0.0, npts=50000):
    df = pd.read_csv(csv_path)
    if df.shape[1] < 2:
        df = pd.read_csv(csv_path, header=None)

    t = df.iloc[:, 0].astype(float).to_numpy()
    v = df.iloc[:, 1].astype(float).to_numpy()

    # 오프셋 적용
    t = t + dt
    v = v + dv

    # 정렬
    idx = np.argsort(t)
    t, v = t[idx], v[idx]

    # 보간
    t_dense = np.linspace(t.min(), t.max(), npts)
    v_dense = np.interp(t_dense, t, v)

    return t_dense, v_dense, t, v

import numpy as np

def negative_pulse_area(t, v, tmin=None, tmax=None):
    t = np.asarray(t)
    v = np.asarray(v)

    # 구간 제한 (원하면)
    if tmin is not None or tmax is not None:
        if tmin is None: tmin = t.min()
        if tmax is None: tmax = t.max()
        m = (t >= tmin) & (t <= tmax)
        t, v = t[m], v[m]

    # 0 아래만 남기기
    v_neg = np.minimum(v, 0.0)

    # 음수 면적을 양수로 반환 (단위: V·ns)
    area = -np.trapz(v_neg, t)
    return area




dt_rf = 10.0     # ns
dv_rf = -70.0     # V

dt_rk = 10.0    # ns (예: rk를 왼쪽으로 이동)
dv_rk = -365.0    # V  (예: rk baseline 보정)
# ======================

t_rf, v_rf, t_rf_raw, v_rf_raw = load_and_interp(
    "rf_1k.csv", dt_rf, dv_rf
)

t_rk, v_rk, t_rk_raw, v_rk_raw = load_and_interp(
    "rf_10k.csv", dt_rk, dv_rk
)

v_rf_smooth = savgol_filter(v_rf, window_length=1001, polyorder=5)
v_rk_smooth = savgol_filter(v_rk, window_length=1001, polyorder=5)

plt.figure(figsize=(9,4))
plt.plot(t_rf, v_rf_smooth, linewidth=2, label=r"$\mathrm{R_{F}} = 1\,\mathrm{k}\Omega$")
plt.plot(t_rk, v_rk_smooth, linewidth=2, label=r"$\mathrm{R_{F}} = 10\,\mathrm{k}\Omega$")

area_rf = negative_pulse_area(t_rf, v_rf_smooth, tmin=0, tmax=120)
area_rk = negative_pulse_area(t_rk, v_rk_smooth, tmin=0, tmax=120)

print("rf area (V*ns) =", area_rf)
print("rk area (V*ns) =", area_rk)

# 원본 포인트 확인용
#plt.plot(t_rf_raw, v_rf_raw, 'o', alpha=0.25)
#plt.plot(t_rk_raw, v_rk_raw, 'o', alpha=0.25)

plt.xlabel("Time (ns)")
plt.ylabel("Voltage (mV)")
plt.xlim(0, 130) 
#plt.title("Waveform comparison with manual offsets")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("rf_waveform.pdf",bbox_inches="tight",transparent=True)
plt.show()
