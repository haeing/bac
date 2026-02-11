import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

df = pd.read_csv("rf_1k.csv")

t = df.iloc[:, 0].to_numpy()
v = df.iloc[:, 1].to_numpy()

# 시간 정렬 (중요)
idx = np.argsort(t)
t, v = t[idx], v[idx]

# spline 보간 (포인트 수 늘리기)
t_smooth = np.linspace(t.min(), t.max(), 1000)
spline = make_interp_spline(t, v, k=3)  # cubic spline
v_smooth = spline(t_smooth)

plt.figure(figsize=(8,4))
plt.plot(t, v, 'o', label="Raw data", alpha=0.6)
plt.plot(t_smooth, v_smooth, '-', label="Cubic spline", linewidth=2)
plt.xlabel("Time (ns)")
plt.ylabel("Voltage (V)")
plt.title("Waveform (smoothed)")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()
