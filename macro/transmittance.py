import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# =========================
# Input
# =========================
fname = "ag22k001.txt"
t_mm = 12.4                 # sample thickness [mm]
lam_eval_nm = 400.0         # wavelength for A_T

# thickness for transmission length
d_cm = t_mm / 10.0          # 10 mm = 1 cm

# =========================
# Load data
# =========================
data = np.loadtxt(fname)
lam_nm = data[:, 0].astype(float)
T = data[:, 1].astype(float) / 100.0   # transmittance (fraction)

# =========================
# Rayleigh model
# T = A * exp(-B / lambda^4),  B = C*t
# =========================
def model(lam_nm, A, B):
    return A * np.exp(-B / (lam_nm**4))

mask = (T > 0) & (T < 1.0)

A0 = np.max(T)
B0 = 1e10

popt, pcov = curve_fit(
    model,
    lam_nm[mask],
    T[mask],
    p0=[A0, B0],
    maxfev=20000
)

A_fit, B_fit = popt
A_err, B_err = np.sqrt(np.diag(pcov))

# =========================
# Rayleigh coefficient C
# =========================
C_nm4_per_mm = B_fit / t_mm
C_nm4_per_mm_err = B_err / t_mm

# unit conversion: nm^4/mm -> um^4/cm
C_um4_per_cm = C_nm4_per_mm * 1e-11
C_um4_per_cm_err = C_nm4_per_mm_err * 1e-11

# =========================
# Transmittance at 400 nm
# =========================
T400 = model(lam_eval_nm, A_fit, B_fit)

# =========================
# Transmission length A_T(400)
# A_T = -d / ln T
# =========================
if 0.0 < T400 < 1.0:
    AT400_cm = -d_cm / np.log(T400)
else:
    AT400_cm = np.nan

# =========================
# Print results (paper-ready)
# =========================
print("==== Rayleigh scattering fit (paper-ready) ====")
print(f"A = {A_fit:.6g} ± {A_err:.2g}")
print(f"C = ({C_um4_per_cm:.4g} ± {C_um4_per_cm_err:.2g})  μm^4 / cm")
print()
print(f"T({lam_eval_nm:.0f} nm) = {T400:.6f}  ({T400*100:.3f} %)")
print(f"A_T({lam_eval_nm:.0f} nm) = {AT400_cm:.6g} cm")
print(f"(using thickness d = {d_cm:.3f} cm = {t_mm} mm)")

# =========================
# Plot
# =========================
lam_dense = np.linspace(lam_nm.min(), lam_nm.max(), 2000)

plt.figure(figsize=(7.5, 4.5))
plt.plot(lam_nm, T, "o", ms=4, label="Data")
plt.plot(
    lam_dense,
    model(lam_dense, A_fit, B_fit),
    "-",
    lw=2,
    label=r"Fit: $A\,\exp(-Ct/\lambda^4)$"
)


plt.xlabel("Wavelength [nm]",fontsize=16)
plt.ylabel("Transmittance",fontsize=16)
plt.grid(alpha=0.3)
#plt.legend()
plt.tight_layout()
plt.savefig("ag22k001_rayleigh_fit.png", dpi=300, bbox_inches="tight")
plt.savefig("ag22k001_rayleigh_fit.pdf", bbox_inches="tight")
plt.show()
