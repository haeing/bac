import numpy as np
import matplotlib.pyplot as plt

# ---- constants / settings ----
alpha = 1/137.035999084
n = 1.115
z = 1
m_pi = 139.57039  # MeV/c^2

# wavelength range [nm] -> [cm]
lam1_nm, lam2_nm = 300.0, 900.0
lam1_cm = lam1_nm * 1e-7
lam2_cm = lam2_nm * 1e-7

# ∫(1/λ^2) dλ = (1/λ1 - 1/λ2)  (λ in cm)
int_1_over_lam2 = (1/lam1_cm - 1/lam2_cm)

def beta_from_p(p_MeVc, m_MeV):
    return p_MeVc / np.sqrt(p_MeVc**2 + m_MeV**2)

def photons_per_cm_pion(p_MeVc):
    beta = beta_from_p(p_MeVc, m_pi)
    factor = 1.0 - 1.0 / (n**2 * beta**2)
    factor = np.where(factor > 0, factor, 0.0)  # below threshold -> 0
    return (2*np.pi*alpha*(z**2)) * factor * int_1_over_lam2  # photons/cm

# momentum grid
p = np.linspace(620.0, 970.0, 400)  # MeV/c
N_per_cm = photons_per_cm_pion(p)

plt.figure(figsize=(8,4))
plt.plot(p, N_per_cm, linewidth=2, label=r"$\pi^\pm$")
plt.xlabel(r"Momentum $p$ (MeV/$c$)")
plt.ylabel(r"Photons per cm")
plt.title(r"Frank–Tamm yield ($n=1.115$, $\lambda=300$–$900$ nm, $z=1$)")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("pion_cherenkov_photons_per_cm_vs_p.pdf", bbox_inches="tight")
plt.savefig("pion_cherenkov_photons_per_cm_vs_p.png", dpi=300, bbox_inches="tight")
plt.show()
