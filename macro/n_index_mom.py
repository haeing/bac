import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle

# -----------------------
# Fixed parameters
# -----------------------
n_aero = 1.115 
pmin, pmax = 250.0, 1500.0         # momentum range [MeV/c]
sep_p1, sep_p2 = 625.0, 970.0    # separation region [MeV/c]

ellipse_height = 0.06            # height in refractive index
ellipse_alpha = 0.25

# particle masses [MeV/c^2]
m_pi = 139.57039
m_K  = 493.677

def n_threshold(p, m):
    """
    Cherenkov threshold:
    n_thr(p) = 1/beta = sqrt(1 + (m/p)^2)
    p, m in MeV units
    """
    p = np.asarray(p)
    return np.sqrt(1.0 + (m/p)**2)

# momentum grid (avoid zero)
p = np.linspace(20.0, pmax, 3000)

# threshold curves
n_pi = n_threshold(p, m_pi)
n_K  = n_threshold(p, m_K)

# -----------------------
# Plot
# -----------------------
fig, ax = plt.subplots(figsize=(8.5, 5.2))

ax.plot(p, n_pi, lw=2, label=r'$\pi$ threshold')
ax.plot(p, n_K,  lw=2, label=r'$K$ threshold')

# BAC refractive index line
ax.axhline(
    n_aero,
    color='black',
    ls='--',
    lw=2.5,
    label='Silica aerogel (n = 1.115)'
)

# separation ellipse

rect_height = 0.01     # 충분히 얇음 (0.005~0.01 추천)

rect = Rectangle(
    (sep_p1, n_aero - rect_height / 2),
    width=sep_p2 - sep_p1,
    height=rect_height,
    facecolor='red',
    alpha=0.25,
    edgecolor='none'
)

ax.add_patch(rect)

# vertical guide lines
ax.axvline(sep_p1, color='black', ls=':', lw=2)
ax.axvline(sep_p2, color='black', ls=':', lw=2)

# axis & style
ax.set_xlim(pmin, pmax)
ax.set_ylim(1.0, 1.35) 
ax.set_xlabel(r'Momentum [MeV/$c$]')
ax.set_ylabel('Refractive index')

ax.grid(True, alpha=0.3)
ax.legend(frameon=False, loc='upper right')


plt.tight_layout()
plt.savefig("pi_k_threshold.pdf")
plt.show()


