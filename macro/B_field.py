import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# -----------------------
# file name
# -----------------------
filename = "ShsFieldMap_E42_20190903_1T"

# -----------------------
# load data
# columns: x y z Bx By Bz   (x,y,z in cm)
# -----------------------
data = np.loadtxt(filename, usecols=(0, 1, 2, 3, 4, 5))

# convert position: cm → mm
x = data[:, 0] * 10.0
y = data[:, 1] * 10.0
z = data[:, 2] * 10.0

By = data[:, 4]   # Tesla

# -----------------------
# selection
# -----------------------

# z = 0 plane (after unit conversion → still 0)
mask_z0 = np.isclose(z, 0.0, atol=1e-6)

# spatial window (mm)
mask_xy = (
    (x >= -56.25) & (x <= 90.72) &
    (y >= -513.0) & (y <= -393.0)
)

mask = mask_z0 & mask_xy

x_sel  = x[mask]
y_sel  = y[mask]
By_sel = By[mask]

print(f"Selected points: {len(x_sel)}")

# -----------------------
# 3D plot: x, y, By
# -----------------------
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

sc = ax.scatter(
    x_sel, y_sel, By_sel,
    c=By_sel,
    cmap='viridis',
    s=20
)

ax.set_xlabel('x [mm]')
ax.set_ylabel('z [mm]')
ax.set_zlabel(r'$B_y$ [T]')

cbar = plt.colorbar(sc, pad=0.12)
cbar.set_label(r'$B_y$ [T]')

ax.set_title(r'$B_y(x, y)$ at $z = 0$ mm')

plt.tight_layout()
plt.show()
