import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# =========================================
# Style
# =========================================

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "stix"

# =========================================
# ROOT histogram
# =========================================

root_file = "t110_graph_344_k_pi.root"
hist_total_name = "hist_bac_btof"
hist_pass_name  = "hist_bac_btof_pass"


with uproot.open(root_file) as f:
    h_total = f[hist_total_name]
    h_pass = f[hist_pass_name]

    values_total, xedges, yedges = h_total.to_numpy()
    values_pass,  _,      _      = h_pass.to_numpy()
values_kaon = values_total - values_pass
values_kaon[values_kaon<0] = 0


# =========================================
# Projections
# =========================================

xproj_total = values_total.sum(axis=1)
yproj_total = values_total.sum(axis=0)

xproj_pion = values_pass.sum(axis=1)
xproj_kaon = values_kaon.sum(axis=1)

xcenters = 0.5 * (xedges[:-1] + xedges[1:])
ycenters = 0.5 * (yedges[:-1] + yedges[1:])


# =========================================
# Mask zero entries -> white
# =========================================

values_masked = np.ma.masked_where(values_total.T <= 0, values_total.T)

# ROOT-like palette
cmap = plt.get_cmap("turbo").copy()
cmap.set_bad("white")

# =========================================
# Figure
# =========================================

fig = plt.figure(figsize=(10, 8))

gs = GridSpec(
    2,
    3,
    width_ratios=[1.3, 4.0, 0.25],
    height_ratios=[1.3, 4.0],
    wspace=0.0,
    hspace=0.0
)

ax_empty = fig.add_subplot(gs[0, 0])

ax_xproj = fig.add_subplot(gs[0, 1])

ax_yproj = fig.add_subplot(gs[1, 0])

ax_2d = fig.add_subplot(
    gs[1, 1],
    sharex=ax_xproj,
    sharey=ax_yproj
)

ax_cbar = fig.add_subplot(gs[1, 2])

ax_empty.axis("off")

# =========================================
# Main 2D histogram
# =========================================

mesh = ax_2d.pcolormesh(
    xedges,
    yedges,
    values_masked,
    cmap=cmap,
    shading="auto",
    vmin=0
)

# =========================================
# Color bar
# =========================================

cbar = fig.colorbar(mesh, cax=ax_cbar)

cbar.set_label(
    "Counts",
    fontsize=16
)

# move colorbar title to center
cbar.ax.yaxis.set_label_coords(2.5, 0.5)

# =========================================
# X projection
# =========================================

ax_xproj.step(
    xcenters,
    xproj_pion,
    where="mid",
    color="blue",
    linewidth=2,
    label="Pion"
)

ax_xproj.step(
    xcenters,
    xproj_kaon,
    where="mid",
    color="red",
    linewidth=2,
    label="Kaon"
)



# make bottom touch exactly at 0
ax_xproj.set_xlim(xedges[0], xedges[-1])
ax_xproj.set_ylim(0, xproj_pion.max() * 1.05)
ax_xproj.set_ylabel("")

ax_xproj.tick_params(labelleft=False)

ax_xproj.set_xlim(
    xedges[0],
    xedges[-1]
)

ax_xproj.tick_params(labelbottom=False)


# =========================================
# Y projection
# =========================================

ax_yproj.step(
    yproj_total,
    ycenters,
    where="mid",
    color="black",
    linewidth=2
)

# reverse counts axis
ax_yproj.set_xlim(
    yproj_total.max() * 1.05,
    0
)

ax_yproj.set_ylim(
    yedges[0],
    yedges[-1]
)


ax_yproj.set_ylabel(
    "Time of Flight [ns]",
    fontsize=16
)

# move y title to top
ax_yproj.yaxis.set_label_coords(-0.2, 0.8)


ax_yproj.tick_params(labelbottom=False)

# =========================================
# Main 2D axis
# =========================================

ax_2d.set_xlim(
    xedges[0],
    xedges[-1]
)

ax_2d.set_ylim(
    yedges[0],
    yedges[-1]
)

ax_2d.set_xlabel(
    r"$\mathrm{N}_\mathrm{p.e.}$",
    fontsize=16
)

# move x title to right edge
ax_2d.xaxis.set_label_coords(0.95, -0.08)

# remove y labels from main histogram
ax_2d.tick_params(labelleft=False)


# Put threshold line
cut_x = 15


# 2D histogram

ax_2d.axvline(
    cut_x,
    color="grey",
    linewidth=1,
    linestyle="--"
)

# X projection

ax_xproj.axvline(
    cut_x,
    color="grey",
    linewidth=1,
    linestyle="--"
)

ax_2d.text(
    cut_x-1.5,
    -2.1,
    "15",
    color="gray",
    fontsize=16,
    ha="center",
    va="top"
)

# pion region

ax_2d.axhspan(
    -0.9,
    1,
    color="deepskyblue",
    alpha=0.10
)

ax_yproj.axhspan(
    -0.9,
    1,
    color="deepskyblue",
    alpha=0.10
)

# kaon region

ax_2d.axhspan(
    4.5,
    5.6,
    color="red",
    alpha=0.08
)

ax_yproj.axhspan(
    4.5,
    5.6,
    color="red",
    alpha=0.08
)

# Put Text
# pion label

ax_2d.text(
    -4.5,
    6.1,
    "Kaon",
    color="red",
    fontsize=16,
    fontweight="bold"
)

# kaon label

ax_2d.text(
    72,
    -1.5,
    "Pion",
    color="blue",
    fontsize=16,
    fontweight="bold"
)

# =========================================
# ROOT-like ticks
# =========================================

for ax in [ax_xproj, ax_yproj, ax_2d, ax_cbar]:

    ax.tick_params(
        direction="in",
        top=True,
        right=True,
        length=12,
        width=1.8,
        labelsize=16
    )

    ax.tick_params(
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=6,
        width=1.2
    )

    ax.minorticks_on()

# thicker frame
for ax in [ax_xproj, ax_yproj, ax_2d]:
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)

plt.savefig(
    "compare_k_pi.pdf",
    bbox_inches="tight"
)
plt.show()
