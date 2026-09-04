import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


# ============================================================
# Global style
# ============================================================
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 13,
    "axes.labelsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "axes.linewidth": 1.0,
})


rng = np.random.default_rng(12)


# ============================================================
# Basic parameters
# ============================================================
beam_width = 100.0       # beam-arrival window [ns]
tdc_threshold = -0.18    # discriminator threshold [V]
dt = 0.2                 # time step [ns]


# Irregular positions of beam windows
beam_starts = np.array([
    20.0,
    180.0,
    355.0,
    535.0
])

beam_ends = beam_starts + beam_width


t_min = 0.0
t_max = beam_ends[-1] + 30.0

t = np.arange(
    t_min,
    t_max + dt,
    dt
)


# ============================================================
# BAC pulse shape
# ============================================================
def bac_pulse(t, t0, amplitude, width):
    """
    Schematic negative-polarity BAC pulse.

    Timing ratio:

        rise : peak : total width
         15  :  40  : 100

    The total width can vary between pulses, while the
    relative timing structure is preserved.
    """

    tau = t - t0

    pulse = np.zeros_like(t)

    # Preserve 15 : 40 : 100 ratio
    rise_time = 0.15 * width
    peak_time = 0.40 * width


    # --------------------------------------------------------
    # Initial rise
    # --------------------------------------------------------
    m1 = (
        (tau >= 0.0)
        & (tau < rise_time)
    )

    if np.any(m1):

        x = tau[m1] / rise_time

        pulse[m1] = (
            -0.35 * amplitude
            * (0.5 - 0.5 * np.cos(np.pi * x))
        )


    # --------------------------------------------------------
    # Approach to peak
    # --------------------------------------------------------
    m2 = (
        (tau >= rise_time)
        & (tau < peak_time)
    )

    if np.any(m2):

        x = (
            (tau[m2] - rise_time)
            / (peak_time - rise_time)
        )

        pulse[m2] = (
            -amplitude
            * (
                0.35
                + 0.65
                * (0.5 - 0.5 * np.cos(np.pi * x))
            )
        )


    # --------------------------------------------------------
    # Recovery to baseline
    # --------------------------------------------------------
    m3 = (
        (tau >= peak_time)
        & (tau <= width)
    )

    if np.any(m3):

        x = (
            (tau[m3] - peak_time)
            / (width - peak_time)
        )

        pulse[m3] = (
            -amplitude
            * (0.5 + 0.5 * np.cos(np.pi * x))
        )


    return pulse


# ============================================================
# Pulse configuration
#
# (relative start [ns], amplitude [V], width [ns])
#
# Window 1 : FAIL
# Window 2 : PASS
# Window 3 : FAIL
# Window 4 : PASS
#
# All pulses are contained within the beam window.
# ============================================================
window_pulses = [

    # Window 1
    # Small pulse
    [
        (27, 0.055, 45),
    ],


    # Window 2
    # Large pulse -> PASS
    [
        (10, 0.245, 80),
    ],


    # Window 3
    # Two smaller pulses
    [
        (7,  0.070, 32),
        (54, 0.125, 38),
    ],


    # Window 4
    # Small + large pulse -> PASS
    [
        (5,  0.055, 23),
        (32, 0.285, 62),
    ],
]


# ============================================================
# Construct BAC waveform
# ============================================================
signal = np.zeros_like(t)


# Baseline noise
noise = rng.normal(
    loc=0.0,
    scale=0.0035,
    size=len(t)
)

noise = np.convolve(
    noise,
    np.ones(5) / 5,
    mode="same"
)

signal += noise


# Add pulses
for iw, pulses in enumerate(window_pulses):

    window_start = beam_starts[iw]

    for rel_start, amplitude, width in pulses:

        pulse_start = (
            window_start
            + rel_start
        )

        signal += bac_pulse(
            t,
            pulse_start,
            amplitude,
            width
        )


# ============================================================
# ADC integration and TDC response
# ============================================================
adc_integrals = []
tdc_hits = []
tdc_signals = []


for iw in range(len(beam_starts)):

    low = beam_starts[iw]
    high = beam_ends[iw]

    mask = (
        (t >= low)
        & (t <= high)
    )


    # --------------------------------------------------------
    # ADC integration
    # --------------------------------------------------------
    adc = -np.trapz(
        signal[mask],
        t[mask]
    )

    adc_integrals.append(adc)


    # --------------------------------------------------------
    # TDC digital response
    #
    # Negative BAC signal:
    #
    # BAC < threshold  --> TDC = 1
    #
    # TDC remains high for the entire time during which
    # the analog signal is beyond threshold.
    # --------------------------------------------------------
    tdc_signal = np.zeros_like(t)

    crossing = (
        mask
        & (signal < tdc_threshold)
    )

    tdc_signal[crossing] = 1.0

    tdc_signals.append(
        tdc_signal
    )

    tdc_hits.append(
        int(np.any(crossing))
    )


adc_integrals = np.array(
    adc_integrals
)

tdc_hits = np.array(
    tdc_hits
)


# ============================================================
# Combined TDC signal
# ============================================================
tdc_combined = np.zeros_like(t)

for tdc_signal in tdc_signals:

    tdc_combined = np.maximum(
        tdc_combined,
        tdc_signal
    )


# ============================================================
# Figure
# ============================================================
fig, (ax1, ax2) = plt.subplots(
    2,
    1,
    figsize=(11, 5.8),
    sharex=True,
    gridspec_kw={
        "height_ratios": [2.7, 1.0],
        "hspace": 0.12
    }
)


# ============================================================
# Beam-window backgrounds
# ============================================================
for low, high in zip(
    beam_starts,
    beam_ends
):

    ax1.axvspan(
        low,
        high,
        facecolor="0.94",
        edgecolor="none",
        zorder=0
    )

    ax2.axvspan(
        low,
        high,
        facecolor="0.94",
        edgecolor="none",
        zorder=0
    )


# ============================================================
# BAC waveform
# ============================================================
ax1.plot(
    t,
    signal,
    color="black",
    lw=1.3,
    zorder=4
)


# ============================================================
# TDC threshold
# ============================================================
ax1.axhline(
    tdc_threshold,
    color="black",
    linestyle="--",
    lw=1.2,
    zorder=3
)

ax1.text(
    1.008,
    tdc_threshold,
    "TDC threshold",
    transform=ax1.get_yaxis_transform(),
    ha="left",
    va="center",
    fontsize=12
)


# ============================================================
# ADC-integrated waveform
#
# The hatched area indicates the BAC signal used for
# determining the ADC value of each event.
# ============================================================
for iw in range(len(beam_starts)):

    low = beam_starts[iw]
    high = beam_ends[iw]

    mask = (
        (t >= low)
        & (t <= high)
    )

    ax1.fill_between(
        t[mask],
        signal[mask],
        0,
        facecolor="none",
        edgecolor="0.45",
        hatch="////",
        linewidth=0.0,
        zorder=2
    )


# ============================================================
# Beam-window label
#
# Show only once.
# ============================================================
low = beam_starts[0]
high = beam_ends[0]

center = 0.5 * (
    low + high
)

ax1.annotate(
    "",
    xy=(low, 0.037),
    xytext=(high, 0.037),
    arrowprops=dict(
        arrowstyle="<->",
        color="black",
        lw=1.1
    )
)

ax1.text(
    center,
    0.047,
    "Beam window",
    ha="center",
    va="bottom",
    fontsize=11
)


# ============================================================
# BAC axis
# ============================================================
ax1.set_ylim(
    -0.34,
    0.075
)

# Standard vertical y-axis title,
# as in the example figure.
ax1.set_ylabel(
    "BAC signal [V]",
    labelpad=10
)


# No x tick labels on upper panel
ax1.tick_params(
    axis="x",
    labelbottom=False
)


# ============================================================
# TDC digital signal
# ============================================================
ax2.plot(
    t,
    tdc_combined,
    color="red",
    lw=1.8,
    drawstyle="steps-post",
    zorder=3
)

ax2.axhline(
    0,
    color="black",
    lw=1.0,
    zorder=1
)


# ============================================================
# TDC axis
# ============================================================
ax2.set_ylim(
    -0.15,
    1.25
)

ax2.set_yticks([
    0,
    1
])

# Standard vertical y-axis title
ax2.set_ylabel(
    "TDC hit",
    labelpad=10
)


# ============================================================
# X axis
#
# Use a normal xlabel and give it sufficient separation
# from the numerical tick labels.
# ============================================================
ax2.set_xlabel(
    "Time [ns]",
    labelpad=12
)

# Move xlabel toward right side, as in the paper style
ax2.xaxis.set_label_coords(
    1.0,
    -0.16
)

ax2.xaxis.label.set_horizontalalignment(
    "right"
)


# ============================================================
# Legend
# ============================================================

legend_elements = [

    Patch(
        facecolor="none",
        edgecolor="0.45",
        hatch="////",
        label="Total: all events"
    ),

    Line2D(
        [0],
        [0],
        color="red",
        lw=1.8,
        label="Pass: events with a TDC hit"
    ),
]


ax1.legend(
    handles=legend_elements,
    loc="lower left",
    bbox_to_anchor=(0.005, 0.015),
    frameon=False,
    fontsize=11
)


# ============================================================
# General cosmetics
# ============================================================
for ax in (ax1, ax2):

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        direction="out",
        length=5,
        width=1
    )


ax2.set_xlim(
    t_min,
    t_max
)


# ============================================================
# Layout
#
# Extra bottom margin prevents "Time [ns]" from overlapping
# with the x-axis numbers.
# ============================================================
plt.subplots_adjust(
    left=0.12,
    right=0.86,
    top=0.95,
    bottom=0.19,
    hspace=0.12
)


# ============================================================
# Save
# ============================================================
plt.savefig(
    "BAC_threshold_schematic.pdf",
    bbox_inches="tight"
)

plt.savefig(
    "BAC_threshold_schematic.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# Check
# ============================================================
print()
print("Beam-window summary")
print("----------------------------------------")
print("Window    ADC integral    TDC")
print("----------------------------------------")

for i in range(len(beam_starts)):

    print(
        f"{i+1:3d}"
        f"       {adc_integrals[i]:8.3f}"
        f"       {tdc_hits[i]}"
    )
