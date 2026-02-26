import pandas as pd
import matplotlib.pyplot as plt

HV_SCAN_CSV = "Q1_pxt_all_kek.csv"   # 다른 셋업 (HV scan)
BASELINE_CSV = "Q1_pxt_all.csv"      # 원래 셋업

# --- read ---
hv = pd.read_csv(HV_SCAN_CSV)
hv["board"] = hv["board"].astype(int)
hv["channel"] = hv["channel"].astype(int)
hv["HV"] = hv["HV"].astype(int)

# ch=1, HV=57만
hv57 = hv[(hv["channel"] == 1) & (hv["HV"] == 58)].sort_values("board")

base = pd.read_csv(BASELINE_CSV)
base["board"] = base["board"].astype(int)
base["channel"] = base["channel"].astype(int)

# ch=1만
base = base[base["channel"] == 1].sort_values("board")

# --- plot ---
fig, ax1 = plt.subplots(figsize=(7.5, 4.8))

# Q1 (left axis)
ax1.errorbar(
    base["board"], base["Q1"], yerr=base["Q1_err"],
    fmt="o", capsize=3, label="Baseline (57 V) Q1"
)
ax1.errorbar(
    hv57["board"], hv57["Q1"], yerr=hv57["Q1_err"],
    fmt="s", capsize=3, label="HV-scan setup (57 V) Q1"
)

ax1.set_xlabel("board")
ax1.set_xticks(sorted(base["board"].unique()))
ax1.set_ylabel("Gain (Q1)")
ax1.set_title("Comparison at 57 V (ch = 1)")

# pxt (right axis)
ax2 = ax1.twinx()
ax2.errorbar(
    base["board"], base["pxt"], yerr=base["pxt_err"],
    fmt="o", capsize=3, color="red", ecolor="red",
    label="Baseline (57 V) pxt"
)
ax2.errorbar(
    hv57["board"], hv57["pxt"], yerr=hv57["pxt_err"],
    fmt="s", capsize=3, color="orange", ecolor="orange",
    label="HV-scan setup (57 V) pxt"
)
ax2.set_ylabel("pxt")

# legend (merge both axes)
h1, l1 = ax1.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=9)

fig.tight_layout()
#fig.savefig("Compare_57V_baseline_vs_HVscan_ch1_Q1_pxt.png", dpi=200)
fig.savefig("Compare_57V_baseline_vs_HVscan_ch1_Q1_pxt.pdf")
plt.close(fig)

print("Saved: Compare_57V_baseline_vs_HVscan_ch1_Q1_pxt.(png/pdf)")
