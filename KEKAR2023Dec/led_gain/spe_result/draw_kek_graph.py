import pandas as pd
import matplotlib.pyplot as plt

HV_SCAN_CSV = "Q1_pxt_all_kek.csv"   # columns: board,channel,HV,Q1,Q1_err,pxt,pxt_err,...
BASELINE_CSV = "Q1_pxt_all.csv"      # columns: board,channel,Q1,Q1_err,pxt,pxt_err,...

# -----------------------------
# 1) HV scan: per-board plots
# -----------------------------
hv = pd.read_csv(HV_SCAN_CSV)
hv["board"] = hv["board"].astype(int)
hv["channel"] = hv["channel"].astype(int)
hv["HV"] = hv["HV"].astype(float)

# channel=1만 사용
hv = hv[hv["channel"] == 1].copy()

for b in sorted(hv["board"].unique()):
    d = hv[hv["board"] == b].sort_values("HV")

    fig, ax1 = plt.subplots(figsize=(7, 4.5))

    # Q1 (Gain) — left axis
    ax1.errorbar(
        d["HV"], d["Q1"], yerr=d["Q1_err"],
        fmt="o", capsize=3
    )
    ax1.set_xlabel("HV (V)")
    ax1.set_ylabel("Gain (Q1)")
    ax1.set_title(f"Board {b}, ch=1: Q1 and pxt vs HV")

    # pxt — right axis (red)
    ax2 = ax1.twinx()
    ax2.errorbar(
        d["HV"], d["pxt"], yerr=d["pxt_err"],
        fmt="s", capsize=3, color="red", ecolor="red"
    )
    ax2.set_ylabel("pxt")

    fig.tight_layout()
    #fig.savefig(f"HVscan_board{b}_ch1_Q1_pxt.png", dpi=200)
    fig.savefig(f"HVscan_board{b}_ch1_Q1_pxt.pdf")
    plt.close(fig)

print("Saved 4 plots: HVscan_board{0..3}_ch1_Q1_pxt.(png/pdf)")

# ----------------------------------------
# 2) Baseline(원래 전체 결과): board-plot
#    board별 channel=1만, x축=board
# ----------------------------------------
base = pd.read_csv(BASELINE_CSV)
base["board"] = base["board"].astype(int)
base["channel"] = base["channel"].astype(int)
base = base[base["channel"] == 1].sort_values("board").copy()

fig, ax1 = plt.subplots(figsize=(7, 4.5))

ax1.errorbar(
    base["board"], base["Q1"], yerr=base["Q1_err"],
    fmt="o", capsize=3
)
ax1.set_xlabel("board")
ax1.set_xticks(sorted(base["board"].unique()))
ax1.set_ylabel("Gain (Q1)")
ax1.set_title("Baseline (all boards, ch=1): Q1 and pxt vs board")

ax2 = ax1.twinx()
ax2.errorbar(
    base["board"], base["pxt"], yerr=base["pxt_err"],
    fmt="s", capsize=3, color="red", ecolor="red"
)
ax2.set_ylabel("pxt")

fig.tight_layout()
#fig.savefig("Baseline_allBoards_ch1_Q1_pxt.png", dpi=200)
fig.savefig("Baseline_allBoards_ch1_Q1_pxt.pdf")
plt.close(fig)

print("Saved 5th plot: Baseline_allBoards_ch1_Q1_pxt.(png/pdf)")
