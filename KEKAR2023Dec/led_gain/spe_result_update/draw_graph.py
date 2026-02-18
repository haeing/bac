import pandas as pd
import matplotlib.pyplot as plt

CSV = "Q1_pxt_all.csv"

df = pd.read_csv(CSV)

# 안전하게 타입 보정
df["board"] = df["board"].astype(int)
df["channel"] = df["channel"].astype(int)

for b in sorted(df["board"].unique()):
    d = df[df["board"] == b].sort_values("channel")

    fig, ax1 = plt.subplots(figsize=(7, 4.5))

    # Gain(Q1): left axis
    ax1.errorbar(
        d["channel"], d["Q1"],
        yerr=d["Q1_err"],
        fmt="o", capsize=3
    )
    ax1.set_xlabel("channel (0-origin)")
    ax1.set_ylabel("Gain (Q1)")
    ax1.set_title(f"Board {b}: Gain(Q1) and pxt vs channel")
    ax1.set_xticks(range(0, 16))

    # pxt: right axis
    ax2 = ax1.twinx()
    ax2.errorbar(
        d["channel"], d["pxt"],
        yerr=d["pxt_err"],
        fmt="s", capsize=3,
        color="red", ecolor="red"  
    )
    ax2.set_ylabel("pxt")

    fig.tight_layout()
    #fig.savefig(f"board{b}_Q1_pxt.png", dpi=200)
    fig.savefig(f"board{b}_Q1_pxt.pdf")
    plt.close(fig)

