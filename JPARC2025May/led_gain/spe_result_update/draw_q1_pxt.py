#!/usr/bin/env python3
"""
Draw Q1 (blue) and PXT (red) vs LED_HV
from spe_summary.csv (no pandas required).

- Filter: mppc_hv == 58
- One PDF
- One page per MPPC board
- Each page: 4x4 grid (ch 0..15)
"""

import csv
import os
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

CSV_PATH = "spe_summary.csv"
TARGET_MPPC_HV = 58.0
TOL = 1e-6

NCH = 16
NROWS, NCOLS = 4, 4

def main():
    if not os.path.isfile(CSV_PATH):
        raise FileNotFoundError(f"CSV not found: {CSV_PATH}")

    # --------------------------------------------------
    # Read CSV (standard library only)
    # data[board][ch] = list of (led_hv, q1, pxt)
    # --------------------------------------------------
    data = defaultdict(lambda: defaultdict(list))

    with open(CSV_PATH, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                mppc_hv = float(row["mppc_hv"])
            except Exception:
                continue

            if abs(mppc_hv - TARGET_MPPC_HV) > TOL:
                continue

            try:
                board = int(row["mppc_board"])
                ch    = int(row["mppc_ch"])
                led_hv = float(row["led_hv"])
                q1  = float(row["q1"])
                pxt = float(row["pxt"])
            except Exception:
                continue

            data[board][ch].append((led_hv, q1, pxt))

    if not data:
        raise RuntimeError(f"No entries with mppc_hv == {TARGET_MPPC_HV}")

    boards = sorted(data.keys())

    out_pdf = f"q1_pxt_vs_ledhv_mppcHV{int(TARGET_MPPC_HV)}_allboards.pdf"

    # --------------------------------------------------
    # Plot
    # --------------------------------------------------
    with PdfPages(out_pdf) as pdf:
        for board in boards:
            fig, axes = plt.subplots(
                NROWS, NCOLS,
                figsize=(12, 10),
                sharex=True,
                sharey=False
            )
            fig.suptitle(
                f"MPPC board {board}  |  MPPC_HV = {TARGET_MPPC_HV} V",
                fontsize=14
            )

            for ch in range(NCH):
                ax = axes[ch // NCOLS][ch % NCOLS]
                ax.set_title(f"ch {ch}", fontsize=10)

                if ch not in data[board]:
                    ax.text(
                        0.5, 0.5, "no data",
                        ha="center", va="center",
                        transform=ax.transAxes,
                        fontsize=9
                    )
                    ax.grid(True, alpha=0.3)
                    continue

                # sort by LED_HV
                rows = sorted(data[board][ch], key=lambda x: x[0])
                x = [r[0] for r in rows]
                q1 = [r[1] for r in rows]
                pxt = [r[2] for r in rows]

                ax.plot(x, q1, "o-", color="blue", label="Q1")
                ax.plot(x, pxt, "s-", color="red", label="PXT")

                ax.grid(True, alpha=0.3)

                if ch == 0:
                    ax.legend(fontsize=9)

            # axis labels
            for ax in axes[-1]:
                ax.set_xlabel("LED HV [V]")
            for r in range(NROWS):
                axes[r][0].set_ylabel("Value")

            fig.tight_layout(rect=[0, 0.02, 1, 0.95])
            pdf.savefig(fig)
            plt.close(fig)

    print(f"[OK] PDF written: {out_pdf}")

if __name__ == "__main__":
    main()
