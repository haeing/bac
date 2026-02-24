#!/usr/bin/env python3
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

# y-range margin (10% of span). If span is too small, use absolute floor margins.
REL_MARGIN = 0.12
MIN_ABS_MARGIN_Q1 = 0.05   # Q1 units
MIN_ABS_MARGIN_PXT = 0.01  # PXT units

def padded_ylim(y, yerr, rel_margin, abs_floor):
    """
    Compute y-limits including error bars, then add margins.
    y, yerr: list of floats (same length)
    """
    vals = []
    for yi, ei in zip(y, yerr):
        if yi is None:
            continue
        if ei is None:
            ei = 0.0
        vals.append(yi - abs(ei))
        vals.append(yi + abs(ei))

    if not vals:
        return None

    ymin = min(vals)
    ymax = max(vals)

    span = ymax - ymin
    pad = max(abs_floor, rel_margin * span) if span > 0 else abs_floor

    return (ymin - pad, ymax + pad)

def main():
    if not os.path.isfile(CSV_PATH):
        raise FileNotFoundError(f"CSV not found: {CSV_PATH}")

    # data[board][ch] = list of (led_hv, q1, q1e, pxt, pxte)
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

            # errors can be empty in CSV; treat as 0
            def get_err(key):
                try:
                    v = row.get(key, "")
                    if v is None or str(v).strip() == "":
                        return 0.0
                    return float(v)
                except Exception:
                    return 0.0

            q1e  = get_err("q1_err")
            pxte = get_err("pxt_err")

            data[board][ch].append((led_hv, q1, q1e, pxt, pxte))

    if not data:
        raise RuntimeError(f"No entries with mppc_hv == {TARGET_MPPC_HV}")

    boards = sorted(data.keys())
    out_pdf = f"q1_left_pxt_right_vs_ledhv_mppcHV{int(TARGET_MPPC_HV)}_allboards.pdf"

    with PdfPages(out_pdf) as pdf:
        for board in boards:
            fig, axes = plt.subplots(
                NROWS, NCOLS,
                figsize=(12, 10),
                sharex=True,
                sharey=False
            )
            fig.suptitle(f"MPPC board {board}  |  MPPC_HV = {TARGET_MPPC_HV} V", fontsize=14)

            for ch in range(NCH):
                ax = axes[ch // NCOLS][ch % NCOLS]
                ax.set_title(f"ch {ch}", fontsize=10)

                if ch not in data[board]:
                    ax.text(0.5, 0.5, "no data", ha="center", va="center",
                            transform=ax.transAxes, fontsize=9)
                    ax.grid(True, alpha=0.3)
                    continue

                rows = sorted(data[board][ch], key=lambda x: x[0])

                x    = [r[0] for r in rows]
                q1   = [r[1] for r in rows]
                q1e  = [r[2] for r in rows]
                pxt  = [r[3] for r in rows]
                pxte = [r[4] for r in rows]

                # left axis: Q1 with error bars
                l1 = ax.errorbar(
                    x, q1, yerr=q1e,
                    fmt="o-", color="blue",
                    elinewidth=1, capsize=2,
                    label="Q1"
                )
                ax.set_ylabel("Q1", color="blue")
                ax.tick_params(axis="y", labelcolor="blue")
                ax.grid(True, alpha=0.3)

                # right axis: PXT with error bars
                ax2 = ax.twinx()
                l2 = ax2.errorbar(
                    x, pxt, yerr=pxte,
                    fmt="s-", color="red",
                    elinewidth=1, capsize=2,
                    label="PXT"
                )
                ax2.set_ylabel("PXT", color="red")
                ax2.tick_params(axis="y", labelcolor="red")

                # set padded y-lims (include errors + margins)
                yl1 = padded_ylim(q1, q1e, REL_MARGIN, MIN_ABS_MARGIN_Q1)
                yl2 = padded_ylim(pxt, pxte, REL_MARGIN, MIN_ABS_MARGIN_PXT)
                if yl1:
                    ax.set_ylim(*yl1)
                if yl2:
                    ax2.set_ylim(*yl2)

                # legend only once to reduce clutter
                if ch == 0:
                    # errorbar returns container; .lines[0] gives the Line2D for legend
                    ax.legend(handles=[l1.lines[0], l2.lines[0]], loc="best", fontsize=9)

            for ax in axes[-1]:
                ax.set_xlabel("LED HV [V]")

            fig.tight_layout(rect=[0, 0.02, 1, 0.95])
            pdf.savefig(fig)
            plt.close(fig)

    print(f"[OK] PDF written: {out_pdf}")

if __name__ == "__main__":
    main()
