#!/usr/bin/env python3
import os, glob
from collections import defaultdict

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# -------------------------
# Config
# -------------------------
INFO_FILE = "../led_run_info.txt"
OUT_PDF   = "spe_summary_all.pdf"

# ROOT files live in current dir (edit if needed)
ROOT_DIR  = "."
FNAME_FMT = "spe_result_ped{ped}_led{led}_board{board}.root"

# objects (exactly as your screenshot)
HIST_FMT  = "ped{ped}_led{led}_board{board}_hist"
FUNC_FMT  = "ped{ped}_led{led}_board{board}_func"

PAR_Q1  = "Q1"
PAR_PXT = "pxt"

NBOARDS = 4
NCH_PER_BOARD = 16  # for 4x4 grid trend pages


# -------------------------
# Helpers
# -------------------------
def read_info(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            led_run = int(parts[0])     # 1st column = LED run
            ped_run = int(parts[1])     # 2nd column = pedestal run
            ch      = int(parts[2])
            mppc_hv = float(parts[3])
            led_hv  = float(parts[4])
            rows.append({
                "led": led_run,
                "ped": ped_run,
                "ch": ch,
                "mppc_hv": mppc_hv,
                "led_hv": led_hv,
            })
    return rows


def get_tf1_par(tf1, name):
    idx = tf1.GetParNumber(name)
    if idx < 0:
        return None
    return tf1.GetParameter(idx), tf1.GetParError(idx)


def text_box(lines, x1=0.12, y1=0.78, x2=0.88, y2=0.93, size=0.030):
    p = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    p.SetFillStyle(0)
    p.SetBorderSize(0)
    p.SetTextAlign(12)
    p.SetTextFont(42)
    p.SetTextSize(size)
    for s in lines:
        p.AddText(s)
    return p


def make_graph(points, title, xlab, ylab):
    # points: list of (x, y, ey)
    points = sorted(points, key=lambda t: t[0])
    g = ROOT.TGraphErrors(len(points))
    for i, (x, y, ey) in enumerate(points):
        g.SetPoint(i, x, y)
        g.SetPointError(i, 0.0, ey if ey is not None else 0.0)
    g.SetTitle(f"{title};{xlab};{ylab}")
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.9)
    return g


def draw_board_grid(canvas, board, data_dict, ylab, out_pdf):
    """
    data_dict key: (board, ch) -> list of (led_hv, val, err)
    """
    canvas.Clear()
    canvas.Divide(4, 4, 0.001, 0.001)
    for ch in range(NCH_PER_BOARD):
        canvas.cd(ch + 1)
        ROOT.gPad.SetGrid(1, 1)
        ROOT.gPad.SetMargin(0.14, 0.05, 0.14, 0.08)

        pts = data_dict.get((board, ch), [])
        if len(pts) == 0:
            ROOT.TLatex().DrawLatexNDC(0.15, 0.55, f"b{board} ch{ch}")
            ROOT.TLatex().DrawLatexNDC(0.15, 0.45, "no data")
            continue

        g = make_graph(pts, f"b{board} ch{ch}", "LED HV", ylab)
        g.Draw("AP")
        text_box([f"b{board} ch{ch}"], x1=0.15, y1=0.82, x2=0.95, y2=0.95, size=0.05).Draw()

    canvas.Print(out_pdf)


# -------------------------
# Main
# -------------------------
def main():
    rows = read_info(INFO_FILE)
    if len(rows) == 0:
        raise RuntimeError("No valid rows in led_run_info.txt")

    # trend: (board,ch) -> list of (led_hv, Q1, Q1err), etc.
    trend_q1  = defaultdict(list)
    trend_pxt = defaultdict(list)

    c = ROOT.TCanvas("c", "c", 1100, 800)
    c.Print(OUT_PDF + "[")  # open multipage

    # -----------------------------------------
    # 1) Per-file pages: loop over INFO first
    # -----------------------------------------
    for r in rows:
        led = r["led"]
        ped = r["ped"]
        ch  = r["ch"]
        mppc_hv = r["mppc_hv"]
        led_hv  = r["led_hv"]

        for board in range(NBOARDS):
            rootfile = os.path.join(ROOT_DIR, FNAME_FMT.format(ped=ped, led=led, board=board))
            if not os.path.exists(rootfile):
                # info-driven: if missing file, just note and continue
                print(f"[miss] {os.path.basename(rootfile)} (no ROOT file)")
                continue

            base_hist = HIST_FMT.format(ped=ped, led=led, board=board)
            base_func = FUNC_FMT.format(ped=ped, led=led, board=board)

            tf = ROOT.TFile.Open(rootfile, "READ")
            if not tf or tf.IsZombie():
                print(f"[skip] cannot open: {rootfile}")
                continue

            h  = tf.Get(base_hist)
            f1 = tf.Get(base_func)

            # update trend (keyed by actual channel from info)
            if f1:
                q1  = get_tf1_par(f1, PAR_Q1)
                pxt = get_tf1_par(f1, PAR_PXT)
                key = (board, ch)
                if q1:
                    trend_q1[key].append((led_hv, float(q1[0]), float(q1[1])))
                if pxt:
                    trend_pxt[key].append((led_hv, float(pxt[0]), float(pxt[1])))

            # draw per-file page (2 pads)
            c.Clear()
            c.Divide(1, 2)

            c.cd(1)
            ROOT.gPad.SetGrid(1, 1)
            if h:
                h.SetLineWidth(2)
                h.Draw("HIST")
            else:
                ROOT.TLatex().DrawLatexNDC(0.2, 0.6, f"Missing histogram: {base_hist}")

            c.cd(2)
            ROOT.gPad.SetGrid(1, 1)
            if h:
                h.Draw("HIST")
                if f1:
                    f1.SetLineWidth(2)
                    f1.Draw("SAME")
            else:
                ROOT.TLatex().DrawLatexNDC(0.2, 0.6, "No overlay available")

            # annotation
            lines = [
                f"file: {os.path.basename(rootfile)}",
                f"led_run={led}, ped_run={ped}, board={board}, ch={ch}",
                f"MPPC_HV={mppc_hv}, LED_HV={led_hv}",
            ]
            if f1:
                q1  = get_tf1_par(f1, PAR_Q1)
                pxt = get_tf1_par(f1, PAR_PXT)
                if q1:  lines.append(f"Q1  = {q1[0]:.6g} ± {q1[1]:.2g}")
                if pxt: lines.append(f"pxt = {pxt[0]:.6g} ± {pxt[1]:.2g}")

            c.cd(1)
            text_box(lines).Draw()

            c.Print(OUT_PDF)
            tf.Close()

    # -----------------------------------------
    # 2) Trend pages: board-wise 4x4 grid
    # -----------------------------------------
    for board in range(NBOARDS):
        c.Clear()
        ROOT.TLatex().SetTextSize(0.05)
        ROOT.TLatex().DrawLatexNDC(0.12, 0.55, f"Trend summary: Q1 vs LED HV (board {board})")
        c.Print(OUT_PDF)
        draw_board_grid(c, board, trend_q1, "Q1", OUT_PDF)

    for board in range(NBOARDS):
        c.Clear()
        ROOT.TLatex().SetTextSize(0.05)
        ROOT.TLatex().DrawLatexNDC(0.12, 0.55, f"Trend summary: pxt vs LED HV (board {board})")
        c.Print(OUT_PDF)
        draw_board_grid(c, board, trend_pxt, "pxt", OUT_PDF)

    c.Print(OUT_PDF + "]")  # close
    print(f"[done] wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
