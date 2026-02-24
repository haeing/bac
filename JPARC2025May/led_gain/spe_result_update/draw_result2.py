#!/usr/bin/env python3
import os
import re
import glob
import math
import csv

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)  # disable ROOT stats box

# ----------------------------
# Robust parser for led_run_info.txt
# ----------------------------
_num_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

def load_led_run_info(path):
    """
    Parse ONLY lines that contain >=5 numeric tokens:
      led_run ped_run ch mppc_hv led_hv
    Ignore:
      headers like "1#runnumber,..."
      separator lines like "5"
    Return dict keyed by led_run.
    """
    by_led = {}
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            # remove inline comment (also handles "1#header" by cutting to "1")
            if "#" in line:
                line = line.split("#", 1)[0].strip()
                if not line:
                    continue

            nums = _num_re.findall(line)
            if len(nums) < 5:
                continue

            led_run = int(float(nums[0]))
            ped_run = int(float(nums[1]))
            ch      = int(float(nums[2]))
            mppc_hv = float(nums[3])
            led_hv  = float(nums[4])

            by_led[led_run] = {
                "led_run": led_run,
                "ped_run": ped_run,
                "mppc_ch": ch,
                "mppc_hv": mppc_hv,
                "led_hv": led_hv,
            }
    return by_led

def parse_from_filename(fname):
    """
    spe_result_ped437_led434_board3.root -> (437, 434, 3)
    """
    m = re.search(r"ped(\d+)_led(\d+)_board(\d+)", os.path.basename(fname))
    if not m:
        return None, None, None
    return int(m.group(1)), int(m.group(2)), int(m.group(3))

# ----------------------------
# ROOT helpers
# ----------------------------
def get_root_objects(tf, ped, led, board):
    """
    Keys expected:
      TFitResult : ped{ped}_led{led}_board{board}
      TH1D      : ped{ped}_led{led}_board{board}_hist
      TF1       : ped{ped}_led{led}_board{board}_func
    """
    tag = f"ped{ped}_led{led}_board{board}"
    h    = tf.Get(f"{tag}_hist")
    f1   = tf.Get(f"{tag}_func")
    fitr = tf.Get(tag)
    return tag, h, f1, fitr

def get_fit_quality(fitr):
    if not fitr:
        return float("nan"), -1
    try:
        return float(fitr.Chi2()), int(fitr.Ndf())
    except Exception:
        return float("nan"), -1

def get_tf1_par_by_name(f1, names):
    """
    Search TF1 parameter by name; returns (value, error) or None.
    """
    if not f1:
        return None
    want = {n.lower() for n in names}
    for i in range(f1.GetNpar()):
        pn = f1.GetParName(i)
        if pn and pn.lower() in want:
            return float(f1.GetParameter(i)), float(f1.GetParError(i))
    return None

def safe_fmt(v, fmt="{:.4g}"):
    if v is None:
        return "N/A"
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return "N/A"
    try:
        return fmt.format(v)
    except Exception:
        return str(v)

def draw_info_box(lines, x1=0.12, y1=0.72, x2=0.62, y2=0.93, text_size=0.04):
    """
    Top-left info box in NDC.
    IMPORTANT: return the box and keep reference to avoid PyROOT GC issue.
    """
    box = ROOT.TPaveText(x1, y1, x2, y2, "NDC")
    box.SetFillStyle(0)   # transparent
    box.SetBorderSize(0)
    box.SetTextAlign(12)  # left
    box.SetTextFont(42)
    box.SetTextSize(text_size)
    for s in lines:
        box.AddText(s)
    box.Draw()
    ROOT.SetOwnership(box, False)
    return box

# ----------------------------
# CSV writer
# ----------------------------
def write_csv(results, out_csv):
    """
    Writes CSV with explicit MPPC fields.
    """
    fields = [
        "rootfile",
        "mppc_board",
        "mppc_ch",
        "led_run",
        "ped_run",
        "mppc_hv",
        "led_hv",
        "chi2",
        "ndf",
        "q1",
        "q1_err",
        "pxt",
        "pxt_err",
    ]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in fields})

# ----------------------------
# Auto-detect (no args)
# ----------------------------
def autodetect():
    """
    If running inside spe_result_update/:
      info: ../led_run_info.txt (fallback: ./led_run_info.txt)
      indir: current dir
    """
    cwd = os.getcwd()
    info1 = os.path.join(cwd, "led_run_info.txt")
    info2 = os.path.join(os.path.dirname(cwd), "led_run_info.txt")
    info = info1 if os.path.isfile(info1) else info2

    indir = cwd
    pattern = "spe_result_ped*_led*_board*.root"
    outpdf = os.path.join(indir, "spe_summary.pdf")
    outcsv = os.path.join(indir, "spe_summary.csv")
    return info, indir, pattern, outpdf, outcsv

# ----------------------------
# Main
# ----------------------------
def main():
    info, indir, pattern, outpdf, outcsv = autodetect()

    if not os.path.isfile(info):
        print("[ERR] led_run_info.txt not found.")
        print("Tried:")
        print("  ./led_run_info.txt")
        print("  ../led_run_info.txt")
        return

    by_led = load_led_run_info(info)

    rootfiles = sorted(glob.glob(os.path.join(indir, pattern)))
    if not rootfiles:
        print(f"[ERR] No ROOT files matched: {os.path.join(indir, pattern)}")
        return

    results = []
    keepalive = []  # keep TPaveText alive across printing (PyROOT GC workaround)

    c = ROOT.TCanvas("c", "", 900, 900)
    c.Print(outpdf + "[")

    for rootfile in rootfiles:
        ped_name, led_name, board = parse_from_filename(rootfile)
        if led_name is None or board is None:
            print(f"[WARN] cannot parse filename: {rootfile}")
            continue

        info_row = by_led.get(led_name, None)

        ped_run = info_row["ped_run"] if info_row else ped_name
        mppc_ch = info_row["mppc_ch"] if info_row else -1
        mppc_hv = info_row["mppc_hv"] if info_row else float("nan")
        led_hv  = info_row["led_hv"]  if info_row else float("nan")

        tf = ROOT.TFile.Open(rootfile, "READ")
        if (not tf) or tf.IsZombie():
            print(f"[WARN] cannot open: {rootfile}")
            continue

        tag, h, f1, fitr = get_root_objects(tf, ped_run, led_name, board)
        chi2, ndf = get_fit_quality(fitr)

        q1  = get_tf1_par_by_name(f1, ["Q1", "q1"])
        pxt = get_tf1_par_by_name(f1, ["pxt", "PXT", "p_xt", "Pxt", "qxt", "Qxt", "QXT"])

        results.append({
            "rootfile": os.path.basename(rootfile),
            "mppc_board": board,
            "mppc_ch": mppc_ch,
            "led_run": led_name,
            "ped_run": ped_run,
            "mppc_hv": mppc_hv,
            "led_hv": led_hv,
            "chi2": chi2,
            "ndf": ndf,
            "q1": (q1[0] if q1 else ""),
            "q1_err": (q1[1] if q1 else ""),
            "pxt": (pxt[0] if pxt else ""),
            "pxt_err": (pxt[1] if pxt else ""),
        })

        # -------- Draw page (2 pads) --------
        c.Clear()
        c.Divide(1, 2)

        q1_str  = f"{q1[0]:.6g}" if q1 else "N/A"
        pxt_str = f"{pxt[0]:.6g}" if pxt else "N/A"

        print(mppc_ch)
        label = f"MPPC board={int(board)}, ch={int(mppc_ch)}"
        lines = [
            str(label),
            f"MPPC_HV={safe_fmt(mppc_hv,'{:.4g}')} V, LED_HV={safe_fmt(led_hv,'{:.4g}')} V",
            f"#chi^{{2}}={safe_fmt(chi2,'{:.4g}')}, ndf={ndf}",
            f"Q1 = {q1_str}",
            f"pxt(Qxt) = {pxt_str}",
        ]

        # pad1
        c.cd(1)
        ROOT.gPad.SetGrid(1, 1)
        if h:
            h.SetLineWidth(2)
            h.SetTitle("")
            h.Draw("HIST")
        else:
            ROOT.TLatex().DrawLatexNDC(0.2, 0.6, f"Missing histogram: {tag}_hist")
        b1 = draw_info_box(lines)
        keepalive.append(b1)
        ROOT.gPad.Update()

        # pad2
        c.cd(2)
        ROOT.gPad.SetGrid(1, 1)
        if h:
            h.Draw("HIST")
            if f1:
                f1.SetLineWidth(2)
                f1.Draw("SAME")
        else:
            ROOT.TLatex().DrawLatexNDC(0.2, 0.6, "No overlay available")
        b2 = draw_info_box(lines)
        keepalive.append(b2)
        ROOT.gPad.Update()

        c.Print(outpdf)
        tf.Close()

    c.Print(outpdf + "]")

    # -------- CSV --------
    write_csv(results, outcsv)

    print(f"[OK] PDF : {outpdf}")
    print(f"[OK] CSV : {outcsv}")
    print(f"[OK] Nfiles: {len(results)}")

if __name__ == "__main__":
    main()
