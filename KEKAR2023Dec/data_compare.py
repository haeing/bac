#!/usr/bin/env python3

import ROOT
import numpy as np
from collections import defaultdict

ROOT_FILE = "kek_result.root"

HV_TARGET = 58

CONFIGS = {
    2: {"thre": 60},
    3: {"thre": 75},
}

xs = np.array([-48, -32, -16, 0, 16, 32, 48], dtype=float)
ys = np.array([-54, -36, -18, 0, 18, 36, 54], dtype=float)


def find_tree(f):
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("TTree"):
            return obj
    raise RuntimeError("No TTree found.")


def find_branch(tree, candidates):
    branches = [b.GetName() for b in tree.GetListOfBranches()]

    for name in candidates:
        if name in branches:
            return name

    print("")
    print("Available branches:")
    for b in branches:
        print("  ", b)

    raise RuntimeError(f"Could not find branch from candidates: {candidates}")


def nearest_scan_position(x, y):
    ix = int(np.argmin(np.abs(xs - x)))
    iy = int(np.argmin(np.abs(ys - y)))

    return xs[ix], ys[iy]


def read_data(tree, thick, thre, br):
    values = []

    for ev in tree:
        hv_val = int(getattr(ev, br["hv"]))
        thick_val = int(getattr(ev, br["thick"]))
        thre_val = int(getattr(ev, br["thre"]))

        if hv_val != HV_TARGET:
            continue
        if thick_val != thick:
            continue
        if thre_val != thre:
            continue

        x = float(getattr(ev, br["x"]))
        y = float(getattr(ev, br["y"]))
        npe = float(getattr(ev, br["npe"]))

        if not np.isfinite(npe):
            continue

        x_scan, y_scan = nearest_scan_position(x, y)

        values.append({
            "x": x,
            "y": y,
            "x_scan": x_scan,
            "y_scan": y_scan,
            "npe": npe,
        })

    return values


def make_position_average(values):
    grouped = defaultdict(list)

    for v in values:
        key = (v["x_scan"], v["y_scan"])
        grouped[key].append(v["npe"])

    pos_mean = {}

    for key, arr in grouped.items():
        arr = np.array(arr, dtype=float)
        pos_mean[key] = {
            "mean": np.mean(arr),
            "std": np.std(arr),
            "n": len(arr),
        }

    return pos_mean


def evaluate_array(arr, label):
    arr = np.array(arr, dtype=float)
    arr = arr[np.isfinite(arr)]

    vmax = np.max(arr)
    vmin = np.min(arr)
    vmean = np.mean(arr)
    vstd = np.std(arr)

    dev_from_max = (vmax - vmin) / vmax * 100.0
    rms_percent = vstd / vmean * 100.0

    print("")
    print("===================================================")
    print(label)
    print(f"N points          = {len(arr)}")
    print(f"max NPE           = {vmax:.3f}")
    print(f"min NPE           = {vmin:.3f}")
    print(f"mean NPE          = {vmean:.3f}")
    print(f"std NPE           = {vstd:.3f}")
    print(f"(max-min)/max     = {dev_from_max:.2f} %")
    print(f"std/mean          = {rms_percent:.2f} %")
    print("===================================================")

    return {
        "n": len(arr),
        "max": vmax,
        "min": vmin,
        "mean": vmean,
        "std": vstd,
        "dev_percent": dev_from_max,
        "rms_percent": rms_percent,
    }


def print_position_table(pos_mean, label):
    print("")
    print("---------------------------------------------------")
    print(label)
    print("position-averaged NPE")
    print("x_scan  y_scan   mean_NPE   std_NPE   n")
    print("---------------------------------------------------")

    for y in ys:
        for x in xs:
            key = (x, y)
            if key not in pos_mean:
                continue

            v = pos_mean[key]
            print(f"{x:6.1f} {y:7.1f} {v['mean']:10.3f} {v['std']:9.3f} {v['n']:4d}")


def main():
    f = ROOT.TFile.Open(ROOT_FILE)

    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {ROOT_FILE}")

    tree = find_tree(f)

    br = {}
    br["hv"] = find_branch(tree, ["bac_HV"])
    br["thick"] = find_branch(tree, ["bac_thick"])
    br["thre"] = find_branch(tree, ["bac_thre"])
    br["x"] = find_branch(tree, ["x_pos"])
    br["y"] = find_branch(tree, ["y_pos"])
    br["npe"] = find_branch(tree, ["bac_npe"])

    print("")
    print("Using branches")
    for k, v in br.items():
        print(f"  {k:8s}: {v}")

    results_raw = {}
    results_pos = {}
    pos_maps = {}

    for thick, cfg in CONFIGS.items():
        thre = cfg["thre"]

        values = read_data(tree, thick, thre, br)

        print("")
        print(f"[NPE] thick={thick}, thre={thre}, selected entries = {len(values)}")

        if len(values) == 0:
            raise RuntimeError(f"No data selected for thick={thick}, thre={thre}")

        raw_npe = [v["npe"] for v in values]

        pos_mean = make_position_average(values)
        pos_npe = [v["mean"] for v in pos_mean.values()]

        results_raw[thick] = evaluate_array(
            raw_npe,
            label=f"[thick={thick}] raw selected entries"
        )

        results_pos[thick] = evaluate_array(
            pos_npe,
            label=f"[thick={thick}] discrete scan-position averages"
        )

        print_position_table(
            pos_mean,
            label=f"[thick={thick}]"
        )

        pos_maps[thick] = pos_mean

    print("")
    print("###################################################")
    print("[THICKNESS COMPARISON: raw selected entries]")
    print("###################################################")
    print(f"mean ratio 3/2 = {results_raw[3]['mean'] / results_raw[2]['mean']:.3f}")
    print(f"max  ratio 3/2 = {results_raw[3]['max']  / results_raw[2]['max'] :.3f}")
    print(f"min  ratio 3/2 = {results_raw[3]['min']  / results_raw[2]['min'] :.3f}")
    print(f"expected thickness ratio = {3.0 / 2.0:.3f}")

    print("")
    print("###################################################")
    print("[THICKNESS COMPARISON: scan-position averages]")
    print("###################################################")
    print(f"mean ratio 3/2 = {results_pos[3]['mean'] / results_pos[2]['mean']:.3f}")
    print(f"max  ratio 3/2 = {results_pos[3]['max']  / results_pos[2]['max'] :.3f}")
    print(f"min  ratio 3/2 = {results_pos[3]['min']  / results_pos[2]['min'] :.3f}")
    print(f"expected thickness ratio = {3.0 / 2.0:.3f}")

    common_keys = sorted(set(pos_maps[2].keys()) & set(pos_maps[3].keys()))

    ratios = []

    for key in common_keys:
        npe2 = pos_maps[2][key]["mean"]
        npe3 = pos_maps[3][key]["mean"]

        if npe2 > 0:
            ratios.append(npe3 / npe2)

    ratios = np.array(ratios, dtype=float)

    print("")
    print("###################################################")
    print("[POINT-BY-POINT THICKNESS RATIO: common positions]")
    print("###################################################")
    print(f"N common positions = {len(ratios)}")
    print(f"mean ratio 3/2     = {np.mean(ratios):.3f}")
    print(f"std ratio          = {np.std(ratios):.3f}")
    print(f"min ratio          = {np.min(ratios):.3f}")
    print(f"max ratio          = {np.max(ratios):.3f}")
    print(f"expected ratio     = {3.0 / 2.0:.3f}")

    print("")
    print("x_scan  y_scan   NPE_thick2   NPE_thick3   ratio_3_over_2")
    print("----------------------------------------------------------")

    for key in common_keys:
        npe2 = pos_maps[2][key]["mean"]
        npe3 = pos_maps[3][key]["mean"]
        ratio = npe3 / npe2 if npe2 > 0 else np.nan
        print(f"{key[0]:6.1f} {key[1]:7.1f} {npe2:11.3f} {npe3:11.3f} {ratio:14.3f}")

    print("")
    print("DONE.")


if __name__ == "__main__":
    main()
