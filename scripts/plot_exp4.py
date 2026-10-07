#!/usr/bin/env python3
"""E4 plots.  Usage: scripts/plot_exp4.py [csv] [outprefix]
Writes results/exp4_s1.png (p=1,2,3, no B0), exp4_s2.png (B0 on), exp4_long.png (n=400)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vp_csv import load

csv_path = sys.argv[1] if len(sys.argv) > 1 else "results/exp4_repeat.csv"
prefix = sys.argv[2] if len(sys.argv) > 2 else "results/exp4"
rows = load(csv_path)
for _r in rows:   # s3b (A_l2, B_l2 long run) is part of the long-run series
    if _r["tag"] == "s3b":
        _r["tag"] = "s3"
OPS = ["A_pt", "A_int", "A_l2", "B_pt", "B_int", "B_l2", "B_l2c"]
COL = dict(A_pt="#D55E00", A_int="#E69F00", A_l2="#CC79A7", B_pt="#56B4E9", B_int="#0072B2",
           B_l2="#009E73", B_l2c="#000000")
MK = dict(A_pt="o", A_int="s", A_l2="^", B_pt="o", B_int="s", B_l2="^", B_l2c="D")
LS = {o: ("-" if o[0] == "A" else "--") for o in OPS}

METRICS = [
    ("E_over_Eref", r"$E_B(n)/E_B^{ref}$ (same mesh)", True),
    ("errB_rel", r"$\|B_n-B\|/\|B\|$", True),
    ("div_rel_L2", r"$\|D_hb_n\|/(\|B\|/h)$", True),
    ("H_dev_scaled", r"$(H_n-H^{ref})/(\|a\|\|B\|)$ (A-route)", False),
    ("a_over_aref", r"$\|a_n\|/\|a^{ref}\|$ (A-route)", True),
    ("grad_frac", r"gradient fraction of $a_n$ (A-route)", True),
]


def series(op, base, mesh="A"):
    s = [r for r in rows if r["op"] == op and r["mesh"] == mesh and all(r.get(k) == v for k, v in base.items())]
    return sorted(s, key=lambda r: r["step"])


def cutoff(op, base, thr=3.0):
    """first step at which the B error exceeds thr (either mesh): the series is drawn up to there and marked with x"""
    allr = sorted([r for r in rows if r["op"] == op and all(r.get(k) == v for k, v in base.items())], key=lambda r: r["step"])
    for r in allr:
        if r["errB_rel"] > thr:
            return r["step"]
    return None


def grid(base_list, titles, fname, suptitle, ops=OPS):
    fig, axs = plt.subplots(len(METRICS), len(base_list), figsize=(4.6 * len(base_list), 2.6 * len(METRICS)), squeeze=False)
    for j, (base, t) in enumerate(zip(base_list, titles)):
        for i, (k, lab, logy) in enumerate(METRICS):
            ax = axs[i, j]
            for op in ops:
                s = series(op, base)
                if not s:
                    continue
                cut = cutoff(op, base)
                if cut is not None:
                    s = [r for r in s if r["step"] <= cut]
                x = [r["step"] for r in s]
                y = [r[k] for r in s]
                if np.all(np.isnan(y)):
                    continue
                if k == "grad_frac":
                    y = [max(v, 1e-12) for v in y]
                (ax.semilogy if logy else ax.plot)(x, y, LS[op], marker=MK[op], ms=3, color=COL[op], label=op, lw=1.3)
                if cut is not None and x:
                    ax.plot([x[-1]], [y[-1]], "x", color=COL[op], ms=8, mew=2)
            ax.grid(True, which="both", alpha=0.3)
            if i == 0: ax.set_title(t, fontsize=10)
            if j == 0: ax.set_ylabel(lab, fontsize=8)
            if i == len(METRICS) - 1: ax.set_xlabel("remap number n (state on mesh A, even n)")
    axs[0, 0].legend(fontsize=6, ncol=2)
    fig.suptitle(suptitle + "  (x: cut where B err > 300%)", fontsize=9); fig.tight_layout(); fig.savefig(fname, dpi=120); plt.close(fig)
    print("wrote", fname)


grid([dict(tag=f"s1_p{p}") for p in (1, 2, 3)], [f"p={p}, N=8, no B0" for p in (1, 2, 3)],
     prefix + "_s1.png", "E4: ping-pong uniform <-> f1 eps 0.3 (solid A-route, dashed B-route)")
if any(r["tag"] == "s2" for r in rows):
    grid([dict(tag="s2")], ["p=2, N=8, with B0"], prefix + "_s2.png", "E4 with mean field B0")
if any(r["tag"] == "s3" for r in rows):
    grid([dict(tag="s3")], ["p=2, n=400"], prefix + "_long.png", "E4 long run",
         ops=["A_pt", "A_int", "A_l2", "B_int", "B_l2", "B_l2c"])
