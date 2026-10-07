#!/usr/bin/env python3
"""Plots for the gauge-fixed A-route (E4-gauge).
Usage: scripts/plot_exp4_gauge.py [new_csv] [old_csv] [outprefix]
Writes  <prefix>_p2.png, <prefix>_p3.png (E_B/E_ref, B error, grad fraction, helicity vs n for A_pt/A_int/A_l2
with Coulomb gauge (solid) vs none (dashed, old CSV) vs B_l2c (black)),  <prefix>_long.png (n=400, p=2),
<prefix>_misc.png (gauge iterations/time, jacobi variant, B0 flux)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vp_csv import load

new_path = sys.argv[1] if len(sys.argv) > 1 else "results/exp4_gauge.csv"
old_path = sys.argv[2] if len(sys.argv) > 2 else "results/exp4_repeat.csv"
prefix = sys.argv[3] if len(sys.argv) > 3 else "results/exp4_gauge"
new, old = load(new_path), load(old_path)
for r in old:
    if r["tag"] == "s3b":
        r["tag"] = "s3"
COL = dict(A_pt="#D55E00", A_int="#E69F00", A_l2="#CC79A7", B_l2c="#000000")
MK = dict(A_pt="o", A_int="s", A_l2="^", B_l2c="D")


def sel(rows, op, tag, mesh="A", gauge=None):
    s = [r for r in rows if r["op"] == op and r["tag"] == tag and r["mesh"] == mesh
         and (gauge is None or r.get("gauge") == gauge)]
    return sorted(s, key=lambda r: r["step"])


def cut(s, thr=3.0):
    out = []
    for r in s:
        if r["errB_rel"] > thr:
            break
        out.append(r)
    return out


METRICS = [("E_over_Eref", r"$E_B(n)/E_B^{ref}$", True),
           ("errB_rel", r"$\|B_n-B\|/\|B\|$", True),
           ("grad_frac", r"gradient fraction of $a_n$", True),
           ("H_over_H0", r"$H(n)/H(0)$", False)]


def panel(ax, key, p, tag_new, tag_old, ops=("A_pt", "A_int", "A_l2")):
    for op in ops:
        s = cut(sel(new, op, tag_new))
        if s:
            ax.plot([r["step"] for r in s], [r[key] for r in s], "-", color=COL[op], marker=MK[op], ms=3,
                    label=f"{op} + Coulomb")
        s = cut(sel(old, op, tag_old))
        if s:
            ax.plot([r["step"] for r in s], [r[key] for r in s], "--", color=COL[op], alpha=0.55, label=f"{op} no gauge")
    if key != "H_over_H0" and key != "grad_frac":
        s = sel(old, "B_l2c", tag_old)
        if s:
            ax.plot([r["step"] for r in s], [r[key] for r in s], "-.", color="k", marker="D", ms=3, label="B_l2c")


for p, tn, to in ((2, "g_p2", "s1_p2"), (3, "g_p3", "s1_p3")):
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    for ax, (k, lab, logy) in zip(axs.ravel(), METRICS):
        panel(ax, k, p, tn, to)
        ax.set_xlabel("remaps n (state on mesh A)")
        ax.set_ylabel(lab)
        if logy:
            ax.set_yscale("log")
        ax.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=7)
    fig.suptitle(f"p={p}, N=8, abc, uniform <-> f1 0.3: Coulomb gauge fix after every A-route transfer (series cut at B err>3)")
    fig.tight_layout()
    fig.savefig(f"{prefix}_p{p}.png", dpi=120)
    plt.close(fig)

# long run
fig, axs = plt.subplots(2, 2, figsize=(10, 7))
for ax, (k, lab, logy) in zip(axs.ravel(), METRICS):
    for op in ("A_int", "A_l2"):
        s = sel(new, op, "g_long")
        if s:
            ax.plot([r["step"] for r in s], [r[k] for r in s], "-", color=COL[op], marker=MK[op], ms=3, label=f"{op} + Coulomb")
        s = cut(sel(old, op, "s3"), 1e3)
        if s:
            ax.plot([r["step"] for r in s], [r[k] for r in s], "--", color=COL[op], alpha=0.55, label=f"{op} no gauge")
    if k in ("E_over_Eref", "errB_rel"):
        s = sel(old, "B_l2c", "s3")
        ax.plot([r["step"] for r in s], [r[k] for r in s], "-.", color="k", marker="D", ms=3, label="B_l2c")
    ax.set_xlabel("remaps n"); ax.set_ylabel(lab); ax.grid(alpha=0.3)
    if logy:
        ax.set_yscale("log")
axs[0, 0].legend(fontsize=7)
fig.suptitle("p=2, n=400")
fig.tight_layout()
fig.savefig(f"{prefix}_long.png", dpi=120)
plt.close(fig)

# misc: gauge iterations, jacobi variants, B0 flux
fig, axs = plt.subplots(1, 4, figsize=(17, 4))
for tn, ls in (("g_p2", "-"), ("g_p3", ":"), ("g_long", "--")):
    for op in ("A_pt", "A_int", "A_l2"):
        s = [r for r in new if r["op"] == op and r["tag"] == tn and r["step"] > 0]
        s = sorted(s, key=lambda r: r["step"])
        if s:
            axs[0].plot([r["step"] for r in s], [r["gauge_iters_last"] for r in s], ls, color=COL[op], label=f"{op} {tn}")
axs[0].set_title("CG iterations of the gauge solve"); axs[0].set_xlabel("step"); axs[0].legend(fontsize=6)
jac = [r for r in new if r["tag"] == "g_jac" and r["mesh"] == "A"]
for g, c in (("jacobi:5", "C0"), ("jacobi:20", "C1")):
    s = sorted([r for r in jac if r["gauge"] == g], key=lambda r: r["step"])
    if s:
        axs[1].plot([r["step"] for r in s], [r["grad_frac"] for r in s], "-", color=c, label=f"A_int {g}")
s = sel(new, "A_int", "g_p2")
axs[1].plot([r["step"] for r in s], [r["grad_frac"] for r in s], "-", color="k", label="A_int coulomb")
s = sel(old, "A_int", "s1_p2")
axs[1].plot([r["step"] for r in s], [r["grad_frac"] for r in s], "--", color="gray", label="A_int none")
axs[1].set_yscale("log"); axs[1].set_title("gradient fraction, local (Chebyshev-Jacobi) gauge"); axs[1].legend(fontsize=7)
axs[2].set_title("errB: jacobi variants")
for g, c in (("jacobi:5", "C0"), ("jacobi:20", "C1")):
    s = sorted([r for r in jac if r["gauge"] == g], key=lambda r: r["step"])
    if s:
        axs[2].plot([r["step"] for r in s], [r["errB_rel"] for r in s], "-", color=c, label=g)
s = sel(new, "A_int", "g_p2"); axs[2].plot([r["step"] for r in s], [r["errB_rel"] for r in s], "-", color="k", label="coulomb")
s = sel(old, "A_int", "s1_p2"); axs[2].plot([r["step"] for r in s], [r["errB_rel"] for r in s], "--", color="gray", label="none")
axs[2].set_yscale("log"); axs[2].legend(fontsize=7)
s = sel(new, "A_int", "g_b0")
s0 = sel(old, "A_int", "s2")
axs[3].set_title("B0 on: net flux error, A_int p=2")
axs[3].plot([r["step"] for r in s], [max(r["flux_err"], 1e-17) for r in s], "-", label="A_int + Coulomb")
axs[3].plot([r["step"] for r in s0], [max(r["flux_err"], 1e-17) for r in s0], "--", label="A_int none")
s1 = sel(old, "B_l2c", "s2")
axs[3].plot([r["step"] for r in s1], [max(r["flux_err"], 1e-17) for r in s1], "-.", color="k", label="B_l2c")
axs[3].set_yscale("log"); axs[3].legend(fontsize=7)
for a in axs:
    a.grid(alpha=0.3); a.set_xlabel("remaps n")
fig.tight_layout()
fig.savefig(f"{prefix}_misc.png", dpi=120)
plt.close(fig)
print("wrote", prefix + "_p2.png", prefix + "_p3.png", prefix + "_long.png", prefix + "_misc.png")
