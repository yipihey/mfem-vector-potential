#!/usr/bin/env python3
"""E5 plots.  Usage: scripts/plot_exp5.py [csv] [outprefix]
Reads results/exp5_ale.csv (main/lag series, b0 = 0 unless noted) and writes
  results/exp5_errB.png   rel. B error vs t   (rows: p/eps, columns: #rezones; lines per operator)
  results/exp5_energy.png E_B/E_exact vs t
  results/exp5_helicity.png  H/H(0) vs t  (A-route)
  results/exp5_div.png    scaled div vs t
  results/exp5_detj.png   min det J vs t
  results/exp5_b0.png     the B0-on series (p=2, eps 0.3, 5 rezones)
At a rezone step two points share the same t: the state before (on the deformed mesh) and after
(on the new uniform mesh) the transfer, so every rezone shows up as a vertical jump."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vp_csv import load

csv_path = sys.argv[1] if len(sys.argv) > 1 else "results/exp5_ale.csv"
prefix = sys.argv[2] if len(sys.argv) > 2 else "results/exp5"
rows = [r for r in load(csv_path) if r["kind"] == "step"]

# (op, gauge) variants
VARS = [("A_int", "none"), ("A_int", "coulomb"), ("A_l2", "coulomb"), ("B_int", "none"), ("B_l2", "none"), ("B_l2c", "none")]
LAB = {("A_int", "none"): "A_int", ("A_int", "coulomb"): "A_int + Coulomb", ("A_l2", "coulomb"): "A_l2 + Coulomb",
       ("A_l2", "none"): "A_l2", ("A_pt", "none"): "A_pt",
       ("B_int", "none"): "B_int", ("B_l2", "none"): "B_l2", ("B_l2c", "none"): "B_l2c"}
COL = {("A_int", "none"): "#E69F00", ("A_int", "coulomb"): "#D55E00", ("A_l2", "coulomb"): "#CC79A7", ("A_l2", "none"): "#999999",
       ("A_pt", "none"): "#7a4a00", ("B_int", "none"): "#56B4E9", ("B_l2", "none"): "#009E73", ("B_l2c", "none"): "#000000"}
LS = {v: ("-" if v[0][0] == "A" else "--") for v in COL}


def sel(**kw):
    out = [r for r in rows if all(r.get(k) == v for k, v in kw.items())]
    return sorted(out, key=lambda r: (r["step"], r["rezone"]))


def lag(p, eps, q=2):
    return sel(p=float(p), eps=eps, nrez=0.0, q=float(q), op="A_int", b0=0.0)


CASES = [(2, 0.3), (2, 0.6), (3, 0.3), (3, 0.6)]
NREZ = [1, 5, 20]


def grid(metric, ylabel, fname, logy=True, ops=VARS, lagref=True, title=None, only_A=False, ylim=None):
    fig, axs = plt.subplots(len(CASES), len(NREZ), figsize=(4.4 * len(NREZ), 2.7 * len(CASES)), squeeze=False, sharex=True)
    for i, (p, eps) in enumerate(CASES):
        for j, nr in enumerate(NREZ):
            ax = axs[i, j]
            if lagref:
                L = lag(p, eps)
                if L:
                    ax.plot([r["t"] for r in L], [r[metric] for r in L], color="#888888", lw=2.0, alpha=0.6, label="Lagrangian (K=0)")
            for v in ops:
                if only_A and v[0][0] != "A":
                    continue
                S = sel(p=float(p), eps=eps, nrez=float(nr), op=v[0], gauge=v[1], b0=0.0, tag="main")
                if not S:
                    continue
                y = [r[metric] for r in S]
                ax.plot([r["t"] for r in S], y, ls=LS[v], color=COL[v], lw=1.3, label=LAB[v])
            if logy:
                ax.set_yscale("log")
            if ylim:
                ax.set_ylim(*ylim)
            ax.grid(alpha=0.25)
            ax.set_title(f"p={p}, eps={eps}, {nr} rezone{'s' if nr > 1 else ''}", fontsize=9)
            if j == 0:
                ax.set_ylabel(ylabel, fontsize=9)
            if i == len(CASES) - 1:
                ax.set_xlabel("t")
    h, l = axs[0, 0].get_legend_handles_labels()
    hh = dict(zip(l, h))
    for r in range(len(CASES)):
        for c in range(len(NREZ)):
            h2, l2 = axs[r, c].get_legend_handles_labels()
            hh.update(dict(zip(l2, h2)))
    fig.legend(hh.values(), hh.keys(), loc="upper center", ncol=7, fontsize=8, frameon=False)
    if title:
        fig.suptitle(title, y=0.995, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(fname, dpi=110)
    plt.close(fig)
    print("wrote", fname)


for r in rows:
    r["E_ratio"] = r["E_over_Eex"]
grid("errB_rel", r"$\|B_h-B_{pf}\|/\|B_{pf}\|$", prefix + "_errB.png")
grid("E_over_Eex", r"$E_B(t)/E_{exact}(t)$", prefix + "_energy.png", logy=False, ylim=(0.9, 1.1))
grid("H_over_H0", r"$H(t)/H(0)$ (A-route)", prefix + "_helicity.png", logy=False, only_A=True, ylim=(0.9, 1.1))
grid("div_scaled", r"$\|\mathrm{div}\,B_h\|/(\|B\|/h)$", prefix + "_div.png")

# min detJ: Lagrangian runs only + one rezoned example
fig, axs = plt.subplots(1, 2, figsize=(9, 3.2))
for k, eps in enumerate([0.3, 0.6]):
    ax = axs[k]
    L = lag(2, eps)
    ax.plot([r["t"] for r in L], [r["minDetJ"] for r in L], color="#888888", lw=2, label="Lagrangian (K=0)")
    for nr, c in zip([1, 5, 20], ["#0072B2", "#009E73", "#D55E00"]):
        S = sel(p=2.0, eps=eps, nrez=float(nr), op="B_l2c", gauge="none", b0=0.0, tag="main")
        if S:
            ax.plot([r["t"] for r in S], [r["minDetJ"] for r in S], color=c, lw=1.2, label=f"{nr} rezone(s)")
    ax.set_title(f"min detJ (relative to uniform cell), eps={eps}", fontsize=9)
    ax.set_xlabel("t")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(prefix + "_detj.png", dpi=110)
plt.close(fig)
print("wrote", prefix + "_detj.png")

# b0 series
b0rows = [r for r in rows if r["b0"] == 1.0]
if b0rows:
    fig, axs = plt.subplots(1, 4, figsize=(17, 3.4))
    metrics = [("errB_rel", "rel. B error", True), ("E_over_Eex", r"$E_B/E_{exact}$", False),
               ("flux_err", "max slice-flux error", True), ("div_scaled", "scaled div", True)]
    for ax, (m, lab, logy) in zip(axs, metrics):
        L = [r for r in rows if r["b0"] == 1.0 and r["nrez"] == 0.0]
        L = sorted(L, key=lambda r: (r["step"], r["rezone"]))
        if L:
            ax.plot([r["t"] for r in L], [r[m] for r in L], color="#888888", lw=2, label="Lagrangian (A_int, K=0)")
        for v in [("B_int", "none"), ("B_l2", "none"), ("B_l2c", "none")]:
            S = sorted([r for r in b0rows if r["op"] == v[0] and r["nrez"] == 5.0], key=lambda r: (r["step"], r["rezone"]))
            if S:
                ax.plot([r["t"] for r in S], [max(r[m], 1e-17) if logy else r[m] for r in S], ls=LS[v], color=COL[v], label=LAB[v])
        if logy:
            ax.set_yscale("log")
        ax.set_xlabel("t")
        ax.set_title(lab + " (B0 on, p=2, eps=0.3, 5 rezones)", fontsize=8)
        ax.grid(alpha=0.25)
    axs[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(prefix + "_b0.png", dpi=110)
    plt.close(fig)
    print("wrote", prefix + "_b0.png")


# gauge comparison and q=3 control (p=2 eps=0.6 5 rezones; p=3 eps=0.6 5 rezones q=2 vs q=3)
allr = [r for r in load(csv_path) if r["kind"] == "step"]


def selx(**kw):
    return sorted([r for r in allr if all(r.get(k) == v for k, v in kw.items())], key=lambda r: (r["step"], r["rezone"]))


GV = [("A_int", "none", "#E69F00", "-", "A_int"), ("A_int", "coulomb", "#D55E00", "-", "A_int + Coulomb (after)"),
      ("A_int", "coulomb_pre", "#8B2500", ":", "A_int + Coulomb (before, on deformed mesh)"),
      ("A_l2", "coulomb", "#CC79A7", "-", "A_l2 + Coulomb (after)"), ("A_l2", "coulomb_pre", "#7B3F65", ":", "A_l2 + Coulomb (before)"),
      ("B_l2c", "none", "#000000", "--", "B_l2c")]
fig, axs = plt.subplots(2, 3, figsize=(15, 6.2), sharex=True)
for i, (p, eps) in enumerate([(2, 0.6), (3, 0.6)]):
    for j, nr in enumerate([2, 5, 20]):
        ax = axs[i, j]
        for op, g, c, ls, lab in GV:
            S_ = [r for r in selx(p=float(p), eps=eps, nrez=float(nr), op=op, gauge=g, b0=0.0) if r["tag"] in ("main", "gpre")]
            if S_:
                ax.plot([r["t"] for r in S_], [r["errB_rel"] for r in S_], color=c, ls=ls, lw=1.3, label=lab)
        ax.set_yscale("log")
        ax.grid(alpha=0.25)
        ax.set_title(f"p={p}, eps={eps}, {nr} rezones (q=2)", fontsize=9)
        if j == 0:
            ax.set_ylabel(r"$\|B_h-B_{pf}\|/\|B_{pf}\|$")
        if i == 1:
            ax.set_xlabel("t")
h, l = axs[0, 1].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=3, fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.92))
fig.savefig(prefix + "_gauge.png", dpi=110)
plt.close(fig)
print("wrote", prefix + "_gauge.png")

fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
for ax, q in zip(axs, [2, 3]):
    for op, g, c, ls, lab in [("B_l2", "none", "#009E73", "--", "B_l2"), ("B_l2c", "none", "#000000", "--", "B_l2c"),
                              ("A_int", "coulomb", "#D55E00", "-", "A_int + Coulomb"), ("A_l2", "coulomb", "#CC79A7", "-", "A_l2 + Coulomb")]:
        tag = "main" if q == 2 else "q3"
        S_ = selx(p=3.0, eps=0.6, nrez=5.0, op=op, gauge=g, b0=0.0, tag=tag)
        if S_:
            ax.plot([r["t"] for r in S_], [r["errB_rel"] for r in S_], color=c, ls=ls, label=lab)
    L = selx(p=3.0, eps=0.6, nrez=0.0, q=float(q), op="A_int", b0=0.0, tag="lag" if q == 2 else "lag_q3")
    if L:
        ax.plot([r["t"] for r in L], [r["errB_rel"] for r in L], color="#888888", lw=2, alpha=0.7, label="Lagrangian (K=0)")
    ax.set_yscale("log")
    ax.grid(alpha=0.25)
    ax.set_title(f"p=3, eps=0.6, 5 rezones, geometric order q={q}", fontsize=9)
    ax.set_xlabel("t")
axs[0].set_ylabel(r"$\|B_h-B_{pf}\|/\|B_{pf}\|$")
axs[0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(prefix + "_q3.png", dpi=110)
plt.close(fig)
print("wrote", prefix + "_q3.png")
