#!/usr/bin/env python3
"""E3 plots.  Usage: scripts/plot_exp3.py [csv] [outprefix]
Writes results/exp3_p_scaling.png, exp3_eps_scaling.png, exp3_gauge.png, exp3_nq.png, exp3_cases.png."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vp_csv import load

csv_path = sys.argv[1] if len(sys.argv) > 1 else "results/exp3_remap.csv"
prefix = sys.argv[2] if len(sys.argv) > 2 else "results/exp3"
rows = load(csv_path)
OPS = ["A_pt", "A_int", "A_l2", "B_pt", "B_int", "B_l2", "B_l2c"]
COL = dict(A_pt="#D55E00", A_int="#E69F00", A_l2="#CC79A7", B_pt="#56B4E9", B_int="#0072B2",
           B_l2="#009E73", B_l2c="#000000")
MK = dict(A_pt="o", A_int="s", A_l2="^", B_pt="o", B_int="s", B_l2="^", B_l2c="D")
LS = {o: ("-" if o[0] == "A" else "--") for o in OPS}


def sel(**kw):
    out = []
    for r in rows:
        if all((r.get(k) == v) for k, v in kw.items()):
            out.append(r)
    return out


def lines(ax, xs_key, y_fn, base, logy=True, ops=OPS):
    for op in ops:
        s = sorted(sel(op=op, **base), key=lambda r: r[xs_key])
        if not s:
            continue
        x = [r[xs_key] for r in s]
        y = [y_fn(r) for r in s]
        (ax.semilogy if logy else ax.plot)(x, y, LS[op], marker=MK[op], ms=4.5, color=COL[op], label=op, lw=1.5)
    ax.grid(True, which="both", alpha=0.3)


METRICS = [
    ("m2_errB_rel", r"$\|B_2-B\|_{L^2}/\|B\|$", lambda r: r["m2_errB_rel"]),
    ("m2_div_rel_L2", r"$\|D_hb_2\|/(\|B\|/h)$  (divergence)", lambda r: max(r["m2_div_rel_L2"], 1e-17)),
    ("dE", r"$|E_2-E_1|/E_1$", lambda r: max(abs(r["m2_dE_rel1"]), 1e-12)),
    ("t", "one-shot time per remap [s]  (setup + apply)", lambda r: r["t_oneshot"]),
]


def metric_fig(xkey, xlabel, base, fname, title, ref_key=True):
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.2))
    for ax, (k, lab, fn) in zip(axs.ravel(), METRICS):
        lines(ax, xkey, fn, base)
        if k == "m2_errB_rel" and ref_key:   # floor: projection of the exact field on M2
            s = sorted(sel(op="A_int", **base), key=lambda r: r[xkey])
            if s:
                ax.semilogy([r[xkey] for r in s], [r["ref_A_errB_rel"] for r in s], ":", color="gray", lw=2,
                            label="floor: projection of exact B on M2")
        ax.set_xlabel(xlabel); ax.set_ylabel(lab)
        if xkey == "p":
            ax.set_xticks([1, 2, 3, 4])
    axs[0, 0].legend(fontsize=7, ncol=2)
    axs[1, 1].text(0.02, 0.02, "timings: shared 4-core machine,\n1 rank, indicative only", transform=axs[1, 1].transAxes,
                   fontsize=7, color="gray")
    fig.suptitle(title); fig.tight_layout(); fig.savefig(fname, dpi=140); plt.close(fig)
    print("wrote", fname)


metric_fig("p", "p", dict(tag="a", N1=8, eps1=0.3), prefix + "_p_scaling.png",
           "E3(a): M1 = f1 eps 0.3 (N=8) -> uniform N=8, abc + B0  (solid: A-route, dashed: B-route)")
metric_fig("eps1", "eps of M1", dict(tag="a", N1=8, p=2), prefix + "_eps_scaling.png",
           "E3(a): p=2, N=8, M1 = f1 eps -> uniform")

# gauge
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for ax, gk in zip(axs, (2.0, 6.0)):
    base = dict(tag="c", gauge_k=gk)
    for op in ("A_pt", "A_int", "A_l2"):
        s = sorted(sel(op=op, **base), key=lambda r: r["p"])
        if s:
            ax.semilogy([r["p"] for r in s], [r["pollution_B"] for r in s], LS[op], marker=MK[op], color=COL[op], label=op)
    ax.set_xlabel("p"); ax.set_xticks([1, 2, 3]); ax.grid(True, which="both", alpha=0.3)
    ax.set_title(f"gauge g=0.5, k={int(gk/2)}pi ({'smooth' if gk == 2 else 'rough'})")
    ax.set_ylabel(r"$\|b_2(A+\nabla\chi)-b_2(A)\|/\|b_2\|$")
axs[0].legend()
fig.tight_layout(); fig.savefig(prefix + "_gauge.png", dpi=140); plt.close(fig); print("wrote gauge")

# nq study
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
for p, c in ((1, "#0072B2"), (2, "#D55E00"), (3, "#009E73")):
    s = sorted(sel(tag="e", op="B_int", p=p), key=lambda r: r["nq"])
    if s:
        axs[0].loglog([r["nq"] for r in s], [r["m2_div_rel_L2"] for r in s], "o--", color=c, label=f"B_int p={p}")
        axs[1].loglog([r["nq"] for r in s], [r["m2_errB_rel"] for r in s], "o--", color=c, label=f"B_int p={p}")
    s = sorted(sel(tag="e", op="A_int", p=p), key=lambda r: r["nq"])
    if s:
        axs[0].loglog([r["nq"] for r in s], [max(r["m2_div_rel_L2"], 1e-17) for r in s], "s-", color=c, label=f"A_int p={p}")
        axs[1].loglog([r["nq"] for r in s], [r["m2_errB_rel"] for r in s], "s-", color=c, label=f"A_int p={p}")
axs[0].set_ylabel(r"$\|D_hb_2\|/(\|B\|/h)$"); axs[1].set_ylabel(r"$\|B_2-B\|/\|B\|$")
for ax in axs:
    ax.set_xlabel("Gauss points per direction nq (sub-edge / sub-face rule)"); ax.grid(True, which="both", alpha=0.3)
axs[0].legend(fontsize=7)
fig.suptitle("E3(e): effect of the quadrature of the integrated transfer (N=8, eps 0.3 -> uniform)")
fig.tight_layout(); fig.savefig(prefix + "_nq.png", dpi=140); plt.close(fig); print("wrote nq")

# cases b, d: error relative to floor per op
cases = [("b", dict(tag="b", N1=8), "(b) f1 0.3 -> f2 0.15, N=8, mod2"),
         ("b", dict(tag="b", N1=16), "(b) same, N=16"),
         ("d", dict(tag="d", N2=12), "(d) N=8 -> N=12"),
         ("d", dict(tag="d", N2=16), "(d) N=8 -> N=16")]
fig, axs = plt.subplots(2, 2, figsize=(11, 7.5), sharey=False)
for ax, (_, base, title) in zip(axs.ravel(), cases):
    w = 0.25
    for i, (p, c) in enumerate(((1, "#56B4E9"), (2, "#E69F00"), (3, "#009E73"))):
        vals = []
        for op in OPS:
            s = sel(op=op, p=p, **base)
            vals.append(s[0]["m2_errB_rel"] / s[0]["ref_A_errB_rel"] if s else np.nan)
        ax.bar(np.arange(len(OPS)) + (i - 1) * w, vals, w * 0.9, color=c, label=f"p={p}")
    ax.set_xticks(range(len(OPS))); ax.set_xticklabels(OPS, fontsize=8)
    ax.set_ylabel("B error / floor (projection of exact B)"); ax.set_title(title, fontsize=9)
    ax.axhline(1, color="gray", lw=1); ax.grid(True, axis="y", alpha=0.3)
    ax.set_yscale("log")
axs[0, 0].legend(fontsize=8)
fig.tight_layout(); fig.savefig(prefix + "_cases.png", dpi=140); plt.close(fig); print("wrote cases")
