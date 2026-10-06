#!/usr/bin/env python3
"""E2 plots/tables.  Usage: scripts/plot_exp2.py [csv] [outprefix]
Writes results/exp2_distort.png (B error, divergence, energy error, min detJ vs
eps, N=8, field abc, f1, q=2) and results/exp2_q_effect.png (effect of q and of
the variant), prints tables."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vp_csv import load, select

csv_path = sys.argv[1] if len(sys.argv) > 1 else "results/exp2_distort.csv"
prefix = sys.argv[2] if len(sys.argv) > 2 else "results/exp2"
allrows = load(csv_path)
rows = [r for r in allrows if r["valid"] == 1]
colors = {1: "#0072B2", 2: "#D55E00", 3: "#009E73", 4: "#CC79A7"}


def series(**kw):
    return sorted(select(rows, **kw), key=lambda r: r["eps"])


def panel(ax, key, ylabel, logy=True, abs_=False, **base):
    for p in (1, 2, 3, 4):
        for proj, ls, mk in (("int", "-", "o"), ("pt", "--", "s")):
            s = series(p=p, proj=proj, **base)
            if not s:
                continue
            e = [r["eps"] for r in s]
            v = [abs(r[key]) if abs_ else r[key] for r in s]
            (ax.semilogy if logy else ax.plot)(e, v, ls, marker=mk, ms=3.5, color=colors[p],
                                              label=f"p={p} {proj}")
    ax.set_xlabel("eps"); ax.set_ylabel(ylabel); ax.grid(True, which="both", alpha=0.3)


base = dict(N=8, q=2, fvariant=1, field="abc", b0=0)
fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
panel(axs[0, 0], "errB_rel", r"$\|B_h-B\|/\|B\|$", **base); axs[0, 0].set_title("B error")
panel(axs[0, 1], "div_rel_L2", r"$\|D_hB_h\|_{L^2}/(\|B\|/h)$", **base); axs[0, 1].set_title("divergence (relative)")
panel(axs[1, 0], "energy_rel_err", "|energy rel. error|", abs_=True, **base); axs[1, 0].set_title("energy error")
for p in (1, 2, 3, 4):   # min detJ depends only on mesh (q, N, variant)
    pass
for q, mk in ((1, "o"), (2, "s")):
    for fv, ls in ((1, "-"), (2, "--"), (3, ":")):
        s = series(p=1, N=8, q=q, fvariant=fv, field="abc", b0=0, proj="int")
        sall = sorted([r for r in allrows if r["p"] == 1 and r["N"] == 8 and r["q"] == q and
                       r["fvariant"] == fv and r["field"] == "abc" and r["b0"] == 0 and r["proj"] == "int"],
                      key=lambda r: r["eps"])
        if sall:
            axs[1, 1].plot([r["eps"] for r in sall], [r["min_detJ"] for r in sall], ls, marker=mk, ms=3.5,
                           label=f"q={q} f{fv}")
axs[1, 1].set_title("min detJ (normalised, quad pts + corners), N=8"); axs[1, 1].set_xlabel("eps")
axs[1, 1].grid(True, alpha=0.3); axs[1, 1].legend(fontsize=7)
axs[0, 0].legend(fontsize=7, ncol=2)
fig.suptitle("E2: ABC field, f1 deformation, N=8, q=2 (solid: integrated proj, dashed: pointwise)")
fig.tight_layout(); fig.savefig(prefix + "_distort.png", dpi=140)
print("wrote", prefix + "_distort.png")

# q / variant effect at fixed p
fig, axs = plt.subplots(1, 3, figsize=(14, 4.5))
for ax, p in zip(axs, (1, 2, 3)):
    for q, c in ((1, "#0072B2"), (2, "#D55E00"), (3, "#009E73")):
        for fv, ls in ((1, "-"), (2, "--"), (3, ":")):
            s = series(p=p, N=8, q=q, fvariant=fv, field="abc", b0=0, proj="int")
            if s:
                ax.semilogy([r["eps"] for r in s], [r["errB_rel"] for r in s], ls, color=c, marker="o",
                            ms=3, label=f"q={q} f{fv}")
    ax.set_title(f"p={p}: rel. B error (N=8, int)"); ax.set_xlabel("eps"); ax.grid(True, which="both", alpha=0.3)
axs[0].legend(fontsize=7, ncol=2)
fig.tight_layout(); fig.savefig(prefix + "_q_effect.png", dpi=140)
print("wrote", prefix + "_q_effect.png")

# ---- tables ----
print("\n### Table A: f1, N=8, q=2, abc, integrated projection\n")
print("| p | eps | min detJ | min/max detJ | errB rel | div max | div rel L2 | E rel err | helicity err | flux err | min flux |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for p in (1, 2, 3, 4):
    for r in series(p=p, N=8, q=2, fvariant=1, field="abc", b0=0, proj="int"):
        print(f"| {p} | {r['eps']:.2f} | {r['min_detJ']:.3f} | {r['min_ratio']:.3f} | {r['errB_rel']:.3e} | {r['div_max']:.2e} | "
              f"{r['div_rel_L2']:.1e} | {r['energy_rel_err']:+.2e} | {r['helicity_abs_err']:.2e} | {r['flux_err_max']:.1e} | |")
print("\n### Invalid / degenerate meshes (valid=0)\n")
for r in allrows:
    if r["valid"] == 0:
        print(f"- p={int(r['p'])} N={int(r['N'])} q={int(r['q'])} f{int(r['fvariant'])} eps={r['eps']}: non-positive detJ (min {r['min_detJ']:.3g}, n_neg={int(r['n_neg_detJ'])})")
print("\n### B0 series (abc + B0): exactness of Pi_RT(B0) and slice flux\n")
print("| proj | p | q | f | eps | b0proj err L2 | b0proj div max | flux err max | flux curl part | errB rel |")
print("|---|---|---|---|---|---|---|---|---|---|")
for proj in ("int", "pt"):
    for p in (1, 2, 3, 4):
        for r in sorted(select(rows, N=8, p=p, b0=1, proj=proj, fvariant=1), key=lambda r: r["eps"]):
            print(f"| {proj} | {p} | {int(r['q'])} | f1 | {r['eps']:.2f} | {r['b0proj_err_L2']:.2e} | {r['b0proj_div_max']:.2e} | "
                  f"{r['flux_err_max']:.2e} | {r['flux_curlpart_max']:.2e} | {r['errB_rel']:.3e} |")
print("\nmax over valid runs: max|C_h G_h| = %.2e, max div_rel_L2 = %.2e, max div_max/(|B|/h) = %.2e, max flux_curlpart = %.2e"
      % (max(r["max_CG"] for r in rows), max(r["div_rel_L2"] for r in rows),
         max(r["div_rel_max"] for r in rows), max(r["flux_curlpart_max"] for r in rows)))
