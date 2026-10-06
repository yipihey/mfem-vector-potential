#!/usr/bin/env python3
"""E1 plots/tables.  Usage: scripts/plot_exp1.py [csv] [png]
Reads results/exp1_static.csv, writes results/exp1_convergence.png and prints
tables of observed convergence rates and divergence maxima."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from vp_csv import load, select, rates

csv_path = sys.argv[1] if len(sys.argv) > 1 else "results/exp1_static.csv"
png = sys.argv[2] if len(sys.argv) > 2 else "results/exp1_convergence.png"
rows = [r for r in load(csv_path) if r["b0"] == 0]
colors = {1: "#0072B2", 2: "#D55E00", 3: "#009E73", 4: "#CC79A7"}
fields = ["abc", "mod2"]

fig, axs = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True)
for j, fld in enumerate(fields):
    for i, (key, lab) in enumerate([("errB_L2", r"$\|B_h-B\|_{L^2}$"),
                                    ("errA_L2", r"$\|A_h-A\|_{L^2}$")]):
        ax = axs[i, j]
        for p in (1, 2, 3, 4):
            for proj, ls, mk in (("pt", "-", "o"), ("int", "--", "s")):
                sel = sorted(select(rows, field=fld, p=p, proj=proj), key=lambda r: r["N"])
                if not sel:
                    continue
                h = [r["h"] for r in sel]
                e = [r[key] for r in sel]
                ax.loglog(h, e, ls, marker=mk, ms=4, color=colors[p],
                          label=f"p={p} {proj}")
            sel = sorted(select(rows, field=fld, p=p, proj="pt"), key=lambda r: r["N"])
            if sel:   # reference slope anchored at the coarsest pt point:
                # h^p for B (RT_{p-1}), h^(p+1) for A (ND_p contains P_p)
                h = np.array([r["h"] for r in sel])
                e0 = sel[0][key]
                sl = p if key == "errB_L2" else p + 1
                ax.loglog(h, e0 * (h / h[0]) ** sl, ":", color=colors[p], lw=1, alpha=0.8)
        ax.set_title(f"{fld}: {lab}")
        ax.grid(True, which="both", alpha=0.3)
        if i == 1:
            ax.set_xlabel("h = 1/N")
        if j == 0:
            ax.set_ylabel("L2 error")
axs[0, 0].legend(fontsize=7, ncol=2)
fig.suptitle("E1: convergence of B_h = C_h a_h (dotted: h^p) and A_h (dotted: h^(p+1)); pt and int coincide on the undeformed mesh")
fig.tight_layout()
fig.savefig(png, dpi=140)
print("wrote", png)

print("\n### Observed convergence rates (between successive N)\n")
print("| field | proj | p | N | h | err B | rate B | err A | rate A | max div_h B | rel div (L2) |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for fld in fields:
    for proj in ("pt", "int"):
        for p in (1, 2, 3, 4):
            sel = sorted(select(rows, field=fld, p=p, proj=proj), key=lambda r: r["N"])
            if not sel:
                continue
            h = [r["h"] for r in sel]
            rb = rates(h, [r["errB_L2"] for r in sel])
            ra = rates(h, [r["errA_L2"] for r in sel])
            for r, x, y in zip(sel, rb, ra):
                print(f"| {fld} | {proj} | {p} | {int(r['N'])} | {r['h']:.4g} | {r['errB_L2']:.3e} | "
                      f"{x:.2f} | {r['errA_L2']:.3e} | {y:.2f} | {r['div_max']:.2e} | {r['div_rel_L2']:.1e} |")
print("\nmax over all runs: div_max = %.3e, div_rel_L2 = %.3e, max|D_h C_h| scaled by h^3 = %.3e, max|C_h G_h| = %.3e, flux_err_max = %.3e"
      % (max(r["div_max"] for r in rows), max(r["div_rel_L2"] for r in rows),
         max(r["max_DC_scaled"] for r in rows), max(r["max_CG"] for r in rows),
         max(r["flux_err_max"] for r in load(csv_path))))
