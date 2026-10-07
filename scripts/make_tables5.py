#!/usr/bin/env python3
"""E5 tables -> results/exp5_tables.md.  Usage: scripts/make_tables5.py [csv] [out.md]"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from vp_csv import load

csv_path = sys.argv[1] if len(sys.argv) > 1 else "results/exp5_ale.csv"
out_path = sys.argv[2] if len(sys.argv) > 2 else "results/exp5_tables.md"
allrows = load(csv_path)
steps = [r for r in allrows if r["kind"] == "step"]
finals = [r for r in allrows if r["kind"] == "final"]

VARS = [("A_pt", "none"), ("A_int", "none"), ("A_int", "coulomb"), ("A_int", "coulomb_pre"), ("A_l2", "none"), ("A_l2", "coulomb"), ("A_l2", "coulomb_pre"),
        ("B_int", "none"), ("B_l2", "none"), ("B_l2c", "none")]
NAME = {("A_int", "coulomb_pre"): "A_int+Coulomb(pre)", ("A_l2", "coulomb_pre"): "A_l2+Coulomb(pre)", ("A_pt", "none"): "A_pt", ("A_int", "none"): "A_int", ("A_int", "coulomb"): "A_int+Coulomb",
        ("A_l2", "none"): "A_l2", ("A_l2", "coulomb"): "A_l2+Coulomb", ("B_int", "none"): "B_int",
        ("B_l2", "none"): "B_l2", ("B_l2c", "none"): "B_l2c"}
CASES = [(2, 0.3), (2, 0.6), (3, 0.3), (3, 0.6)]
NREZ = [1, 2, 5, 20]


def fmt(x, f="{:.2e}"):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "--"
    return f.format(x)


def F(p, eps, nr, op, gauge, b0=0.0):
    for r in finals:
        if (r["p"], r["eps"], r["nrez"], r["op"], r["gauge"], r["b0"]) == (float(p), eps, float(nr), op, gauge, b0) and r["tag"] in ("main", "apt", "extra", "b0", "gpre"):
            return r
    return None


def S(p, eps, nr, op, gauge, b0=0.0):
    return sorted([r for r in steps if (r["p"], r["eps"], r["nrez"], r["op"], r["gauge"], r["b0"]) == (float(p), eps, float(nr), op, gauge, b0)
                   and r["tag"] in ("main", "apt", "extra", "b0", "gpre")], key=lambda r: (r["step"], r["rezone"]))


out = []
w = out.append

# ---------------------------------------------------------------- Lagrangian references
w("### Table 5.1: Lagrangian reference (no rezone), N=8, abc, K=100 steps")
w("")
w("| tag | p | q | eps | B0 | err B(t=0) | peak err B (t) | peak/initial | err B(T) | err A(0) | peak err A | (E-E_ex)/E_ex at 0 | worst over t | min detJ | max &#124;H/H0-1&#124; | max flux err | max scaled div | &#124;&#124;a_T-a_0&#124;&#124;_inf |")
w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for tag in ["lag", "lag_q3", "lag_b0"]:
    for p in [2, 3]:
        for eps in [0.3, 0.6]:
            R = sorted([r for r in steps if r["tag"] == tag and r["p"] == p and r["eps"] == eps], key=lambda r: (r["step"], r["rezone"]))
            if not R:
                continue
            fin = [r for r in finals if r["tag"] == tag and r["p"] == p and r["eps"] == eps][0]
            pk = max(R, key=lambda r: r["errB_rel"])
            Hdev = max([abs(r["H_over_H0"] - 1) for r in R if not math.isnan(r["H_over_H0"])] or [float("nan")])
            w(f"| {tag} | {p} | {int(R[0]['q'])} | {eps} | {int(R[0]['b0'])} | {fmt(R[0]['errB_rel'])} | {fmt(pk['errB_rel'])} ({pk['t']:.2f}) | {pk['errB_rel']/R[0]['errB_rel']:.2f} | "
              f"{fmt(R[-1]['errB_rel'])} | {fmt(R[0]['errA_rel'])} | {fmt(max(r['errA_rel'] for r in R))} | {fmt(R[0]['E_err_rel'])} | "
              f"{fmt(max(abs(r['E_err_rel']) for r in R))} | {min(r['minDetJ'] for r in R):.3f} | {fmt(Hdev)} | {fmt(max(r['flux_err'] for r in R))} | "
              f"{fmt(max(r['div_scaled'] for r in R))} | {fmt(fin['a_change_inf'])} |")
w("")

# ---------------------------------------------------------------- final comparison, headline
def pivot(title, fn, f="{:.2e}"):
    w(f"### {title}")
    w("")
    hdr = "| operator |" + "".join(f" p={p} e={e} r={n} |" for (p, e) in CASES for n in NREZ)
    w(hdr)
    w("|---|" + "---|" * (len(CASES) * len(NREZ)))
    for v in VARS:
        cells = []
        any_ = False
        for (p, e) in CASES:
            for n in NREZ:
                r = F(p, e, n, v[0], v[1])
                x = fn(r) if r is not None else None
                if x is not None and not (isinstance(x, float) and math.isnan(x)):
                    any_ = True
                cells.append(fmt(x, f))
        if any_:
            w(f"| {NAME[v]} |" + "".join(f" {c} |" for c in cells))
    w("")


pivot("Table 5.2: final ||b_T - b_0||/||b_0|| on the uniform mesh at T (rezones r = 1, 2, 5, 20)", lambda r: r["final_dB"])
pivot("Table 5.3: final error vs the exact field ||b_T - B(T)||/||B(T)|| (the floor at T is the t=0 error: p=2 2.28e-2, p=3 1.52e-3)", lambda r: r["final_errB"])
pivot("Table 5.4: energy ratio E_B(T)/E_B(0) - 1", lambda r: r["final_E_ratio"] - 1.0)
pivot("Table 5.5: helicity ratio H(T)/H(0) - 1 (A-route only)", lambda r: r["final_H_ratio"] - 1.0)
pivot("Table 5.6: slice-flux error at T (max over the 6 slice families, B0 = 0)", lambda r: r["final_flux_err"])
pivot("Table 5.7: scaled divergence ||div b_h||/(||B||/h) at T", lambda r: r["final_div"])
pivot("Table 5.8: total transfer wall time (s) of all rezones (evaluator + plan setup + apply + gauge + compose; machine shared)", lambda r: r["t_transfer_total"], "{:.1f}")

# ---------------------------------------------------------------- per-rezone jumps
w("### Table 5.9: error added by one rezone (jump of the relative B error across the transfer; mean over the rezones of a run)")
w("")
w("Rows: operator; columns: (p, eps) at r = 5 and 20 rezones.  Entry: mean(err_post - err_pre) with err = ||b_h - B_pf||/||B_pf|| measured on the source mesh right before and on the uniform target right after the transfer; in parentheses: mean err_post/err_pre; in brackets: sum of the jumps over all rezones of the run (note that for r = 5 the rezone at t=0.6 is an exact identity because s(0.6)=s(0.4): 4 genuine transfers; r=20 has none, r=2 none).")
w("")
cols = [(p, e, n) for (p, e) in CASES for n in (2, 5, 20)]
VARS_J = VARS
w("| operator |" + "".join(f" p={p} e={e} r={n} |" for (p, e, n) in cols))
w("|---|" + "---|" * len(cols))
for v in VARS:
    cells = []
    any_ = False
    for (p, e, n) in cols:
        R = S(p, e, n, v[0], v[1])
        pre = {r["step"]: r for r in R if r["rezone"] == 0}
        jumps = []
        ratios = []
        for r in R:
            if r["rezone"] == 1 and r["step"] in pre:
                jumps.append(r["errB_rel"] - pre[r["step"]]["errB_rel"])
                ratios.append(r["errB_rel"] / pre[r["step"]]["errB_rel"])
        if jumps:
            any_ = True
            cells.append(f"{np.mean(jumps):+.1e} ({np.mean(ratios):.2f}) [sum {np.sum(jumps):+.1e}]")
        else:
            cells.append("--")
    if any_:
        w(f"| {NAME[v]} |" + "".join(f" {c} |" for c in cells))
w("")

# ---------------------------------------------------------------- scaling with the number of rezones
w("### Table 5.10: does the accumulated difference scale with the number of rezones? (final ||b_T-b_0||/||b_0||, ratio to r = 2, and fitted exponent over r = 2, 5, 20)")
w("")
w("A rezone at t = T from the (uniform) Lagrangian mesh is the identity (r = 1 column of 5.2 is the roundoff of the transfer), so r = 2 is the smallest case with genuine transfers.  Linear growth in r would give ratios 2.5 (r=5) and 10 (r=20) and exponent 1; random-walk growth (errors adding in quadrature) gives exponent 0.5; a fixed per-rezone error that does not accumulate (projection idempotence) gives exponent 0.")
w("")
w("| (p, eps) | operator | r=1 | r=2 | r=5 | r=20 | r5/r2 | r20/r2 | exponent (log-log fit r=2,5,20) |")
w("|---|---|---|---|---|---|---|---|---|")
for (p, e) in CASES:
    for v in VARS:
        vals = {n: F(p, e, n, v[0], v[1]) for n in NREZ}
        if any(vals[n] is None for n in (2, 5, 20)):
            continue
        d = {n: vals[n]["final_dB"] for n in NREZ if vals[n] is not None}
        ex = np.polyfit(np.log([2, 5, 20]), np.log([d[2], d[5], d[20]]), 1)[0]
        w(f"| p={p}, eps={e} | {NAME[v]} | {fmt(d.get(1))} | {fmt(d[2])} | {fmt(d[5])} | {fmt(d[20])} | {d[5]/d[2]:.2f} | {d[20]/d[2]:.2f} | {ex:.2f} |")
w("")

# ---------------------------------------------------------------- gauge
w("### Table 5.11: Coulomb gauge fixing at each rezone (A_int/A_l2 + Coulomb after the transfer; (pre) = on the deformed source mesh before the transfer): mean over rezones")
w("")
w("| p | eps | r | op | ||P_grad a||/||a|| before fix (for (pre): on the deformed source mesh) | Coulomb residual before | after | gauge CG iters | gauge time per rezone (s) | rezone total (s) |")
w("|---|---|---|---|---|---|---|---|---|---|")
for (p, e) in CASES:
    for n in (2, 5, 20):
        for v in [("A_int", "coulomb"), ("A_l2", "coulomb"), ("A_int", "coulomb_pre"), ("A_l2", "coulomb_pre")]:
            R = [r for r in S(p, e, n, v[0], v[1]) if r["rezone"] == 1 and r["gauge_iters"] >= 0]
            if R:
                w(f"| {p} | {e} | {n} | {NAME[v]} | {np.mean([r['gauge_grad_before'] for r in R]):.2e} | {np.mean([r['gauge_resid_before'] for r in R]):.2e} | "
                  f"{np.mean([r['gauge_resid_after'] for r in R]):.2e} | {np.mean([r['gauge_iters'] for r in R]):.0f} | {np.mean([r['t_gauge'] for r in R]):.2f} | {np.mean([r['t_rezone'] for r in R]):.2f} |")
w("")

# ---------------------------------------------------------------- b0 series
w("### Table 5.12: mean field B0 = (0.3,-0.2,0.5) on (p=2, eps=0.3, 5 rezones; B-route; A_int Lagrangian reference with B0 for comparison)")
w("")
w("| operator | rezones | ||b_T-b_0||/||b_0|| | err B(T) | E(T)/E(0)-1 | slice-flux error at T | scaled div | max flux err over t | transfer time (s) |")
w("|---|---|---|---|---|---|---|---|---|")
for r in finals:
    if r["b0"] == 1.0:
        R = [x for x in steps if (x["p"], x["eps"], x["nrez"], x["op"], x["gauge"], x["b0"], x["tag"]) == (r["p"], r["eps"], r["nrez"], r["op"], r["gauge"], r["b0"], r["tag"])]
        w(f"| {r['op']} | {int(r['nrez'])} | {fmt(r['final_dB'])} | {fmt(r['final_errB'])} | {fmt(r['final_E_ratio']-1)} | {fmt(r['final_flux_err'])} | {fmt(r['final_div'])} | {fmt(max(x['flux_err'] for x in R))} | {fmt(r['t_transfer_total'], '{:.1f}')} |")
w("")

# ---------------------------------------------------------------- peak during cycle (pre-rezone state at s=1 etc.)
w("### Table 5.13: B error during the cycle (max over all logged rows of the run), relative to the Lagrangian peak of the same (p, eps)")
w("")
w("| (p, eps) | Lagrangian peak | " + " | ".join(f"{NAME[v]} r=5 / r=20" for v in VARS if v != ("A_pt", "none") and v != ("A_l2", "none")) + " |")
w("|---|---|" + "---|" * (len([v for v in VARS if v != ("A_pt", "none") and v != ("A_l2", "none")])))
for (p, e) in CASES:
    L = [r for r in steps if r["tag"] == "lag" and r["p"] == p and r["eps"] == e]
    if not L:
        continue
    lp = max(r["errB_rel"] for r in L)
    cells = []
    for v in VARS:
        if v in (("A_pt", "none"), ("A_l2", "none")):
            continue
        c = []
        for n in (5, 20):
            R = S(p, e, n, v[0], v[1])
            c.append(fmt(max(r["errB_rel"] for r in R) if R else None))
        cells.append(" / ".join(c))
    w(f"| p={p}, eps={e} | {fmt(lp)} | " + " | ".join(cells) + " |")
w("")

# ---------------------------------------------------------------- q = 3 control
w("### Table 5.14: geometry control at p=3, eps=0.6, 5 rezones: q=2 (main runs) versus q=3 (tag q3)")
w("")
w("| operator | q | ||b_T-b_0||/||b_0|| | err B(T) | max err B over t | mean jump per rezone | E(T)/E(0)-1 | transfer time (s) |")
w("|---|---|---|---|---|---|---|---|")
for v in [("B_l2", "none"), ("B_l2c", "none"), ("A_int", "coulomb"), ("A_l2", "coulomb")]:
    for q in (2, 3):
        if q == 2:
            fr = F(3, 0.6, 5, v[0], v[1])
            R = S(3, 0.6, 5, v[0], v[1])
        else:
            fl = [r for r in finals if r["tag"] == "q3" and r["op"] == v[0] and r["gauge"] == v[1]]
            fr = fl[0] if fl else None
            R = sorted([r for r in steps if r["tag"] == "q3" and r["op"] == v[0] and r["gauge"] == v[1]], key=lambda r: (r["step"], r["rezone"]))
        if fr is None:
            continue
        pre = {r["step"]: r for r in R if r["rezone"] == 0}
        jumps = [r["errB_rel"] - pre[r["step"]]["errB_rel"] for r in R if r["rezone"] == 1 and r["step"] in pre]
        w(f"| {NAME[v]} | {q} | {fmt(fr['final_dB'])} | {fmt(fr['final_errB'])} | {fmt(max(r['errB_rel'] for r in R))} | {np.mean(jumps):+.1e} | {fmt(fr['final_E_ratio']-1)} | {fmt(fr['t_transfer_total'], '{:.1f}')} |")
w("")

# ---------------------------------------------------------------- floors
floor_path = "results/exp5_floors.csv"
if os.path.exists(floor_path):
    FL = load(floor_path)
    w("### Table 5.15: error right after each rezone versus the projection floor of the exact push-forward on the new uniform mesh")
    w("")
    w("floor = ||Pi_RT^int B_pf(t_r) - B_pf(t_r)||/||B_pf|| on the new uniform mesh (best any transfer can do; operator independent, `-floors` mode), floor(src) = same on the deformed source mesh just before the rezone. Columns: relative B error right after the transfer (and in brackets error/floor). r=5: the rezone at t=0.6 is an identity (s(0.6)=s(0.4)).")
    w("")
    opsT = [("B_l2c", "none"), ("B_l2", "none"), ("B_int", "none"), ("A_int", "coulomb"), ("A_int", "none"), ("A_l2", "coulomb")]
    for (p, e) in CASES:
        for n in (2, 5):
            rows_f = sorted([r for r in FL if r["p"] == p and r["eps"] == e and r["nrez"] == n and r["q"] == 2], key=lambda r: r["rezone_index"])
            if not rows_f:
                continue
            w(f"**p={p}, eps={e}, r={n} rezones (q=2)**")
            w("")
            w("| rezone | t | min detJ (source) | floor(src) | floor(new) | " + " | ".join(NAME[v] for v in opsT) + " |")
            w("|---|---|---|---|---|" + "---|" * (len(opsT)))
            for r in rows_f:
                cells = []
                for v in opsT:
                    R = [x for x in S(p, e, n, v[0], v[1]) if x["rezone"] == 1 and x["step"] == r["step"]]
                    cells.append(f"{R[0]['errB_rel']:.2e} ({R[0]['errB_rel']/r['floor_post']:.2f})" if R else "--")
                w(f"| {int(r['rezone_index'])} | {r['t']:.2f} | {r['minDetJ_src']:.3f} | {r['floor_pre']:.2e} | {r['floor_post']:.2e} | " + " | ".join(cells) + " |")
            w("")
    rows_q3 = sorted([r for r in FL if r["p"] == 3 and r["eps"] == 0.6 and r["nrez"] == 5 and r["q"] == 3], key=lambda r: r["rezone_index"])
    if rows_q3:
        w("**p=3, eps=0.6, r=5 rezones, q=3 control**")
        w("")
        w("| rezone | t | min detJ (source) | floor(src) | floor(new) | B_l2c | B_l2 | A_int+Coulomb | A_l2+Coulomb |")
        w("|---|---|---|---|---|---|---|---|---|")
        for r in rows_q3:
            cells = []
            for op, g in [("B_l2c", "none"), ("B_l2", "none"), ("A_int", "coulomb"), ("A_l2", "coulomb")]:
                R = [x for x in steps if x["tag"] == "q3" and x["op"] == op and x["gauge"] == g and x["rezone"] == 1 and x["step"] == r["step"]]
                cells.append(f"{R[0]['errB_rel']:.2e} ({R[0]['errB_rel']/r['floor_post']:.2f})" if R else "--")
            w(f"| {int(r['rezone_index'])} | {r['t']:.2f} | {r['minDetJ_src']:.3f} | {r['floor_pre']:.2e} | {r['floor_post']:.2e} | " + " | ".join(cells) + " |")
        w("")

with open(out_path, "w") as f:
    f.write("\n".join(out) + "\n")
print("wrote", out_path)
