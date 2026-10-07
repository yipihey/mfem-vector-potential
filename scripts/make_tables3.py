#!/usr/bin/env python3
"""Markdown tables for E3.  Usage: scripts/make_tables3.py [csv]  (prints to stdout)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vp_csv import load
csv_path = sys.argv[1] if len(sys.argv) > 1 else "results/exp3_remap.csv"
rows = load(csv_path)
OPS = ["A_pt", "A_int", "A_l2", "B_pt", "B_int", "B_l2", "B_l2c"]


def sel(**kw):
    return [r for r in rows if all(r.get(k) == v for k, v in kw.items())]


def e(v, f="%.2e"):
    return "n/a" if v != v else f % v


def table(head, body):
    print("| " + " | ".join(head) + " |")
    print("|" + "---|" * len(head))
    for b in body:
        print("| " + " | ".join(b) + " |")
    print()


print("### Table 3.1 (a): M1 = f1 eps -> uniform, N=8, abc + B0\n")
head = ["eps", "p", "op", "B err", "floor", "err/floor", "div (scaled L2)", "max&#124;D b&#124; abs", "dE/E1", "(E2-Eref)/Eref", "flux err", "t one-shot [s]", "iters (outer/inner)"]
body = []
for eps in (0.1, 0.3, 0.5):
    for p in (1, 2, 3, 4):
        for op in OPS:
            s = sel(tag="a", N1=8, eps1=eps, p=p, op=op)
            if not s: continue
            r = s[0]
            body.append([str(eps), str(p), op, e(r["m2_errB_rel"]), e(r["ref_A_errB_rel"]), "%.2f" % (r["m2_errB_rel"] / r["ref_A_errB_rel"]),
                         e(r["m2_div_rel_L2"]), e(r["m2_div_max"]), e(r["m2_dE_rel1"], "%+.2e"),
                         e((r["m2_energy"] - r["ref_A_energy"]) / r["ref_A_energy"], "%+.2e"), e(r["m2_flux_err"]), "%.2f" % r["t_oneshot"],
                         ("%d/%d" % (r["iters"], r["inner_iters"])) if r["iters"] else "-"])
table(head, body)

print("### Table 3.2 (a): N=16 (p<=2), M1 = f1 eps -> uniform N=16\n")
body = []
for eps in (0.1, 0.3, 0.5):
    for p in (1, 2):
        for op in OPS:
            s = sel(tag="a", N1=16, eps1=eps, p=p, op=op)
            if not s: continue
            r = s[0]
            body.append([str(eps), str(p), op, e(r["m2_errB_rel"]), e(r["ref_A_errB_rel"]), "%.2f" % (r["m2_errB_rel"] / r["ref_A_errB_rel"]),
                         e(r["m2_div_rel_L2"]), e(r["m2_dE_rel1"], "%+.2e"), "%.2f" % r["t_oneshot"],
                         ("%d/%d" % (r["iters"], r["inner_iters"])) if r["iters"] else "-"])
table(["eps", "p", "op", "B err", "floor", "err/floor", "div (scaled L2)", "dE/E1", "t one-shot [s]", "iters"], body)

print("### Table 3.3 (b): f1 0.3 -> f2 0.15 (mod2, no B0)\n")
body = []
for N in (8, 16):
    for p in (1, 2, 3):
        for op in OPS:
            s = sel(tag="b", N1=N, p=p, op=op)
            if not s: continue
            r = s[0]
            body.append([str(N), str(p), op, e(r["m2_errB_rel"]), e(r["ref_A_errB_rel"]), "%.2f" % (r["m2_errB_rel"] / r["ref_A_errB_rel"]),
                         e(r["m2_div_rel_L2"]), e(r["m2_dE_rel1"], "%+.2e"), e((r["m2_energy"] - r["ref_A_energy"]) / r["ref_A_energy"], "%+.2e"),
                         e(r["m2_helicity"], "%+.2e") if r["route"] == "A" else "n/a",
                         e(r["m2_dH_scaled"], "%+.2e"), e(r["m2_a_L2"] / r["m1_a_L2"] if r["route"] == "A" else float("nan"), "%.3f"),
                         "%.2f" % r["t_oneshot"]])
table(["N", "p", "op", "B err", "floor", "err/floor", "div (scaled L2)", "dE/E1", "(E2-Eref)/Eref", "H2 (exact 0)", "(H2-H1)/(&#124;a&#124;&#124;B&#124;)", "&#124;a2&#124;/&#124;a1&#124;", "t [s]"], body)

print("### Table 3.4 (c): gauge pollution, g=0.5 (A1 += grad chi on M1), eps 0.3 -> uniform, N=8\n")
body = []
for gk in (2.0, 6.0):
    for p in (1, 2, 3):
        for op in ("A_pt", "A_int", "A_l2"):
            s = sel(tag="c", gauge_k=gk, p=p, op=op)
            if not s: continue
            r = s[0]
            body.append([("2pi" if gk == 2 else "6pi"), str(p), op, e(r["pollution_B"]), e(r["pollution_A"], "%.3f"), e(r["m2_errB_rel"]),
                         e(r["ref_A_errB_rel"]), e(r["m2_grad_frac"], "%.3f"), e(r["m1_grad_frac"], "%.3f"), e(r["m2_div_rel_L2"])])
table(["k", "p", "op", "&#124;&#124;b2(g)-b2(0)&#124;&#124;/&#124;&#124;b2&#124;&#124;", "&#124;&#124;a2(g)-a2(0)&#124;&#124;/&#124;&#124;a2(0)&#124;&#124;", "B err (with gauge)", "B err floor", "grad frac a2", "grad frac a1", "div"], body)

print("### Table 3.5 (d): resolution change M1 N=8 (f1 0.3) -> uniform N2\n")
body = []
for N2 in (12, 16):
    for p in (1, 2, 3):
        for op in OPS:
            s = sel(tag="d", N2=N2, p=p, op=op)
            if not s: continue
            r = s[0]
            body.append([str(N2), str(p), op, e(r["m2_errB_rel"]), e(r["ref_A_errB_rel"]), e(r["m1_errB_rel"]), "%.2f" % (r["m2_errB_rel"] / r["m1_errB_rel"]),
                         e(r["m2_div_rel_L2"]), e(r["m2_dE_rel1"], "%+.2e"), "%.2f" % r["t_oneshot"]])
table(["N2", "p", "op", "B err (M2)", "floor on M2", "B err on M1", "err2/err1", "div (scaled L2)", "dE/E1", "t [s]"], body)

print("### Table 3.6 (e): quadrature of the integrated transfers (N=8, eps 0.3 -> uniform, abc+B0)\n")
body = []
for p in (1, 2, 3):
    for nq in (2, 3, 4, 6, 8):
        for op in ("A_int", "B_int"):
            s = sel(tag="e", p=p, nq=nq, op=op)
            if not s: continue
            r = s[0]
            body.append([str(p), str(nq), op, e(r["m2_errB_rel"]), e(r["m2_div_rel_L2"]), e(r["m2_div_max"]), e(r["m2_dE_rel1"], "%+.2e"),
                         str(int(r["npoints"])), "%.2f" % r["t_oneshot"]])
table(["p", "nq", "op", "B err", "div (scaled L2)", "max&#124;D b&#124; abs", "dE/E1", "#points", "t [s]"], body)

print("### Table 3.7: cost of the global variants (apply = t_apply, setup = point location etc.), N=8 eps 0.3 -> uniform\n")
body = []
for p in (1, 2, 3, 4):
    for op in OPS:
        s = sel(tag="a", N1=8, eps1=0.3, p=p, op=op)
        if not s: continue
        r = s[0]
        body.append([str(p), op, str(int(r["npoints"])), "%.3f" % r["t_setup"], "%.3f" % r["t_locate"], "%.3f" % r["t_eval"], "%.3f" % r["t_dofs"], "%.3f" % r["t_solve"],
                     "%.3f" % r["t_apply"], ("%d" % r["iters"]) if r["iters"] else "-", ("%d" % r["inner_iters"]) if r["inner_iters"] else "-",
                     e(r["div_before_clean"]) if r["div_before_clean"] == r["div_before_clean"] else "-"])
table(["p", "op", "#points", "t_setup", "(locate)", "t_eval", "t_dofs", "t_solve", "t_apply", "CG iters", "inner iters", "max&#124;D b&#124; before cleaning"], body)
