#!/usr/bin/env python3
"""Markdown tables for E4.  Usage: scripts/make_tables4.py [csv] (stdout)"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vp_csv import load
csv_path = sys.argv[1] if len(sys.argv) > 1 else "results/exp4_repeat.csv"
rows = load(csv_path)
for _r in rows:   # s3b (A_l2, B_l2 long run) is part of the long-run series
    if _r["tag"] == "s3b":
        _r["tag"] = "s3"
OPS = ["A_pt", "A_int", "A_l2", "B_pt", "B_int", "B_l2", "B_l2c"]


def ser(tag, op, mesh="A"):
    return sorted([r for r in rows if r["tag"] == tag and r["op"] == op and r["mesh"] == mesh], key=lambda r: r["step"])


def at(s, n):
    for r in s:
        if r["step"] == n:
            return r
    return None


def f(v, fmt="%.3e"):
    if v != v:
        return "n/a"
    return fmt % v


def first_fail(tag, op, thr=0.5):
    allr = sorted([r for r in rows if r["tag"] == tag and r["op"] == op], key=lambda r: r["step"])
    for r in allr:
        if r["errB_rel"] > thr:
            return str(int(r["step"]))
    return ">%d" % int(allr[-1]["step"]) if allr else "n/a"


def block(tag, title, ns):
    print(f"### {title}\n")
    head = ["op", "n", "E/E_ref", "B err", "B err / floor", "div (scaled L2)", "&#124;&#124;B_n-B_first&#124;&#124;/&#124;&#124;B&#124;&#124;", "flux err",
            "(H-H_ref)/(&#124;a&#124;&#124;B&#124;)", "&#124;a&#124;/&#124;a_ref&#124;", "grad frac", "Coulomb res.", "first n with B err > 0.5"]
    print("| " + " | ".join(head) + " |")
    print("|" + "---|" * len(head))
    for op in OPS:
        s = ser(tag, op)
        if not s:
            continue
        for n in ns:
            r = at(s, n)
            if r is None:
                continue
            print("| " + " | ".join([op, str(n), f(r["E_over_Eref"], "%.4g"), f(r["errB_rel"]), f(r["errB_over_ref"], "%.3g"), f(r["div_rel_L2"], "%.1e"), f(r["dB_first"], "%.2e"),
                                     f(r["flux_err"], "%.1e"), f(r["H_dev_scaled"], "%+.2e"), f(r["a_over_aref"], "%.4g"), f(r["grad_frac"], "%.3e"), f(r["coulomb"], "%.2e"),
                                     first_fail(tag, op) if n == ns[-1] else ""]) + " |")
    print()


for p in (1, 2, 3):
    if ser(f"s1_p{p}", "A_int"):
        block(f"s1_p{p}", f"Table 4.{p}: p={p}, N=8, abc, no B0, uniform <-> f1 0.3, state on mesh A (even n)", [10, 20, 50, 100])
if ser("s2", "A_int"):
    block("s2", "Table 4.4: p=2, N=8, abc WITH mean field B0 = (0.3,-0.2,0.5)", [10, 20, 50, 100])
if ser("s3", "A_int"):
    block("s3", "Table 4.5: long run, p=2, N=8, no B0", [100, 200, 300, 400])

print("### Table 4.6: cost per remap (transfer + composition, without diagnostics; contended machine) and CG iterations per remap\n")
print("| series | op | s/remap | outer or L2 iterations / remap | inner iterations / remap | plan setup [s] |")
print("|---|---|---|---|---|---|")
for tag in ("s1_p1", "s1_p2", "s1_p3", "s3"):
    for op in OPS:
        s = ser(tag, op)
        if not s:
            continue
        r = s[-1]
        n = r["step"]
        print("| %s | %s | %.3f | %s | %s | %.2f |" % (tag, op, r["t_cum"] / n, f(r["iters_total"] / n, "%.1f") if r["iters_total"] else "-",
                                                     f(r["inner_total"] / n, "%.0f") if r["inner_total"] else "-", r["t_setup"]))
print()

print("### Table 4.7: drift per 100 remaps (state on mesh A), from the long run (n=400, p=2) and the n=100 runs\n")
print("| series | op | (E_n/E_0-1) per 100 | B err growth per 100 | grad frac at end | &#124;a_n&#124;/&#124;a_0&#124; at end |")
print("|---|---|---|---|---|---|")
for tag in ("s1_p1", "s1_p2", "s1_p3", "s2", "s3"):
    for op in OPS:
        s = ser(tag, op)
        if len(s) < 3:
            continue
        r0, r1 = s[0], s[-1]
        n = r1["step"]
        if r1["E_over_E0"] > 1e3 or r1["errB_rel"] > 2:
            print("| %s | %s | diverged (E/E_0 = %.2g, B err = %.2g at n=%d) | | | |" % (tag, op, r1["E_over_E0"], r1["errB_rel"], n))
            continue
        print("| %s | %s | %+.2e | %+.2e | %s | %s |" % (tag, op, (r1["E_over_E0"] - 1) * 100 / n, (r1["errB_rel"] - r0["errB_rel"]) * 100 / n,
                                                         f(r1["grad_frac"], "%.3e"), f(r1["a_over_a0"], "%.4g")))
