#!/usr/bin/env python3
"""Markdown tables for results/exp1_summary.md and results/exp2_summary.md.
Usage: scripts/make_tables.py exp1|exp2 > tables.md"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vp_csv import load, select, rates


def exp1():
    rows = load(sys.argv[2] if len(sys.argv) > 2 else "results/exp1_static.csv")
    base = [r for r in rows if r["b0"] == 0]
    print("### Table 1.1 Convergence (field, p, N): errors, observed order (between successive N)\n")
    print("pt and int rows are identical to all printed digits on the undeformed mesh (A_i is independent of x_i for both fields, so line integrals along coordinate edges equal point values times h); only `pt` is listed.\n")
    print("| field | p | N | h | ND dofs | err A (L2) | order A | err B (L2) | order B | rel err B |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for fld in ("abc", "mod2"):
        for p in (1, 2, 3, 4):
            sel = sorted(select(base, field=fld, p=p, proj="pt"), key=lambda r: r["N"])
            h = [r["h"] for r in sel]
            ra = rates(h, [r["errA_L2"] for r in sel]); rb = rates(h, [r["errB_L2"] for r in sel])
            for r, a, b in zip(sel, ra, rb):
                print(f"| {fld} | {p} | {int(r['N'])} | {r['h']:.4g} | {int(r['ndofs_nd'])} | {r['errA_L2']:.3e} | {a:.2f} | "
                      f"{r['errB_L2']:.3e} | {b:.2f} | {r['errB_rel']:.2e} |")
    print("\n### Table 1.2 Exactness of the discrete structure (field abc, proj pt)\n")
    print("max_DC_scaled = max|D_h C_h| * hmin^3 (D_h returns L2 *values* of div and so carries a 1/detJ = N^3 factor); div columns are for b = C_h a_h; flux error = max over the 6 slice families |flux - B0_i|.\n")
    print("| p | N | max|C_h G_h| | max|D_h C_h| | max_DC_scaled | max&#124;D_h b&#124; | &#124;D_h b&#124;_L2 | rel div (L2) | rel div (max) | flux err | B0-proj err (L2) | B0-proj div max |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for p in (1, 2, 3, 4):
        for r in sorted(select(base, field="abc", p=p, proj="pt"), key=lambda r: r["N"]):
            print(f"| {p} | {int(r['N'])} | {r['max_CG']:.1e} | {r['max_DC']:.1e} | {r['max_DC_scaled']:.1e} | {r['div_max']:.1e} | {r['div_L2']:.1e} | "
                  f"{r['div_rel_L2']:.1e} | {r['div_rel_max']:.1e} | {r['flux_err_max']:.1e} | {r['b0proj_err_L2']:.1e} | {r['b0proj_div_max']:.1e} |")
    print("\n### Table 1.3 Energy and helicity (abc: E_exact = 0.5 k^2 (a^2+b^2+c^2) = %.6f, H_exact = k(a^2+b^2+c^2) = %.6f; mod2: E_exact = %.6f, H_exact = 0)\n"
          % (select(base, field="abc")[0]["energy_exact"], select(base, field="abc")[0]["helicity_exact"],
             select(base, field="mod2")[0]["energy_exact"]))
    print("| field | p | N | energy rel err | order | helicity abs err | order | matrix vs matrix-free energy diff |")
    print("|---|---|---|---|---|---|---|---|")
    for fld in ("abc", "mod2"):
        for p in (1, 2, 3, 4):
            sel = sorted(select(base, field=fld, p=p, proj="pt"), key=lambda r: r["N"])
            h = [r["h"] for r in sel]
            re = rates(h, [abs(r["energy_rel_err"]) for r in sel]); rh = rates(h, [max(r["helicity_abs_err"], 1e-300) for r in sel])
            for r, a, b in zip(sel, re, rh):
                md = r["energy_mass_diff"]
                print(f"| {fld} | {p} | {int(r['N'])} | {r['energy_rel_err']:+.2e} | {a:.2f} | {r['helicity_abs_err']:.2e} | {b:.2f} | "
                      f"{('%.1e' % md) if md == md else 'n/a'} |")
    print("\n### Table 1.4 Cost (field abc, proj pt, 1 rank; times in seconds, mean of repeated C_h applications)\n")
    print("| p | N | H1 dofs | ND dofs | RT dofs | L2 dofs | t assemble G,C,D | t C_h apply | CG its (Jacobi, ND mass) | RSS MB | t total |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for p in (1, 2, 3, 4):
        for r in sorted(select(base, field="abc", p=p, proj="pt"), key=lambda r: r["N"]):
            print(f"| {p} | {int(r['N'])} | {int(r['ndofs_h1'])} | {int(r['ndofs_nd'])} | {int(r['ndofs_rt'])} | {int(r['ndofs_l2'])} | "
                  f"{r['t_ops']:.3f} | {r['t_curl']:.2e} | {int(r['cg_iters'])} | {r['rss_mb']:.0f} | {r['t_total']:.1f} |")
    b0 = [r for r in rows if r["b0"] == 1]
    if b0:
        print("\n### Table 1.5 Mean field B0 = (0.3,-0.2,0.5), abc (flux error = max |flux - B0_i| over 6 slice families; curl part = max |flux of C_h a|)\n")
        print("| p | N | proj | flux err | flux of curl part | max|D_h b| | B0-proj L2 err |")
        print("|---|---|---|---|---|---|---|")
        for r in sorted(b0, key=lambda r: (r["p"], r["N"], r["proj"])):
            print(f"| {int(r['p'])} | {int(r['N'])} | {r['proj']} | {r['flux_err_max']:.1e} | {r['flux_curlpart_max']:.1e} | {r['div_max']:.1e} | {r['b0proj_err_L2']:.1e} |")


def exp2():
    allrows = load(sys.argv[2] if len(sys.argv) > 2 else "results/exp2_distort.csv")
    rows = [r for r in allrows if r["valid"] == 1]
    print("### Table 2.1 Distortion scan (f1, N=8, q=2, field abc, no B0, integrated projection)\n")
    print("min detJ is normalised to the undeformed value (quad points + corners); ratio = min over elements of (min detJ)/(max detJ) inside one element.\n")
    print("| p | eps | min detJ | min ratio | rel err B | rel err A | max&#124;D_h b&#124; | rel div L2 | energy rel err | helicity abs err | flux err | max&#124;D_hC_h&#124; scaled |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for p in (1, 2, 3, 4):
        for r in sorted(select(rows, N=8, q=2, fvariant=1, field="abc", b0=0, proj="int", p=p), key=lambda r: r["eps"]):
            print(f"| {p} | {r['eps']:.2f} | {r['min_detJ']:.3f} | {r['min_ratio']:.3f} | {r['errB_rel']:.3e} | {r['errA_rel']:.3e} | {r['div_max']:.1e} | {r['div_rel_L2']:.1e} | "
                  f"{r['energy_rel_err']:+.2e} | {r['helicity_abs_err']:.1e} | {r['flux_err_max']:.1e} | {r['max_DC_scaled']:.1e} |")
    print("\n### Table 2.2 Invalid / degenerate meshes (valid=0: some detJ <= 0)\n")
    print("| p | N | q | variant | eps | min detJ | #points detJ<=0 |")
    print("|---|---|---|---|---|---|---|")
    seen = set()
    for r in allrows:
        if r["valid"] == 0:
            key = (r["N"], r["q"], r["fvariant"], r["eps"])
            if key in seen: continue
            seen.add(key)
            print(f"| {int(r['p'])} | {int(r['N'])} | {int(r['q'])} | f{int(r['fvariant'])} | {r['eps']} | {r['min_detJ']:.3f} | {int(r['n_neg_detJ'])} |")
    print("\n### Table 2.3 Geometric order and variant (N=8, abc, int): relative B error at selected eps\n")
    print("| p | q | variant | eps | min detJ | rel err B | div rel L2 | energy rel err |")
    print("|---|---|---|---|---|---|---|---|")
    for p in (1, 2, 3, 4):
        for q in sorted(set(int(r['q']) for r in rows if r['p'] == p and r['N'] == 8)):
            for fv, e in ((1, 0.3), (2, 0.3), (3, 0.3), (1, 0.6)):
                for r in select(rows, N=8, p=p, q=q, fvariant=fv, eps=e, field="abc", b0=0, proj="int"):
                    print(f"| {p} | {q} | f{fv} | {e} | {r['min_detJ']:.3f} | {r['errB_rel']:.3e} | {r['div_rel_L2']:.1e} | {r['energy_rel_err']:+.2e} |")
    print("\n### Table 2.4 Pointwise vs integrated projection (N=8, q=2, f1, abc)\n")
    print("| p | eps | rel err B (pt) | rel err B (int) | rel err A (pt) | rel err A (int) | E err (pt) | E err (int) |")
    print("|---|---|---|---|---|---|---|---|")
    for p in (1, 2, 3, 4):
        for e in (0.1, 0.2, 0.3, 0.5):
            a = select(rows, N=8, q=2, p=p, fvariant=1, eps=e, field="abc", b0=0, proj="pt")
            b = select(rows, N=8, q=2, p=p, fvariant=1, eps=e, field="abc", b0=0, proj="int")
            if a and b:
                a, b = a[0], b[0]
                print(f"| {p} | {e} | {a['errB_rel']:.3e} | {b['errB_rel']:.3e} | {a['errA_rel']:.3e} | {b['errA_rel']:.3e} | {a['energy_rel_err']:+.1e} | {b['energy_rel_err']:+.1e} |")
    b0 = [r for r in rows if r["b0"] == 1]
    print("\n### Table 2.5 Mean field B0 = (0.3,-0.2,0.5), abc, N=8, q=2 (slice flux; B0 projection exactness)\n")
    print("| proj | p | variant | eps | b0-proj err L2 | b0-proj max&#124;D_h&#124; | flux err (max of 6 slices) | flux of curl part | max&#124;D_h b&#124; | rel err B |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for r in sorted(b0, key=lambda r: (r["proj"], r["p"], r["fvariant"], r["eps"])):
        print(f"| {r['proj']} | {int(r['p'])} | f{int(r['fvariant'])} | {r['eps']:.1f} | {r['b0proj_err_L2']:.1e} | {r['b0proj_div_max']:.1e} | {r['flux_err_max']:.1e} | {r['flux_curlpart_max']:.1e} | {r['div_max']:.1e} | {r['errB_rel']:.2e} |")
    m2 = select(rows, field="mod2")
    if m2:
        print("\n### Table 2.6 mod2 (multi-mode, non-Beltrami), N=8, q=2, f1, int\n")
        print("| p | eps | rel err B | rel div L2 | energy rel err | helicity abs err (exact 0) |")
        print("|---|---|---|---|---|---|")
        for r in sorted(m2, key=lambda r: (r["p"], r["eps"])):
            print(f"| {int(r['p'])} | {r['eps']:.1f} | {r['errB_rel']:.3e} | {r['div_rel_L2']:.1e} | {r['energy_rel_err']:+.2e} | {r['helicity_abs_err']:.1e} |")
    print("\n### Table 2.7 Cost vs eps (field abc, q=2, f1, int; N=8)\n")
    print("| p | eps | ND dofs | t assemble | t C_h apply | RSS MB | t total |")
    print("|---|---|---|---|---|---|---|")
    for p in (1, 2, 3, 4):
        for e in (0.0, 0.3):
            for r in select(rows, N=8, q=2, p=p, fvariant=1, eps=e, field="abc", b0=0, proj="int"):
                print(f"| {p} | {e} | {int(r['ndofs_nd'])} | {r['t_ops']:.3f} | {r['t_curl']:.2e} | {r['rss_mb']:.0f} | {r['t_total']:.1f} |")
    print("\n### Table 2.8 h-convergence on deformed meshes (f1, abc, int, no B0, N=8 -> 16; p <= 2)\n")
    print("| p | q | eps | min detJ (N=16) | err B N=8 | err B N=16 | order B | err A N=8 | err A N=16 | order A | max&#124;D_h b&#124; N=16 |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    import math
    for p in (1, 2):
        for q in (1, 2):
            for e in (0.0, 0.1, 0.2, 0.3, 0.5):
                a = select(rows, N=8, q=q, p=p, fvariant=1, eps=e, field="abc", b0=0, proj="int")
                b = select(rows, N=16, q=q, p=p, fvariant=1, eps=e, field="abc", b0=0, proj="int")
                if a and b:
                    a, b = a[0], b[0]
                    ob = math.log(a["errB_L2"] / b["errB_L2"]) / math.log(2); oa = math.log(a["errA_L2"] / b["errA_L2"]) / math.log(2)
                    print(f"| {p} | {q} | {e} | {b['min_detJ']:.3f} | {a['errB_L2']:.3e} | {b['errB_L2']:.3e} | {ob:.2f} | {a['errA_L2']:.3e} | {b['errA_L2']:.3e} | {oa:.2f} | {b['div_max']:.1e} |")
    good = [r for r in rows if not (r["b0"] == 1 and r["proj"] == "pt")]
    bad = [r for r in rows if (r["b0"] == 1 and r["proj"] == "pt")]
    print("\nOverall, all valid runs except the pointwise-B0 series (%d rows): max|C_h G_h| = %.1e; max rel div L2 = %.1e; max div_max/(|B|/h) = %.1e; max max|D_h b| = %.1e; max flux of curl part = %.1e; max flux err = %.1e; max (max_DC_scaled) = %.1e"
          % (len(good), max(r["max_CG"] for r in good), max(r["div_rel_L2"] for r in good), max(r["div_rel_max"] for r in good),
             max(r["div_max"] for r in good), max(r["flux_curlpart_max"] for r in good), max(r["flux_err_max"] for r in good),
             max(r["max_DC_scaled"] for r in good)))
    print("Pointwise-B0 series (%d rows): max rel div L2 = %.1e, max max|D_h b| = %.1e, max flux err = %.1e"
          % (len(bad), max(r["div_rel_L2"] for r in bad), max(r["div_max"] for r in bad), max(r["flux_err_max"] for r in bad)))

if __name__ == "__main__":
    {"exp1": exp1, "exp2": exp2}[sys.argv[1]]()
