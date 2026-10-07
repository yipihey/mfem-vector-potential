#!/usr/bin/env python3
"""E7 analysis: figures + tables for results/exp7_summary.md.
Usage: scripts/plot_exp7.py [perf.csv] [baseline.csv] [outdir]   (tables -> <outdir>/exp7_tables.md)"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

perf = sys.argv[1] if len(sys.argv) > 1 else "results/exp7_perf.csv"
base = sys.argv[2] if len(sys.argv) > 2 else "results/exp7_baseline.csv"
out = sys.argv[3] if len(sys.argv) > 3 else "results"

df = pd.read_csv(perf)
bdf = pd.read_csv(base) if os.path.exists(base) else None
key = ["opname", "N", "p", "np"]
# best repetition = the row (tag) with the smallest one-shot remap time; breakdown taken from that row
best = df.loc[df.groupby(key)["t_remap_oneshot"].idxmin()].copy()
best["nrep"] = df.groupby(key)["t_remap_oneshot"].transform("size").loc[best.index]
# min over repetitions of the two headline numbers independently
best["t_steady_min"] = df.groupby(key)["t_remap_steady"].transform("min").loc[best.index]
best["t_oneshot_min"] = df.groupby(key)["t_remap_oneshot"].transform("min").loc[best.index]
best["t_oneshot_max"] = df.groupby(key)["t_remap_oneshot"].transform("max").loc[best.index]
best["t_curl_min"] = df.groupby(key)["t_curl"].transform("min").loc[best.index]
best["t_build_min"] = df.groupby(key)["t_build"].transform("min").loc[best.index]
if bdf is not None:
    bb = bdf.groupby(["N", "p", "np"]).agg(base_peak=("peak_build_mb", "min"), base_peak_sum=("peak_build_sum_mb", "min"),
                                          base_build=("t_build", "min")).reset_index()
    best = best.merge(bb, on=["N", "p", "np"], how="left")
best["t_findpts_grp"] = best["t_ev_setup"] + best["t_findpoints"]
best["t_eval_grp"] = best["t_collect"] + best["t_eval"]
best["t_func_grp"] = best["t_dofs"] + best["t_plan_other"]
best["t_solve_grp"] = best["t_solve"]
best["t_gauge_grp"] = best["t_gauge"]
best["t_comp_grp"] = best["t_compose"]
best["t_sum_grp"] = best[["t_findpts_grp", "t_eval_grp", "t_func_grp", "t_solve_grp", "t_gauge_grp", "t_comp_grp"]].sum(axis=1)

OPS = ["A_pt", "A_int", "A_int+coulomb", "A_int+jacobi:20", "A_l2+coulomb", "B_l2", "B_l2c"]
OPS = [o for o in OPS if o in set(best.opname)]
col = dict(zip(["A_pt", "A_int", "A_int+coulomb", "A_int+jacobi:20", "A_l2+coulomb", "B_l2", "B_l2c"],
               ["#8c8c8c", "#1f77b4", "#2ca02c", "#17becf", "#9467bd", "#ff7f0e", "#d62728"]))
mk = {1: "o", 2: "s", 3: "^", 4: "D"}
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
f1 = lambda x: f"{x:.3g}"
L = []   # markdown lines


def slope(x, y):
    x = np.log(np.asarray(x, float)); y = np.log(np.asarray(y, float))
    if len(x) < 2:
        return np.nan
    return np.polyfit(x, y, 1)[0]


def sel(op, p=None, npr=1):
    s = best[(best.opname == op) & (best.np == npr)]
    if p is not None:
        s = s[s.p == p]
    return s.sort_values("nd")

# ---------------------------------------------------------------- 1. time vs ND with slopes
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
for ax, (cname, title) in zip(axs, [("t_oneshot_min", "one-shot remap (evaluator setup + FindPoints + plan + apply + gauge + compose)"),
                                    ("t_steady_min", "steady remap (plan reused: eval + dofs/solves + gauge + compose)"),
                                    ("t_build_min", "mesh + spaces + operators build (M1 and M2)")]):
    for o in OPS if cname != "t_build_min" else ["A_int"]:
        for p in (1, 2, 3, 4):
            s = sel(o, p)
            if len(s):
                ax.loglog(s.nd, s[cname], marker=mk[p], color=col[o], ls="-" if p == 2 else ":", ms=5, lw=1,
                          label=f"{o}" if p == 2 else None)
    ax.set_xlabel("ND true dofs"); ax.set_ylabel("time [s] (1 rank)"); ax.set_title(title, fontsize=8)
    ref = np.array([1e3, 1e6])
    ax.loglog(ref, ref * (best[cname].min() / 1e3) * 0.5, "k--", lw=0.6, alpha=0.5)
axs[0].legend(fontsize=7)
fig.text(0.5, 0.005, "markers: o p=1, s p=2, ^ p=3, D p=4 (solid p=2, dotted others); dashed black: slope 1", ha="center", fontsize=8)
fig.tight_layout(); fig.savefig(f"{out}/exp7_time_vs_nd.png", dpi=130); plt.close(fig)

L.append("### Table 7.1 Time per remap (s, 1 rank, min of repetitions) and log-log slope vs ND\n")
L.append("`oneshot` = RemoteEvaluator setup + collect + FindPoints + plan setup + apply (min of the timed repetitions, after one untimed warm-up apply that absorbs lazy AMG setup) + gauge fix + compose (mesh/space/operator build excluded); `steady` = plan reused (the E4 convention). Slope = d log t / d log ND from a least-squares fit over the N values at fixed p (blank if <3 N values, `[2pt]` = two points).\n")
Ns = sorted(best[best.np == 1].groupby(["p", "N"]).size().index.tolist())
hdr = "| op | p | " + " | ".join(f"N={n}" for n in sorted(set(n for _, n in Ns))) + " | slope oneshot | slope steady |"
L.append(hdr); L.append("|" + "---|" * (hdr.count("|") - 1))
allN = sorted(set(n for _, n in Ns))
slopes = {}
for o in OPS:
    for p in (1, 2, 3, 4):
        s = sel(o, p)
        if not len(s):
            continue
        cells = []
        for n in allN:
            r = s[s.N == n]
            cells.append(f"{r.t_oneshot_min.iloc[0]:.3g} / {r.t_steady_min.iloc[0]:.3g}" if len(r) else "")
        so = slope(s.nd, s.t_oneshot_min); ss = slope(s.nd, s.t_steady_min)
        tag = "" if len(s) >= 3 else " [2pt]"
        slopes[(o, p)] = (so, ss, len(s))
        L.append(f"| {o} | {p} | " + " | ".join(cells) + f" | {so:.2f}{tag} | {ss:.2f}{tag} |" if len(s) >= 2 else
                 f"| {o} | {p} | " + " | ".join(cells) + " | | |")
L.append("\n(each cell: oneshot / steady, seconds)\n")

# ---------------------------------------------------------------- 2. fractions
L.append("### Table 7.2 Where the one-shot remap time goes (1 rank; % of the sum of the parts)\n")
L.append("FindPoints = RemoteEvaluator setup + GSLIB search (incl. periodic-image retries); eval = collecting the target points + interpolating the source at them; functionals = dof functionals / rhs assembly + plan setup (mass operator, AMG for B_l2c); solves = CG solves (L2 mass, div-clean Schur); gauge = Coulomb/Jacobi gauge fix; compose = b = b0 + C a with b0 precomputed (the one-off projection of the constant B0 is excluded from all totals; see Table 7.4).\n")
L.append("| op | p | N | ND | total [s] | FindPoints | eval | functionals | solves | gauge | compose |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|")
frac_cases = [(2, 16), (3, 16), (2, 24), (2, 32), (4, 8), (1, 32)]
frac_cases = [(p, n) for p, n in frac_cases if len(best[(best.p == p) & (best.N == n) & (best.np == 1)])]
if not frac_cases:   # (smoke data) largest N of each p
    frac_cases = [(p, int(best[(best.p == p) & (best.np == 1)].N.max())) for p in sorted(set(best.p))]
fig, axs = plt.subplots(1, len(frac_cases), figsize=(3.4 * len(frac_cases), 4.2), squeeze=False)
grp = [("t_findpts_grp", "FindPoints", "#4c72b0"), ("t_eval_grp", "eval", "#8fb4e3"), ("t_func_grp", "functionals", "#55a868"),
       ("t_solve_grp", "solves", "#c44e52"), ("t_gauge_grp", "gauge", "#8172b2"), ("t_comp_grp", "compose", "#ccb974")]
for ax, (p, n) in zip(axs[0], frac_cases):
    ops_here = [o for o in OPS if len(sel(o, p)[sel(o, p).N == n])]
    bottoms = np.zeros(len(ops_here))
    for g, lab, c in grp:
        vals = np.array([sel(o, p)[sel(o, p).N == n][g].iloc[0] / sel(o, p)[sel(o, p).N == n].t_sum_grp.iloc[0] for o in ops_here])
        ax.bar(range(len(ops_here)), vals, bottom=bottoms, color=c, label=lab)
        bottoms += vals
    ax.set_xticks(range(len(ops_here))); ax.set_xticklabels(ops_here, rotation=60, ha="right", fontsize=7)
    ax.set_title(f"p={p} N={n}"); ax.set_ylim(0, 1)
    for o in ops_here:
        r = sel(o, p)[sel(o, p).N == n].iloc[0]
        t = r.t_sum_grp
        L.append(f"| {o} | {p} | {n} | {int(r.nd)} | {t:.3g} | " + " | ".join(f"{100 * r[g] / t:.0f}%" for g, _, _ in grp) + " |")
axs[0][0].set_ylabel("fraction of remap time"); axs[0][-1].legend(fontsize=7, loc="upper right")
fig.tight_layout(); fig.savefig(f"{out}/exp7_fractions.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 3. curl apply and memory
fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
cs = best[(best.opname == "A_int") & (best.np == 1)]
for p in (1, 2, 3, 4):
    s = cs[cs.p == p].sort_values("nd")
    if len(s):
        axs[0].loglog(s.nd, s.t_curl_min / s.nd * 1e9, marker=mk[p], label=f"p={p}")
        axs[1].loglog(s.nd, s.t_curl_min, marker=mk[p], label=f"p={p}")
axs[0].set_ylabel("curl apply [ns per ND dof]"); axs[0].set_xlabel("ND dofs"); axs[0].legend()
axs[1].set_ylabel("curl apply C_h a [s]"); axs[1].set_xlabel("ND dofs")
for o in OPS:
    for p in (1, 2, 3, 4):
        s = sel(o, p)
        if len(s):
            axs[2].semilogx(s.nd, s.peak_sum_mb * 1024 / s.nd, marker=mk[p], color=col[o], ls="-" if p == 2 else ":", ms=5, lw=1,
                            label=o if p == 2 else None)
axs[2].set_xlabel("ND dofs"); axs[2].set_ylabel("peak RSS (1 rank) [kB per ND dof]"); axs[2].legend(fontsize=7)
axs[2].set_ylim(bottom=0)
fig.tight_layout(); fig.savefig(f"{out}/exp7_curl_memory.png", dpi=130); plt.close(fig)

L.append("\n### Table 7.3 Operators, DOFs, memory (1 rank; peak RSS = VmHWM of the whole process, includes the two meshes, both operator sets, the exact source projection and GSLIB)\n")
L.append("nnz/dof: C_h (ND->RT) per ND dof, D_h per RT dof (the discrete operators are H(curl)/H(div) incidence-type matrices with small constant rows; at p=1 they are exactly the signed incidence matrices). 'increment' = peak(op) - peak(baseline: meshes+spaces+operators only).\n")
L.append("| p | N | NE | ND | RT | L2 | H1 | nnz C/ND | nnz D/RT | baseline peak [MB] | A_int peak [MB] | B_l2 peak [MB] | B_l2c peak [MB] | A_int+coulomb peak [MB] | B_l2c kB/ND (increment) |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for (p, n), _ in best[best.np == 1].groupby(["p", "N"]):
    g = best[(best.p == p) & (best.N == n) & (best.np == 1)]
    r0 = g.iloc[0]
    def pk(o):
        r = g[g.opname == o]
        return f1(r.peak_mb.iloc[0]) if len(r) else ""
    bp = f1(r0.base_peak) if "base_peak" in g and not np.isnan(r0.get("base_peak", np.nan)) else ""
    inc = ""
    r = g[g.opname == "B_l2c"]
    if len(r) and "base_peak" in g:
        inc = f1((r.peak_mb.iloc[0] - r.base_peak.iloc[0]) * 1024 / r.nd.iloc[0])
    L.append(f"| {p} | {n} | {int(r0["ne"])} | {int(r0.nd)} | {int(r0.rt)} | {int(r0.l2)} | {int(r0.h1)} | {r0.nnz_C_per_nd:.2f} | {r0.nnz_D_per_rt:.2f} | {bp} | {pk('A_int')} | {pk('B_l2')} | {pk('B_l2c')} | {pk('A_int+coulomb')} | {inc} |")

L.append("\n### Table 7.4 Curl apply C_h a (mean of 10, 1 rank) and composition b = b0 + C a\n")
L.append("| p | N | ND | t_curl [ms] | ns per ND dof | t_compose [ms] | b0 projection (one-off per mesh) [s] | mesh+spaces+operators build [s] |")
L.append("|---|---|---|---|---|---|---|---|")
for (p, n), _ in best[best.np == 1].groupby(["p", "N"]):
    g = best[(best.p == p) & (best.N == n) & (best.np == 1) & (best.opname == "A_int")]
    if len(g):
        r = g.iloc[0]
        L.append(f"| {p} | {n} | {int(r.nd)} | {1e3 * r.t_curl_min:.3g} | {1e9 * r.t_curl_min / r.nd:.3g} | {1e3 * r.t_compose:.3g} | {r.t_b0proj:.3g} | {r.t_build_min:.3g} |")

# ---------------------------------------------------------------- 4. iteration counts
L.append("\n### Table 7.5 Iteration counts (1 rank)\n")
L.append("L2 = CG iterations of the (Jacobi-PCG, PA) mass solve, rtol 1e-12; clean = outer (AMG-preconditioned CG on the Schur complement, summed over defect-correction rounds) / inner (M_RT CG, rtol 1e-14) / rounds; gauge = AMG-CG iterations (coulomb) or Chebyshev-Jacobi sweeps.\n")
L.append("| p | N | A_l2 (ND mass CG) | B_l2 (RT mass CG) | B_l2c outer / inner / rounds | A_int+coulomb | A_l2+coulomb | A_int+jacobi:20 |")
L.append("|---|---|---|---|---|---|---|---|")
for (p, n), _ in best[best.np == 1].groupby(["p", "N"]):
    g = best[(best.p == p) & (best.N == n) & (best.np == 1)]
    def one(o, c):
        r = g[g.opname == o]
        return str(int(r[c].iloc[0])) if len(r) else ""
    c3 = g[g.opname == "B_l2c"]
    L.append(f"| {p} | {n} | {one('A_l2+coulomb', 'iters_l2')} | {one('B_l2', 'iters_l2')} | "
             + (f"{int(c3.iters.iloc[0])} / {int(c3.inner_iters.iloc[0])} / {int(c3.clean_rounds.iloc[0])}" if len(c3) else "")
             + f" | {one('A_int+coulomb', 'gauge_iters')} | {one('A_l2+coulomb', 'gauge_iters')} | {one('A_int+jacobi:20', 'gauge_iters')} |")

# ---------------------------------------------------------------- 5. strong scaling
sc_cases = [(2, 16), (3, 16), (2, 24)]
sc_cases = [(p, n) for p, n in sc_cases if len(best[(best.p == p) & (best.N == n) & (best.np > 1)])]
if not sc_cases:
    sc_cases = sorted(set((int(r.p), int(r.N)) for r in best[best.np > 1].itertuples()))
L.append("\n### Table 7.6 Strong scaling (4 cores, oversubscribed/contended by the harness: indicative only)\n")
L.append("Speedup = T(1 rank)/T(np ranks), of the one-shot remap / of the plan+FindPoints part / of the apply (steady) part / of the solves only. Max-over-ranks times.\n")
L.append("| p | N | op | T1 oneshot [s] | S(2) oneshot | S(4) oneshot | S(4) FindPoints | S(4) eval+functionals | S(4) solves | S(4) gauge | S(4) curl |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|")
fig, axs = plt.subplots(1, max(1, len(sc_cases)), figsize=(4.6 * max(1, len(sc_cases)), 4.2), squeeze=False)
for ax, (p, n) in zip(axs[0], sc_cases):
    for o in OPS:
        s = best[(best.opname == o) & (best.p == p) & (best.N == n)].sort_values("np")
        if len(s) < 2 or s.np.iloc[0] != 1:
            continue
        ax.plot(s.np, s.t_oneshot_min.iloc[0] / s.t_oneshot_min, marker="o", color=col[o], label=o)
        t = lambda c, k: s[s.np == k][c].iloc[0] if len(s[s.np == k]) else np.nan
        sp = lambda c, k: t(c, 1) / t(c, k) if t(c, k) > 0 else np.nan
        def spf(c, k):
            v = sp(c, k)
            return "" if np.isnan(v) else f"{v:.2f}"
        s = s.copy()
        L.append(f"| {p} | {n} | {o} | {t('t_oneshot_min', 1):.3g} | {spf('t_oneshot_min', 2)} | {spf('t_oneshot_min', 4)} | {spf('t_findpts_grp', 4)} | "
                 f"{'' if np.isnan(t('t_eval_grp', 4)) else f'{(t('t_eval_grp', 1) + t('t_func_grp', 1)) / (t('t_eval_grp', 4) + t('t_func_grp', 4)):.2f}'} | {spf('t_solve_grp', 4)} | {spf('t_gauge_grp', 4)} | {spf('t_curl_min', 4)} |")
    ax.plot([1, 2, 4], [1, 2, 4], "k--", lw=0.6)
    ax.set_xscale("log", base=2); ax.set_yscale("log", base=2)
    ax.set_xticks([1, 2, 4]); ax.set_xticklabels([1, 2, 4]); ax.set_yticks([1, 2, 4]); ax.set_yticklabels([1, 2, 4])
    ax.set_title(f"p={p} N={n}: speedup of one-shot remap"); ax.set_xlabel("MPI ranks")
axs[0][0].legend(fontsize=7); axs[0][0].set_ylabel("speedup")
fig.tight_layout(); fig.savefig(f"{out}/exp7_strong_scaling.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 6. cost per dof vs p
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
L.append("\n### Table 7.7 Cost per ND dof vs p at fixed N (1 rank; microseconds per ND dof, one-shot / steady)\n")
pn = sorted(set(best[best.np == 1].N))
Nfix = [n for n in pn if len(set(best[(best.N == n) & (best.np == 1)].p)) >= 3]
L.append("| N | op | " + " | ".join(f"p={p}" for p in (1, 2, 3, 4)) + " |")
L.append("|---|---|---|---|---|---|")
for ax, n in zip(axs, Nfix):
    for o in OPS:
        s = best[(best.opname == o) & (best.N == n) & (best.np == 1)].sort_values("p")
        if len(s):
            ax.semilogy(s.p, s.t_oneshot_min / s.nd * 1e6, marker="o", color=col[o], label=o)
            cells = []
            for p in (1, 2, 3, 4):
                r = s[s.p == p]
                cells.append(f"{1e6 * r.t_oneshot_min.iloc[0] / r.nd.iloc[0]:.3g} / {1e6 * r.t_steady_min.iloc[0] / r.nd.iloc[0]:.3g}" if len(r) else "")
            L.append(f"| {n} | {o} | " + " | ".join(cells) + " |")
    ax.set_title(f"N={n}"); ax.set_xlabel("p"); ax.set_ylabel("one-shot remap [us per ND dof]"); ax.set_xticks([1, 2, 3, 4])
if len(Nfix):
    axs[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{out}/exp7_cost_vs_p.png", dpi=130); plt.close(fig)

# ---------------------------------------------------------------- 7. FindPoints / misc diagnostics
L.append("\n### Table 7.8 Points searched, FindPoints cost, accuracy sanity (1 rank, A_int and B_l2 rows)\n")
L.append("| p | N | op | #points | image-shifted | not found | t_FindPoints [s] | us per point | t_ev_setup [s] | rel err B (vs exact) | max|D_h b| |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|")
for o in ("A_int", "B_l2", "B_l2c"):
    for r in best[(best.opname == o) & (best.np == 1)].sort_values(["p", "N"]).itertuples():
        L.append(f"| {r.p} | {r.N} | {o} | {int(r.npoints)} | {int(r.nshifted)} | {int(r.nnotfound)} | {r.t_findpoints:.3g} | {1e6 * r.t_findpoints / r.npoints:.3g} | {r.t_ev_setup:.3g} | {r.errB_rel:.3e} | {r.div_max:.1e} |")


# ---------------------------------------------------------------- 8. gauge object setup (per-mesh cost, NOT in the remap totals above)
L.append("\n### Table 7.10 Gauge-fix object: per-mesh setup vs per-remap application (1 rank, seconds)\n")
L.append("The CoulombGauge object (assembled ND mass, L = G^T M_ND G by sparse triple product, AMG) is built once per mesh; on a moving (ALE) mesh it would be rebuilt every remap. The setup is NOT included in the one-shot/steady totals of Table 7.1 (it is an implementation cost of this code: L equals the H1 stiffness matrix and could be assembled directly, the Jacobi variant could be matrix-free). 'first' = first application (lazy AMG setup), 'apply' = min of the later ones.\n")
L.append("| op | p | N | ND | setup | first apply | apply | iters | setup / (steady remap incl. gauge apply) |")
L.append("|---|---|---|---|---|---|---|---|---|")
for o in ("A_int+coulomb", "A_l2+coulomb", "A_int+jacobi:20"):
    for r in best[(best.opname == o) & (best.np == 1)].sort_values(["p", "N"]).itertuples():
        L.append(f"| {o} | {r.p} | {r.N} | {int(r.nd)} | {r.t_gauge_setup:.3g} | {r.t_gauge_first:.3g} | {r.t_gauge:.3g} | {int(r.gauge_iters)} | {r.t_gauge_setup / r.t_steady_min:.1f} |")

L.append("\n### Table 7.9 Slopes (log-log, 1 rank) per op and p\n")
L.append("| op | p | #N values | slope oneshot | slope steady |")
L.append("|---|---|---|---|---|")
for (o, p), (so, ss, k) in slopes.items():
    if k < 2:
        continue
    L.append(f"| {o} | {p} | {k} | {so:.2f} | {ss:.2f} |")
# global fit over all p>=2 data (cost per dof is roughly flat in p -> mixes p dependence into the slope; informative only)
open(f"{out}/exp7_tables.md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
