#!/usr/bin/env python3
"""Analysis of the E6 MHD-ALE runs.  Usage: python3 analyze.py [outroot] [tag]
Reads <outroot>/<run>/{step_diag.dat,remap_diag.dat,MHD.out,walltime.txt}, writes
results/mhdale_<tag>_*.png and results/mhdale_summary.md (section appended per tag)."""
import sys, os, re, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/home/user/mfem-vector-potential"
OUT = sys.argv[1] if len(sys.argv) > 1 else REPO + "/external/MHD-ALE/output/E6"
TAG = sys.argv[2] if len(sys.argv) > 2 else "3d"
RES = REPO + "/results"
MU = 4.0 * np.pi   # mu0 of testcase_TaylorGreen.cpp; E_B = 0.5 int B^2 / mu


def load_table(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        hdr = f.readline().lstrip("# ").split()
    d = np.loadtxt(path, skiprows=1, ndmin=2)
    if d.size == 0:
        return {k: np.zeros(0) for k in hdr}
    return {k: d[:, i] for i, k in enumerate(hdr)}


def final_errors(path):
    """errors printed by the code at the end of the run (last block of 'error of' lines)"""
    lines = [l for l in open(path) if "error of" in l]
    last = lines[-12:]
    out = {}
    for l in last:
        m = re.match(r"L_(inf|1|2)\s+error of (\w+):\s*(\S+)", l)
        if m:
            out[(m.group(2), m.group(1))] = float(m.group(3))
    return out


def parse_name(n):
    m = re.match(r"(noale|fsri(\d+)_rma(\d))", n)
    if n.startswith("noale"):
        return "no ALE (Lagrangian)"
    m = re.match(r"fsri(\d+)_rma(\d)", n)
    return "remap/%s steps, -rma %s" % (m.group(1), m.group(2))


runs = {}
for d in sorted(glob.glob(OUT + "/*/")):
    n = os.path.basename(d.rstrip("/"))
    if not os.path.exists(d + "step_diag.dat") or not os.path.exists(d + "MHD.out"):
        continue
    r = dict(name=n, label=parse_name(n))
    r["step"] = load_table(d + "step_diag.dat")
    r["remap"] = load_table(d + "remap_diag.dat")
    r["err"] = final_errors(d + "MHD.out")
    r["has_final"] = len(r["err"]) == 12
    wt = d + "walltime.txt"
    r["wall"] = float(open(wt).read().split()[0]) if os.path.exists(wt) and open(wt).read().split() else np.nan
    runs[n] = r
if not runs:
    sys.exit("no runs found in " + OUT)

# group by (tf, rs) suffix so that plots compare like with like
groups = {}
for n, r in runs.items():
    key = re.search(r"rs\d_tf[\d.]+", n).group(0)
    groups.setdefault(key, []).append(n)

colors = {"noale": "k", "fsri10_rma0": "tab:blue", "fsri10_rma2": "tab:orange",
          "fsri2_rma0": "tab:green", "fsri2_rma2": "tab:red",
          "fsri5_rma0": "tab:cyan", "fsri5_rma2": "tab:brown"}


def col(n):
    for k, c in colors.items():
        if n.startswith(k):
            return c
    return None


def sty(n):
    return "--" if "rma2" in n else "-"


def new_axes(nrow, ncol, w=5.2, h=3.6):
    fig, ax = plt.subplots(nrow, ncol, figsize=(w * ncol, h * nrow), squeeze=False)
    return fig, ax


summary_rows = []
for gk, names in groups.items():
    nm = [n for n in names]
    # ---- energies, errors
    fig, ax = new_axes(2, 2)
    for n in nm:
        s = runs[n]["step"]
        EB = s["EB"] / MU
        ax[0, 0].plot(s["t"], EB / EB[0] - 1, sty(n), color=col(n), label=runs[n]["label"])
        ax[0, 1].plot(s["t"], s["KE"] / s["KE"][0] - 1, sty(n), color=col(n))
        tot = EB + s["KE"] + s["IE"]
        ax[1, 0].plot(s["t"], s["errB_L2"], sty(n), color=col(n))
        ax[1, 1].plot(s["t"], (s["IE"]) / s["IE"][0] - 1, sty(n), color=col(n))
    ax[0, 0].set_title("magnetic energy E_B = 0.5 int B^2 / mu: relative change")
    ax[0, 1].set_title("kinetic energy: relative change")
    ax[1, 0].set_title("L2 error of B vs exact (steady) B")
    ax[1, 1].set_title("internal energy: relative change")
    for a in ax.ravel():
        a.set_xlabel("t"); a.grid(alpha=.3)
    ax[0, 0].legend(fontsize=7)
    fig.suptitle(TAG.upper() + " Taylor-Green, " + gk); fig.tight_layout()
    fig.savefig("%s/mhdale_%s_energy_%s.png" % (RES, TAG, gk), dpi=110); plt.close(fig)

    # ---- divergence and helicity
    fig, ax = new_axes(1, 3)
    for n in nm:
        s = runs[n]["step"]
        ax[0, 0].semilogy(s["t"], np.maximum(s["divL2"], 1e-18), sty(n), color=col(n), label=runs[n]["label"])
        ax[0, 1].semilogy(s["t"], np.maximum(s["divMax"], 1e-18), sty(n), color=col(n))
        ax[0, 2].plot(s["t"], s["helicity"], sty(n), color=col(n))
        r = runs[n]["remap"]
        if r is not None and len(r["t"]):
            ax[0, 0].semilogy(r["t"], np.maximum(r["divL2_s"], 1e-18), "x", color=col(n), ms=6)
            ax[0, 1].semilogy(r["t"], np.maximum(r["divMax_s"], 1e-18), "x", color=col(n), ms=6)
    ax[0, 0].set_title("L2 norm of div B (lines: A-route; x: shadow B)")
    ax[0, 1].set_title("max |dof| of div B (x: shadow)")
    ax[0, 2].set_title("helicity int A.B (A_z e_z, B in-plane: =0)")
    for a in ax.ravel():
        a.set_xlabel("t"); a.grid(alpha=.3)
    ax[0, 0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig("%s/mhdale_%s_div_helicity_%s.png" % (RES, TAG, gk), dpi=110); plt.close(fig)

    # ---- per-remap quantities
    rn = [n for n in nm if runs[n]["remap"] is not None and len(runs[n]["remap"]["t"])]
    if rn:
        fig, ax = new_axes(2, 3)
        for n in rn:
            r = runs[n]["remap"]
            jA = (r["EB_a"] - r["EB_b"]) / r["EB_b"]
            jS = (r["EB_s"] - r["EB_b"]) / r["EB_b"]
            ax[0, 0].plot(r["step"], jA, "o" + sty(n), color=col(n), label=runs[n]["label"] + " (A-route)")
            ax[0, 0].plot(r["step"], jS, "x:", color=col(n), label=runs[n]["label"] + " (shadow B)")
            ax[0, 1].semilogy(r["step"], np.maximum(r["divL2_a"], 1e-18), "o" + sty(n), color=col(n), label="A-route")
            ax[0, 1].semilogy(r["step"], np.maximum(r["divL2_s"], 1e-18), "x:", color=col(n), label="shadow")
            ax[0, 2].semilogy(r["step"], r["errB_a"], "o" + sty(n), color=col(n))
            ax[0, 2].semilogy(r["step"], r["errB_s"], "x:", color=col(n))
            ax[1, 0].semilogy(r["step"], np.maximum(r["L2diff_sA"], 1e-18), "o" + sty(n), color=col(n))
            ax[1, 1].plot(r["step"], r["wall_remap_total"], "o" + sty(n), color=col(n), label=runs[n]["label"] + " remap")
            s = runs[n]["step"]
            ax[1, 1].plot(s["step"][1:], s["wall_ode"][1:], ".", color=col(n), alpha=.4)
            ax[1, 2].plot(r["step"], r["maxdisp"] / r["hmin"], "o" + sty(n), color=col(n))
        ax[0, 0].set_title("per-remap E_B jump (E_after-E_before)/E_before")
        ax[0, 0].legend(fontsize=6)
        ax[0, 1].set_title("L2 norm div B after remap (o: A, x: shadow)")
        ax[0, 2].set_title("L2 error of B vs exact after remap (o: A, x: shadow)")
        ax[1, 0].set_title("|| B_shadow - B_A-route ||_L2 on new mesh")
        ax[1, 1].set_title("wall time: remap (markers) vs ODE step (dots) [s]")
        ax[1, 1].set_yscale("log"); ax[1, 1].legend(fontsize=6)
        ax[1, 2].set_title("max |newnodes - nodes| / h_min at remap")
        for a in ax.ravel():
            a.set_xlabel("step"); a.grid(alpha=.3)
        fig.suptitle(TAG.upper() + " Taylor-Green per-remap diagnostics, " + gk); fig.tight_layout()
        fig.savefig("%s/mhdale_%s_remap_%s.png" % (RES, TAG, gk), dpi=110); plt.close(fig)

# ---------------------------------------------------------------- summary table
def fmt(x, f="%.3e"):
    return f % x if np.isfinite(x) else "n/a"

lines = []
lines.append("### %s runs (output root `%s`)\n" % (TAG, os.path.relpath(OUT, REPO)))
lines.append("E_B = 0.5 int B^2/mu (mu = 4 pi). dE_B total = (E_B(end)-E_B(0))/E_B(0). 'remap sum' = sum over remaps of (E_after-E_before)/E_B(0) for the A-route "
             "(what the code does) and 'shadow sum' the same for the direct-B (RT GSLIB interpolation) shadow transfer, which was NOT used for evolution. "
             "max divB columns: L2 norm of div B over the whole run for the A-route (every step, from step_diag.dat), and the max over remap events of the L2 norm / max|dof| of div B of the shadow B. "
             "'remap overhead' = sum of (A + rho + v + e + mesh) remap wall time / (sum ODE-step + remap wall time), instrumentation excluded; 'A-remap' = only the A-remap part.\n")
hdr = ["run", "steps", "#remaps", "final L2 err B", "final L2 err v", "final L2 err e", "dE_B total", "remap sum (A)", "shadow sum", "dKE total",
       "max abs helicity", "max divB L2 (A-route)", "max divB L2 (shadow)", "max abs div dof (A / shadow)", "mean L2diff(shadow,A)/normB",
       "mean abs E jump A / shadow", "mean rel. B-err change at remap A / shadow", "remap overhead", "A-remap frac", "shadow cost / A-remap", "total wall [s]"]
lines.append("| " + " | ".join(hdr) + " |")
lines.append("|" + "|".join(["---"] * len(hdr)) + "|")
for gk, names in sorted(groups.items()):
    for n in sorted(names):
        r = runs[n]; s = r["step"]; rm = r["remap"]
        nrem = 0 if rm is None else len(rm["t"])
        EB = s["EB"] / MU
        dEB = EB[-1] / EB[0] - 1
        dKE = s["KE"][-1] / s["KE"][0] - 1
        if nrem:
            sumA = np.sum((rm["EB_a"] - rm["EB_b"]) / MU) / EB[0]
            sumS = np.sum((rm["EB_s"] - rm["EB_b"]) / MU) / EB[0]
            mdS = np.max(rm["divL2_s"]); mdSm = np.max(rm["divMax_s"]); mdAm = np.max(rm["divMax_a"])
            l2d = np.mean(rm["L2diff_sA"] / rm["normB_a"])
            wr = np.sum(rm["wall_remap_total"]); wA = np.sum(rm["wall_Aremap"]); wS = np.sum(rm["wall_shadow"])
        else:
            sumA = sumS = mdS = mdSm = mdAm = l2d = np.nan; wr = wA = wS = 0.0
        wtot = np.sum(s["wall_ode"]) + wr
        if nrem:
            jA = np.mean(np.abs((rm["EB_a"] - rm["EB_b"]) / rm["EB_b"])); jS = np.mean(np.abs((rm["EB_s"] - rm["EB_b"]) / rm["EB_b"]))
            eA = np.mean((rm["errB_a"] - rm["errB_b"]) / rm["errB_b"]); eS = np.mean((rm["errB_s"] - rm["errB_b"]) / rm["errB_b"])
            jtxt = "%.1e / %.1e" % (jA, jS); etxt = "%+.2f / %+.2f" % (eA, eS)
        else:
            jtxt = etxt = "-"
        e = r["err"]
        lines.append("| " + " | ".join([
            n, "%d" % s["step"][-1], "%d" % nrem,
            fmt(e.get(("B", "2"), np.nan)), fmt(e.get(("v", "2"), np.nan)), fmt(e.get(("e", "2"), np.nan)),
            fmt(dEB), fmt(sumA), fmt(sumS), fmt(dKE),
            fmt(np.max(np.abs(s["helicity"]))), fmt(np.max(s["divL2"])), fmt(mdS),
            ("%.1e / %.1e" % (np.max(s["divMax"]), mdSm)) if nrem else "%.1e / -" % np.max(s["divMax"]),
            fmt(l2d), jtxt, etxt, fmt(wr / wtot, "%.3f") if nrem else "0", fmt(wA / wtot, "%.3f") if nrem else "0",
            fmt(wS / wA, "%.3f") if nrem and wA > 0 else "n/a", fmt(r["wall"], "%.0f")]) + " |")
lines.append("")
open("%s/mhdale_summary_%s.md" % (RES, TAG), "w").write("\n".join(lines))
print("\n".join(lines))
