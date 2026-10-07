# mfem-vector-potential

Numerical study: is transferring the magnetic vector potential A (H(curl))
across ALE rezoning, with B = B0 + curl_h A in H(div), better than
transferring B directly with divergence cleaning? Built on MFEM; compared
against MHD-ALE.

Start with `docs/results.md` (findings and recommendation), then
`docs/design_note.md` (discrete representation) and `PLAN.md`.

| what | where |
|---|---|
| build provenance (commits, compilers, commands) | `BUILD.md`, `scripts/build_all.sh` |
| core library: periodic meshes, deformations, de Rham operators, fields, diagnostics | `src/vp_core.{hpp,cpp}` |
| mesh-to-mesh transfers (A_pt, A_int, A_l2, B_pt, B_int, B_l2, B_l2c) | `src/vp_transfer.{hpp,cpp}` |
| Coulomb gauge projection and local gauge smoothing | `src/vp_gauge.{hpp,cpp}` |
| Experiment drivers (static, distorted, one remap, repeated, ALE cycle, performance) | `experiments/exp{1,2,3,4,5,7}_*.cpp`, `experiments/selftest.cpp` |
| run scripts, plots, tables | `scripts/` |
| CSV data, figures, per-experiment summaries | `results/` |
| MHD-ALE instrumentation patch, run scripts, analysis | `mhdale/`, `results/mhdale_summary.md` |

Build: `scripts/build_all.sh` (external dependencies into `external/`), then
`make`. Self-test: `scripts/selftest.sh`. Each experiment has
`scripts/run_expN.sh` and `scripts/plot_expN.py`.
