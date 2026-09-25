"""Quantify the effect of each departure listed in DEPARTURES.md, one at a time.

For every solver flag the model is re-solved with that single flag switched to
the textbook behaviour (everything else faithful) and re-simulated with the
faithful simulator.  For every simulator/statistics flag the faithful solution
is re-simulated with that single flag switched.  Finally everything is
switched.  Output: output/departure_impacts.csv and a printed table.
"""
import sys, os, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd
from dataclasses import fields
from keam import Params, Options, solve, simulate, SimConfig
from keam.moments import moments
warnings.simplefilter("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
CODE26 = os.path.join(HERE, "..", "..")
p = Params.from_paras_mat(os.path.join(CODE26, "Solution", "Baseline", "paras.mat"))
cfg = SimConfig(NiLFdef=0.3)          # value used for the stored outputs / slides

solver_flags = ["aging_linear_index", "stale_h_nonemployed", "roundup_continuation",
                "derivative_at_current_state", "alpha_h_typo", "kapbar_indexing_bug",
                "young_cost_all_ages", "loose_vf_tolerance"]
sim_flags = ["reversed_e_weights", "reversed_type_weights_gQ", "same_seed_kappa_draws",
             "no_bcwage_in_sim", "all_employed_at_birth", "lagged_z_in_sim", "quit_before_loss",
             "expansion_includes_recessions", "drop_undefined_careers", "phantom_last_cohort"]
labels = {"aging_linear_index": "S1 ageing index", "stale_h_nonemployed": "S2 stale h (non-emp e')",
          "roundup_continuation": "S3 round-up e'", "derivative_at_current_state": "S4 dV/de at (y,z)",
          "alpha_h_typo": "S5 alpha_h typo", "kapbar_indexing_bug": "S6 kapbar indexing",
          "young_cost_all_ages": "S7 young cost all ages", "loose_vf_tolerance": "S8 VF tol 0.1",
          "reversed_e_weights": "M1 reversed e-weights", "reversed_type_weights_gQ": "M2 reversed gQ type-weights",
          "same_seed_kappa_draws": "M3 same seed kappa draws", "no_bcwage_in_sim": "M4 no BCwage in sim",
          "all_employed_at_birth": "M5 all employed at birth", "lagged_z_in_sim": "M6 lagged z (husband)",
          "quit_before_loss": "M7 quit before loss", "expansion_includes_recessions": "R1 ExpI=1 always",
          "drop_undefined_careers": "R2 drop undefined careers", "phantom_last_cohort": "R3 phantom last cohort"}

rows = {}
t0 = time.time()
base_opts = Options.faithful()
sol_f = solve(p, base_opts)
r = simulate(p, sol_f, base_opts, cfg)
rows["faithful (MATLAB)"] = moments(r, sol_f)
print(f"faithful done {time.time()-t0:.0f}s")

for f in solver_flags:
    o = base_opts.with_(**{f: False})
    sol = solve(p, o)
    r = simulate(p, sol, base_opts, cfg)       # simulator stays faithful
    rows[labels[f]] = moments(r, sol)
    dV = np.abs(sol.V[:, :3] - sol_f.V[:, :3]).max()
    rows[labels[f]]["max|dV| vs faithful"] = dV
    print(f"{labels[f]:28s} done {time.time()-t0:.0f}s  (max|dV|={dV:.3f}, iters type0={list(sol.n_iter[0,:3])})")

for f in sim_flags:
    o = base_opts.with_(**{f: False})
    r = simulate(p, sol_f, o, cfg)
    rows[labels[f]] = moments(r, sol_f)
    print(f"{labels[f]:28s} done {time.time()-t0:.0f}s")

o_solver_all = base_opts.with_(**{f: False for f in solver_flags})
sol_c = solve(p, o_solver_all)
r = simulate(p, sol_c, base_opts, cfg)
rows["all solver fixes, sim faithful"] = moments(r, sol_c)
rows["all solver fixes, sim faithful"]["max|dV| vs faithful"] = np.abs(sol_c.V[:, :3] - sol_f.V[:, :3]).max()
o_all = Options.corrected()
r = simulate(p, sol_c, o_all, cfg)
rows["all fixes"] = moments(r, sol_c)
print(f"all fixes done {time.time()-t0:.0f}s")

df = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_rows", 60)
out = os.path.join(HERE, "..", "output", "departure_impacts.csv")
df.to_csv(out)
print(df.T.to_string(float_format=lambda x: f"{x:.4f}"))
print("saved", out)
