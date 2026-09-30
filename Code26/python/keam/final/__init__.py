"""Final model: the KEAM life-cycle model as written in the slides, monthly, with assets.

Specification decisions (see FINAL_MODEL.md at the repository root):
  * monthly period, beta = 0.99
  * wage w = phi(Z) tau_w omega (1 + gam_e e^xi); experience e' = min(e_max, (1-delta) e + theta e h^psi)
  * home production f(omega) (1-h)^nu_h, f = ybar + z omega^alpha
  * fixed cost of work kappa_bar * kappa_m at ages 25-39, kappa_bar afterwards
  * husband: states E / R (re-employed with scar) / U, income shares (1, 0.85, 0.30)
  * assets a >= 0, gross return R = 1, retirement with a pension
  * every simulated woman's type is solved explicitly (no interpolation across types)
"""
from .params import FinalParams, make_types
from .solve import solve_all, FinalSolution
from .simulate import simulate_final, FinalSim, SimConfigFinal
from .moments import moments_final
