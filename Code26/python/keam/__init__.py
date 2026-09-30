"""Python translation of the KEAM 'simpler model' MATLAB code (Code26).

Modules
-------
params     : parameter container reproducing SimplerMod_May17_splines.m lines 28-113
functions  : utility, wage, home production and experience technology
solve      : value function iteration (SimplerMod_May17_splines.m lines 133-393)
simulate   : life-cycle panel simulation (SimplerMod_May17_sim.m lines 48-554)
stats      : aggregate statistics (SimplerMod_May17_sim.m lines 559-1118)
matlab_rng : exact replication of MATLAB rng(seed,'twister') streams
matio      : loaders for the saved MATLAB solution files

Every known departure of the MATLAB computation from the model written in
the paper/slides is controlled by a flag in `params.Options`.  With
`Options.faithful()` (the default) the Python code reproduces the MATLAB
output; with `Options.corrected()` all flags are switched to the textbook
behaviour.  See DEPARTURES.md at the repository root for the catalogue.
"""
from .params import Params, Options
from .solve import solve, Solution
from .simulate import simulate, SimResult, SimConfig
from . import stats

__all__ = ["Params", "Options", "solve", "Solution", "simulate", "SimResult", "SimConfig", "stats"]
