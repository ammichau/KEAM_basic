"""Compact set of moments used to compare scenarios (faithful vs corrected)."""
from __future__ import annotations
import numpy as np
from . import stats
from .simulate import SimResult
from .solve import Solution


def moments(r: SimResult, sol: Solution | None = None, nilfdef_alt: float | None = None) -> dict:
    cs = stats.cross_section(r)
    ca = stats.careers(r)
    cy_all = stats.cycle(r, expansion_includes_recessions=True)
    cy = stats.cycle(r, expansion_includes_recessions=False)
    alive = r.Alive[:, : r.Tsim] == 1
    emp = (r.EmpI[:, : r.Tsim] == 1) & alive
    e_over = float(np.mean(r.exp[:, : r.Tsim][alive] > r.params.egrid[-1]))
    m = {
        "E/pop (all ages)": float(r.EmpI[:, : r.Tsim][alive].mean()),
        "E/pop young": cs["Erate"][0], "E/pop middle": cs["Erate"][1], "E/pop old": cs["Erate"][2],
        "U/pop (search>NiLFbar)": float(r.UnempI[:, : r.Tsim][alive].mean()),
        "hours | E": float(r.Hours[:, : r.Tsim][emp].mean()),
        "quit rate/qtr (exp)": cy["Quit"]["exp"], "quit rate/qtr (rec)": cy["Quit"]["rec"],
        "E->nonE/qtr (exp)": cy["EnonE"]["exp"], "E->nonE/qtr (rec)": cy["EnonE"]["rec"],
        "dE/pop rec-exp (pts)": 100 * (cy["Emp"]["rec"] - cy["Emp"]["exp"]),
        "dE/pop rec-all (pts, MATLAB)": 100 * (cy_all["Emp"]["rec"] - cy_all["Emp"]["exp"]),
        "dHours rec-exp (%)": 100 * (cy["Hours"]["rec"] / cy["Hours"]["exp"] - 1),
        "wife inc share (exp)": cy["Wshare"]["exp"], "wife inc share (rec)": cy["Wshare"]["rec"],
        "share PT": ca["PT"], "share Lifecycle": ca["Cycle"], "share Career": ca["Career"],
        "share NiLF": ca["NiLF"], "share undefined (of all)": ca["Undefined"],
        "share person-qtrs e > grid max": e_over,
    }
    if sol is not None:
        m["mean gQ (working ages)"] = float(sol.gQ[:, :3].mean())
        m["mean gH (working ages)"] = float(sol.gH[:, :3].mean())
        m["mean gS (working ages)"] = float(sol.gS[:, :3].mean())
        m["VFI iters type0 (y,m,o)"] = str([int(x) for x in sol.n_iter[0, :3]])
    return m
