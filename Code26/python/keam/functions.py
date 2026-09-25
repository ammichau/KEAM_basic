"""Utility, wage, home-production and experience functions.

Translations of utilC.m, utilL.m, du_dC.m, du_dL.m, wage.m, hprod.m, exp2.m,
de2_dh.m.  All functions are vectorised in their first argument(s).
"""
from __future__ import annotations
import numpy as np
from .params import Params


def utilC(p: Params, c):
    """phi_c * c^(1-crra)/(1-crra)."""
    return p.phi_c * np.power(c, 1.0 - p.crra) / (1.0 - p.crra)


def du_dC(p: Params, c):
    return p.phi_c * np.power(c, -p.crra)


def utilL(p: Params, h):
    """-mu * h^(1+eta)/(1+eta)."""
    return -p.mu * np.power(h, 1.0 + p.eta) / (1.0 + p.eta)


def du_dL(p: Params, h):
    return -p.mu * np.power(h, p.eta)


def wage(p: Params, ftt, e):
    """tau_wf * (ftt + gam_e * e^xi)   [wage.m]

    NOTE: additive in experience.  The slides write tau_w * omega * (1 + gam_e e^xi),
    which is multiplicative in the fixed type (DEPARTURES.md, D-P7).
    """
    return p.tau_wf_eff * (ftt + p.gam_e_eff * np.power(e, p.xi))


def hprod(p: Params, ftt):
    """ybar_h + z_h * ftt^alpha_h   [hprod.m]"""
    return p.ybar_h + p.z_h * np.power(ftt, p.alpha_h)


def exp2(p: Params, e, h):
    """e' = (1-delta_e) e + alpha_e e h^psi   [exp2.m]"""
    return (1.0 - p.delta_e) * e + p.alpha_e * e * np.power(h, p.psi)


def de2_dh(p: Params, e, h):
    """d e'/d h = alpha_e psi e h^(psi-1)   [de2_dh.m]"""
    return p.alpha_e * p.psi * e * np.power(h, p.psi - 1.0)
