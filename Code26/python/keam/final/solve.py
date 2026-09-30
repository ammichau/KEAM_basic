"""Value function iteration for the final model (monthly, assets, per-type solve).

Per type the state is (age tau, experience e, assets a, husband y, aggregate z, cost-shock state j)
and there are two value functions: employed V^E (chooses hours h and savings a', may quit) and
non-employed V^N (chooses search s and savings a').  Choices are on grids; continuation values are
linearly interpolated in e', and in a' when the savings choice grid is finer than the state grid
(`nAc`; with a' restricted to the state grid, saving is suppressed when beta (1 + r) is close to
one because a' can only move in lumps of the grid spacing).  Modified policy iteration
(maximisation step followed by `howard_steps` evaluation steps) is used for speed.

The cost-of-work shock kappa_T is realised at the start of the period, before the quit decision:
an employed woman (or a job finder) with shock node j' gets max{V^E_j' - kappa_j', V^N_j'}.  With an
iid shock (rho_kT = 0) the values do not depend on the shock and a single shock state is carried
(last axis of length 1); with persistence the state carries the current node (last axis n_kT).
"""
from __future__ import annotations
from dataclasses import dataclass
import time
import numpy as np
from .params import FinalParams, make_types, make_types4


@dataclass
class FinalSolution:
    params: FinalParams
    omega: np.ndarray; kbar: np.ndarray; km: np.ndarray
    VE: np.ndarray; VN: np.ndarray          # (nK, nT, nE, nA, nY, nZ, nJ); nJ = 1 for an iid shock
    gH: np.ndarray; gAE: np.ndarray         # hours, savings when employed (values, not indices)
    gS: np.ndarray; gAN: np.ndarray         # search, savings when non-employed
    VR: np.ndarray                          # retirement value on the asset grid
    info: dict
    zh: np.ndarray = None                  # home-productivity multiplier per type (ones if absent)


def u(c, gamma):
    if abs(gamma - 1.0) < 1e-12:          # log utility (balanced-growth preferences)
        return np.log(c)
    return np.power(c, 1.0 - gamma) / (1.0 - gamma)


def U_kpr(x, gamma):
    """King-Plosser-Rebelo outer function of the composite x = log c - v(h) - kappa."""
    if abs(gamma - 1.0) < 1e-12:
        return x
    return np.exp((1.0 - gamma) * x) / (1.0 - gamma)


def _bracket(grid, x):
    k = np.clip(np.searchsorted(grid, x, side="right") - 1, 0, grid.size - 2)
    w = np.clip((x - grid[k]) / (grid[k + 1] - grid[k]), 0.0, 1.0)
    return k, w


def _interp_choice(grid, choice):
    """Brackets and weights to interpolate a function on `grid` at the choice points (identity when equal)."""
    if choice.size == grid.size and np.array_equal(choice, grid):
        return None, None
    return _bracket(grid, choice)


def solve_retirement(p: FinalParams):
    a = p.agrid; nA = a.size; ac = p.agrid_c
    yR = p.pension * p.yH_age[1]
    c = yR + (1.0 + p.r_a) * a[:, None] - ac[None, :]       # (a, a')
    flow = np.where(c > 1e-8, u(np.maximum(c, 1e-8), p.gamma), -1e10)
    V = flow.max(axis=1) / (1 - p.beta * (1 - p.death))
    kc, wc = _interp_choice(a, ac)
    # mortality acts as extra discounting (no certainty equivalent over death: with log utility the
    # level of V is arbitrary, so a "death value" of 0 would not be innocuous under risk sensitivity)
    for _ in range(5000):
        Vc = V if kc is None else (1 - wc) * V[kc] + wc * V[kc + 1]     # continuation on the choice grid
        Vn = (flow + p.beta * (1 - p.death) * Vc[None, :]).max(axis=1)
        if np.max(np.abs(Vn - V)) < 1e-9:
            V = Vn; break
        V = Vn
    return V


def solve_type(p: FinalParams, omega: float, kbar: float, km: float, VR: np.ndarray, verbose=False, zh: float = 1.0):
    nE, nA, nH, nS = p.nE, p.nA, p.nH, p.nS
    nY, nZ, nT = 3, 2, 3
    eg, ag, hg, sg = p.egrid, p.agrid, p.hgrid, p.sgrid
    ac = p.agrid_c; nAc = ac.size                                                   # savings choice grid
    kc, wc = _interp_choice(ag, ac)                                                 # None: choice on the state grid
    if kc is not None:
        wc1 = wc[None, None, None, :, None, None]

    def on_choice(X):
        """Continuation (j, e, h|s, a', y, z) on the state grid -> on the choice grid (axis 3)."""
        if kc is None:
            return X
        return (1 - wc1) * X[:, :, :, kc] + wc1 * X[:, :, :, kc + 1]
    beta = p.beta
    phi = np.array([1.0, p.phi_rec])
    w = phi[None, :] * p.tau_w * omega * (1 + p.gam_e * eg[:, None] ** p.xi)      # (nE, nZ)
    f = p.ybar_h + p.z_h * zh * omega ** p.alpha_h
    yH = p.y_husband()                                                              # (nT, nY, nZ)
    lamH = p.lamH()                                                                 # (nZ, nY, nY)
    lam_u = np.asarray(p.lam_u); lam_f = np.asarray(p.lam_f); lam_n = np.asarray(p.lam_n, float)
    pi_f = np.minimum(1.0, lam_n[None, :] + lam_f[None, :] * sg[:, None] ** p.nu)  # (nS, nZ): offers without search + search
    if p.lam_ue[0] > 0:                 # independent arrival rates: lam_ue in U (s >= s_bar), lam_n in N (s < s_bar)
        pi_f = np.where(sg[:, None] >= p.s_bar, np.asarray(p.lam_ue, float)[None, :], lam_n[None, :])
    p_age = p.p_age
    kT, nJ, PJ = p.kT_chain()                                                       # nodes, states, (nJ, n_kT)
    st = np.arange(kT.size) if nJ > 1 else np.zeros(kT.size, int)                   # state reached at node j'
    if p.kT_mult and not p.kpr:
        raise ValueError("kT_mult requires kpr preferences")
    st_of_j = np.arange(nJ)                                                         # node of state j (nJ = n_kT)
    kTdec = np.zeros_like(kT) if p.kT_mult else kT                                  # additive shock at the decision
    # iid proportional shock: V^E carries today's node (the flow depends on it) but the continuation and
    # V^N do not, so they are computed once (nJc = 1) and V^N is broadcast over the nodes when stored
    iidm = p.kT_mult and nJ > 1 and p.rho_kT <= 0
    nJc = 1 if iidm else nJ
    PJc = PJ[:1] if iidm else PJ
    stN = np.zeros(kT.size, int) if iidm else st
    jjc = np.zeros((nJ, 1, 1, 1, 1), int) if iidm else np.arange(nJ)[:, None, None, None, None]
    jjn = np.arange(nJc)[:, None, None, None, None]
    jj = np.arange(nJ)[:, None, None, None, None]
    # joint transition of (y, z): T[z, y, z', y'] = piz[z, z'] lamH[z, y, y']
    T = p.piz[:, None, :, None] * lamH[:, :, None, :]                               # (nZ, nY, nZ', nY')

    def expect(X):
        """E over (y', z') of X(j, e, a, y', z') given today's (y, z): returns (j, e, a, y, z)."""
        return np.einsum("zyxw,jeaxw->jeayz", T, X.transpose(0, 1, 2, 4, 3))       # X indexed (j,e,a,y',z') -> (j,e,a,z',y')

    # experience transitions
    eE = np.minimum(p.e_max, (1 - p.delta_e) * eg[:, None] + p.theta_e * eg[:, None] * hg[None, :] ** p.psi_e)  # (nE, nH)
    kE, wE = _bracket(eg, eE)
    eN = (1 - p.delta_e) * eg
    kN, wN = _bracket(eg, eN)

    VE = np.zeros((nT, nJ, nE, nA, nY, nZ)); VN = np.zeros((nT, nJ, nE, nA, nY, nZ))
    gH = np.zeros((nT, nJ, nE, nA, nY, nZ), np.int16); gAE = np.zeros_like(gH)
    gS = np.zeros_like(gH); gAN = np.zeros_like(gH)
    iters = np.zeros(nT, int)

    for tau in range(nT - 1, -1, -1):
        kap = kbar * (km if tau == 0 else 1.0)
        f_tau = f * (p.home_young_mult if tau == 0 else 1.0)
        fh = f_tau * (1 - hg) ** p.nu_h                                             # (nH,)
        fs = f_tau * (1 - sg) ** p.nu_h                                             # (nS,)
        # flow utilities (independent of V): employed (e, a, y, z, h, a'), non-employed (e, a, y, z, s, a')
        Ra = 1.0 + p.r_a                                                            # gross return on assets
        cE = (w[:, None, None, :, None, None] * hg[None, None, None, None, :, None]
              + fh[None, None, None, None, :, None] + yH[tau][None, None, :, :, None, None]
              + Ra * ag[None, :, None, None, None, None] - ac[None, None, None, None, None, :])
        kap_h = kap * (hg / 0.4) ** p.kappa_h_power if p.kappa_h_power > 0 else np.full(nH, kap)
        vh = p.mu * hg[None, None, None, None, :, None] ** (1 + p.eta) / (1 + p.eta) + kap_h[None, None, None, None, :, None]
        cN = (fs[None, None, None, None, :, None] + yH[tau][None, None, :, :, None, None]
              + Ra * ag[None, :, None, None, None, None] - ac[None, None, None, None, None, :])
        if p.kpr and p.kT_mult:     # transitory shock inside the aggregator: flow by today's node j
            lcv = np.log(np.maximum(cE, 1e-8)) - vh
            flowE = np.stack([np.where(cE > 1e-8, U_kpr(lcv - kT[st_of_j[j]], p.gamma), -1e10) for j in range(nJ)])
            flowN = np.where(cN > 1e-8, U_kpr(np.log(np.maximum(cN, 1e-8)), p.gamma), -1e10)
        elif p.kpr:    # non-separable: the composite log c - v(h) - kappa raised to the CRRA curvature
            flowE = np.where(cE > 1e-8, U_kpr(np.log(np.maximum(cE, 1e-8)) - vh, p.gamma), -1e10)
            flowN = np.where(cN > 1e-8, U_kpr(np.log(np.maximum(cN, 1e-8)), p.gamma), -1e10)
        else:
            flowE = np.where(cE > 1e-8, u(np.maximum(cE, 1e-8), p.gamma), -1e10) - vh
            flowN = np.where(cN > 1e-8, u(np.maximum(cN, 1e-8), p.gamma), -1e10)
        if p.liq < 1.0:      # partial liquidity: a' below the illiquid remainder (1 - liq) a is infeasible
            illiq = np.where(ac[None, :] < (1.0 - p.liq) * ag[:, None] - 1e-12, -1e10, 0.0)      # (a, a')
            flowE = flowE + (illiq[None, None, :, None, None, None, :] if flowE.ndim == 7 else illiq[:, None, None, None, :])
            flowN = flowN + illiq[:, None, None, None, :]
        flowN = np.ascontiguousarray(np.broadcast_to(flowN, (nE, nA, nY, nZ, nS, nAc)))
        # next-age continuation (fixed during the iteration): start-of-period values at age tau+1 by
        # shock node j' (quit decision at node j'), then by today's shock state through PJ
        if tau == nT - 1:
            VRb = np.broadcast_to(VR[None, :, None, None], (nE, nA, nY, nZ))
            Qnext = np.broadcast_to(VRb, (kT.size, nE, nA, nY, nZ)); VNnext = np.broadcast_to(VRb, (nJc, nE, nA, nY, nZ))
        else:
            Qnext = np.maximum(VE[tau + 1][st] - kTdec[:, None, None, None, None], VN[tau + 1][stN]); VNnext = VN[tau + 1][:nJc]
        th = (p.ez_rra - 1.0) * (1.0 - p.beta) if p.ez_rra > 1.0 else 0.0     # risk sensitivity in V units
        vref = float(max(np.max(Qnext), np.max(VNnext))) if th > 0 else 0.0
        Tr = (lambda X: np.exp(-th * (X - vref))) if th > 0 else (lambda X: X)    # to the risk-sensitive domain
        Trinv = (lambda M: vref - np.log(np.maximum(M, 1e-300)) / th) if th > 0 else (lambda M: M)
        EVmax_next = expect(np.tensordot(PJc, Tr(Qnext), axes=(1, 0)))               # (j, e, a, y, z)
        EVN_next = expect(np.tensordot(PJc if nJc > 1 else np.ones((1, 1)), Tr(VNnext), axes=(1, 0)))
        pa = p_age[tau]
        # initial guess
        VEc = (np.maximum(VE[tau + 1], VN[tau + 1]) if tau < nT - 1
               else np.broadcast_to(VR[None, None, :, None, None], (nJ, nE, nA, nY, nZ))).copy()
        VNc = VEc[:nJc].copy()
        hE = np.zeros((nJ, nE, nA, nY, nZ), int); aE = np.zeros_like(hE); sN = np.zeros_like(hE); aN = np.zeros_like(hE)

        def continuation(VEc, VNc):
            # value at the start of a period by shock node j' (before the quit decision), then the
            # expectation over j' given today's state j; all arrays (j, e', a', y, z)
            Q = np.maximum(VEc[st] - kTdec[:, None, None, None, None], VNc[stN])          # (j', e, a, y, z)
            Vmax = np.tensordot(PJc, Tr(Q), axes=(1, 0))          # with ez_rra > 1 everything below is in
            VNj = np.tensordot(PJc if nJc > 1 else np.ones((1, 1)), Tr(VNc), axes=(1, 0))  # the exp(-theta V) domain
            Wmax = (1 - pa) * expect(Vmax) + pa * EVmax_next
            WN = (1 - pa) * expect(VNj) + pa * EVN_next
            WE = (1 - lam_u)[None, None, None, None, :] * Wmax + lam_u[None, None, None, None, :] * WN
            # employed: interpolate WE at e'(e, h): (j, e, h, a', y, z)
            contE = (1 - wE)[None, :, :, None, None, None] * WE[:, kE] + wE[None, :, :, None, None, None] * WE[:, kE + 1]
            # non-employed: interpolate at e'_N(e): (j, e, a', y, z)
            WNs = (1 - wN)[None, :, None, None, None] * WN[:, kN] + wN[None, :, None, None, None] * WN[:, kN + 1]
            WNf = (1 - wN)[None, :, None, None, None] * Wmax[:, kN] + wN[None, :, None, None, None] * Wmax[:, kN + 1]
            # (j, e, s, a', y, z)
            contN = ((1 - pi_f)[None, None, :, None, None, :] * WNs[:, :, None]
                     + pi_f[None, None, :, None, None, :] * WNf[:, :, None])
            return on_choice(Trinv(contE)), on_choice(Trinv(contN))   # certainty equivalents (identity if theta = 0)

        ie, ia, iy, iz = np.ogrid[:nE, :nA, :nY, :nZ]
        for it in range(p.max_iter):
            contE, contN = continuation(VEc, VNc)
            # maximisation
            totE = (flowE if flowE.ndim == 7 else flowE[None]) + beta * contE.transpose(0, 1, 4, 5, 2, 3)[:, :, None]   # (j, e, a, y, z, h, a')
            flatE = totE.reshape(nJ, nE, nA, nY, nZ, -1)
            idxE = flatE.argmax(axis=-1)
            hE, aE = np.divmod(idxE, nAc)
            VEn = np.take_along_axis(flatE, idxE[..., None], axis=-1)[..., 0]
            totN = flowN[None] + beta * contN.transpose(0, 1, 4, 5, 2, 3)[:, :, None]   # (j, e, a, y, z, s, a')
            flatN = totN.reshape(nJc, nE, nA, nY, nZ, -1)
            idxN = flatN.argmax(axis=-1)
            sN, aN = np.divmod(idxN, nAc)
            VNn = np.take_along_axis(flatN, idxN[..., None], axis=-1)[..., 0]
            err = max(np.max(np.abs(VEn - VEc)), np.max(np.abs(VNn - VNc)))
            VEc, VNc = VEn, VNn
            if err < p.vf_tol:
                break
            # Howard evaluation with fixed policies
            fE = (flowE[jj, ie, ia, iy, iz, hE, aE] if flowE.ndim == 7 else flowE[ie, ia, iy, iz, hE, aE]); fN = flowN[ie, ia, iy, iz, sN, aN]
            for _ in range(p.howard_steps):
                contE, contN = continuation(VEc, VNc)
                VEc = fE + beta * contE[jjc, ie, hE, aE, iy, iz]
                VNc = fN + beta * contN[jjn, ie, sN, aN, iy, iz]
        iters[tau] = it + 1
        VE[tau], VN[tau] = VEc, VNc                         # V^N broadcast over the nodes when nJc = 1
        gH[tau], gAE[tau], gS[tau], gAN[tau] = hE, aE, sN, aN
        if verbose:
            print(f"    age {tau}: {it + 1} iterations, err {err:.2e}")
    return VE, VN, gH, gAE, gS, gAN, iters


def _solve_one(args):
    p, om, kb, km_, VR, zh = args
    return solve_type(p, om, kb, km_, VR, zh=zh)


def solve_all(p: FinalParams, verbose=False, n_jobs: int | None = None) -> FinalSolution:
    import os
    if n_jobs is None:
        n_jobs = int(os.environ.get("KEAM_NJOBS", os.cpu_count() or 1))
    omega, kbar, km, zh = make_types4(p)
    nK = omega.size
    VR = solve_retirement(p)
    nJ = p.kT_chain()[1]
    shape = (nK, 3, p.nE, p.nA, 3, 2, nJ)
    VE = np.zeros(shape); VN = np.zeros(shape)
    gH = np.zeros(shape, np.float32); gAE = np.zeros(shape, np.float32)
    gS = np.zeros(shape, np.float32); gAN = np.zeros(shape, np.float32)
    t0 = time.time(); iters = np.zeros((nK, 3), int)
    jobs = [(p, omega[k], kbar[k], km[k], VR, zh[k]) for k in range(nK)]
    if n_jobs > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(n_jobs) as pool:
            results = pool.map(_solve_one, jobs, chunksize=max(1, nK // (4 * n_jobs)))
    else:
        results = [_solve_one(j) for j in jobs]
    for k, (ve, vn, h, aE, s, aN, its) in enumerate(results):
        # per-type arrays are (nT, nJ, nE, nA, nY, nZ); the shock state is stored last
        VE[k], VN[k], iters[k] = np.moveaxis(ve, 1, -1), np.moveaxis(vn, 1, -1), its
        gH[k] = np.moveaxis(p.hgrid[h], 1, -1); gAE[k] = np.moveaxis(p.agrid_c[aE], 1, -1)
        gS[k] = np.moveaxis(p.sgrid[s], 1, -1); gAN[k] = np.moveaxis(p.agrid_c[aN], 1, -1)
    if verbose:
        print(f"solved {nK} types in {time.time() - t0:.0f}s; max iterations {iters.max()}")
    return FinalSolution(params=p, omega=omega, kbar=kbar, km=km, zh=zh, VE=VE, VN=VN, gH=gH, gAE=gAE,
                         gS=gS, gAN=gAN, VR=VR, info=dict(seconds=time.time() - t0, iters=iters))
