"""Value function iteration for the final model (monthly, assets, per-type solve).

Per type the state is (age tau, experience e, assets a, husband y, aggregate z, cost-shock state j)
and there are two value functions: employed V^E (chooses hours h and savings a', may quit) and
non-employed V^N (chooses search s and savings a').  Choices are on grids; continuation values are
linearly interpolated in e'.  Modified policy iteration (maximisation step followed by
`howard_steps` evaluation steps) is used for speed.

The cost-of-work shock kappa_T is realised at the start of the period, before the quit decision:
an employed woman (or a job finder) with shock node j' gets max{V^E_j' - kappa_j', V^N_j'}.  With an
iid shock (rho_kT = 0) the values do not depend on the shock and a single shock state is carried
(last axis of length 1); with persistence the state carries the current node (last axis n_kT).
"""
from __future__ import annotations
from dataclasses import dataclass
import time
import numpy as np
from .params import FinalParams, make_types


@dataclass
class FinalSolution:
    params: FinalParams
    omega: np.ndarray; kbar: np.ndarray; km: np.ndarray
    VE: np.ndarray; VN: np.ndarray          # (nK, nT, nE, nA, nY, nZ, nJ); nJ = 1 for an iid shock
    gH: np.ndarray; gAE: np.ndarray         # hours, savings when employed (values, not indices)
    gS: np.ndarray; gAN: np.ndarray         # search, savings when non-employed
    VR: np.ndarray                          # retirement value on the asset grid
    info: dict


def u(c, gamma):
    return np.power(c, 1.0 - gamma) / (1.0 - gamma)


def _bracket(grid, x):
    k = np.clip(np.searchsorted(grid, x, side="right") - 1, 0, grid.size - 2)
    w = np.clip((x - grid[k]) / (grid[k + 1] - grid[k]), 0.0, 1.0)
    return k, w


def solve_retirement(p: FinalParams):
    a = p.agrid; nA = a.size
    yR = p.pension * p.yH_age[1]
    c = yR + a[:, None] - a[None, :]                        # (a, a')
    flow = np.where(c > 1e-8, u(np.maximum(c, 1e-8), p.gamma), -1e10)
    V = flow.max(axis=1) / (1 - p.beta * (1 - p.death))
    for _ in range(5000):
        Vn = (flow + p.beta * (1 - p.death) * V[None, :]).max(axis=1)
        if np.max(np.abs(Vn - V)) < 1e-9:
            V = Vn; break
        V = Vn
    return V


def solve_type(p: FinalParams, omega: float, kbar: float, km: float, VR: np.ndarray, verbose=False):
    nE, nA, nH, nS = p.nE, p.nA, p.nH, p.nS
    nY, nZ, nT = 3, 2, 3
    eg, ag, hg, sg = p.egrid, p.agrid, p.hgrid, p.sgrid
    beta = p.beta
    phi = np.array([1.0, p.phi_rec])
    w = phi[None, :] * p.tau_w * omega * (1 + p.gam_e * eg[:, None] ** p.xi)      # (nE, nZ)
    f = p.ybar_h + p.z_h * omega ** p.alpha_h
    yH = p.y_husband()                                                              # (nT, nY, nZ)
    lamH = p.lamH()                                                                 # (nZ, nY, nY)
    lam_u = np.asarray(p.lam_u); lam_f = np.asarray(p.lam_f)
    pi_f = lam_f[None, :] * sg[:, None] ** p.nu                                     # (nS, nZ)
    p_age = p.p_age
    kT, nJ, PJ = p.kT_chain()                                                       # nodes, states, (nJ, n_kT)
    st = np.arange(kT.size) if nJ > 1 else np.zeros(kT.size, int)                   # state reached at node j'
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
        cE = (w[:, None, None, :, None, None] * hg[None, None, None, None, :, None]
              + fh[None, None, None, None, :, None] + yH[tau][None, None, :, :, None, None]
              + ag[None, :, None, None, None, None] - ag[None, None, None, None, None, :])
        kap_h = kap * (hg / 0.4) ** p.kappa_h_power if p.kappa_h_power > 0 else np.full(nH, kap)
        flowE = np.where(cE > 1e-8, u(np.maximum(cE, 1e-8), p.gamma), -1e10) \
            - p.mu * hg[None, None, None, None, :, None] ** (1 + p.eta) / (1 + p.eta) - kap_h[None, None, None, None, :, None]
        cN = (fs[None, None, None, None, :, None] + yH[tau][None, None, :, :, None, None]
              + ag[None, :, None, None, None, None] - ag[None, None, None, None, None, :])
        flowN = np.where(cN > 1e-8, u(np.maximum(cN, 1e-8), p.gamma), -1e10)
        flowN = np.ascontiguousarray(np.broadcast_to(flowN, (nE, nA, nY, nZ, nS, nA)))
        # next-age continuation (fixed during the iteration): start-of-period values at age tau+1 by
        # shock node j' (quit decision at node j'), then by today's shock state through PJ
        if tau == nT - 1:
            VRb = np.broadcast_to(VR[None, :, None, None], (nE, nA, nY, nZ))
            Qnext = np.broadcast_to(VRb, (kT.size, nE, nA, nY, nZ)); VNnext = np.broadcast_to(VRb, (nJ, nE, nA, nY, nZ))
        else:
            Qnext = np.maximum(VE[tau + 1][st] - kT[:, None, None, None, None], VN[tau + 1][st]); VNnext = VN[tau + 1]
        EVmax_next = expect(np.tensordot(PJ, Qnext, axes=(1, 0)))                   # (j, e, a, y, z)
        EVN_next = expect(np.tensordot(PJ if nJ > 1 else np.ones((1, 1)), VNnext, axes=(1, 0)))
        pa = p_age[tau]
        # initial guess
        VEc = (np.maximum(VE[tau + 1], VN[tau + 1]) if tau < nT - 1
               else np.broadcast_to(VR[None, None, :, None, None], (nJ, nE, nA, nY, nZ))).copy()
        VNc = VEc.copy()
        hE = np.zeros((nJ, nE, nA, nY, nZ), int); aE = np.zeros_like(hE); sN = np.zeros_like(hE); aN = np.zeros_like(hE)

        def continuation(VEc, VNc):
            # value at the start of a period by shock node j' (before the quit decision), then the
            # expectation over j' given today's state j; all arrays (j, e', a', y, z)
            Q = np.maximum(VEc[st] - kT[:, None, None, None, None], VNc[st])           # (j', e, a, y, z)
            Vmax = np.tensordot(PJ, Q, axes=(1, 0))
            VNj = np.tensordot(PJ if nJ > 1 else np.ones((1, 1)), VNc, axes=(1, 0))
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
            return contE, contN

        ie, ia, iy, iz = np.ogrid[:nE, :nA, :nY, :nZ]
        for it in range(p.max_iter):
            contE, contN = continuation(VEc, VNc)
            # maximisation
            totE = flowE[None] + beta * contE.transpose(0, 1, 4, 5, 2, 3)[:, :, None]   # (j, e, a, y, z, h, a')
            flatE = totE.reshape(nJ, nE, nA, nY, nZ, -1)
            idxE = flatE.argmax(axis=-1)
            hE, aE = np.divmod(idxE, nA)
            VEn = np.take_along_axis(flatE, idxE[..., None], axis=-1)[..., 0]
            totN = flowN[None] + beta * contN.transpose(0, 1, 4, 5, 2, 3)[:, :, None]   # (j, e, a, y, z, s, a')
            flatN = totN.reshape(nJ, nE, nA, nY, nZ, -1)
            idxN = flatN.argmax(axis=-1)
            sN, aN = np.divmod(idxN, nA)
            VNn = np.take_along_axis(flatN, idxN[..., None], axis=-1)[..., 0]
            err = max(np.max(np.abs(VEn - VEc)), np.max(np.abs(VNn - VNc)))
            VEc, VNc = VEn, VNn
            if err < p.vf_tol:
                break
            # Howard evaluation with fixed policies
            fE = flowE[ie, ia, iy, iz, hE, aE]; fN = flowN[ie, ia, iy, iz, sN, aN]
            for _ in range(p.howard_steps):
                contE, contN = continuation(VEc, VNc)
                VEc = fE + beta * contE[jj, ie, hE, aE, iy, iz]
                VNc = fN + beta * contN[jj, ie, sN, aN, iy, iz]
        iters[tau] = it + 1
        VE[tau], VN[tau] = VEc, VNc
        gH[tau], gAE[tau], gS[tau], gAN[tau] = hE, aE, sN, aN
        if verbose:
            print(f"    age {tau}: {it + 1} iterations, err {err:.2e}")
    return VE, VN, gH, gAE, gS, gAN, iters


def _solve_one(args):
    p, om, kb, km_, VR = args
    return solve_type(p, om, kb, km_, VR)


def solve_all(p: FinalParams, verbose=False, n_jobs: int | None = None) -> FinalSolution:
    import os
    if n_jobs is None:
        n_jobs = int(os.environ.get("KEAM_NJOBS", os.cpu_count() or 1))
    omega, kbar, km = make_types(p)
    nK = omega.size
    VR = solve_retirement(p)
    nJ = p.kT_chain()[1]
    shape = (nK, 3, p.nE, p.nA, 3, 2, nJ)
    VE = np.zeros(shape); VN = np.zeros(shape)
    gH = np.zeros(shape, np.float32); gAE = np.zeros(shape, np.float32)
    gS = np.zeros(shape, np.float32); gAN = np.zeros(shape, np.float32)
    t0 = time.time(); iters = np.zeros((nK, 3), int)
    jobs = [(p, omega[k], kbar[k], km[k], VR) for k in range(nK)]
    if n_jobs > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(n_jobs) as pool:
            results = pool.map(_solve_one, jobs, chunksize=max(1, nK // (4 * n_jobs)))
    else:
        results = [_solve_one(j) for j in jobs]
    for k, (ve, vn, h, aE, s, aN, its) in enumerate(results):
        # per-type arrays are (nT, nJ, nE, nA, nY, nZ); the shock state is stored last
        VE[k], VN[k], iters[k] = np.moveaxis(ve, 1, -1), np.moveaxis(vn, 1, -1), its
        gH[k] = np.moveaxis(p.hgrid[h], 1, -1); gAE[k] = np.moveaxis(p.agrid[aE], 1, -1)
        gS[k] = np.moveaxis(p.sgrid[s], 1, -1); gAN[k] = np.moveaxis(p.agrid[aN], 1, -1)
    if verbose:
        print(f"solved {nK} types in {time.time() - t0:.0f}s; max iterations {iters.max()}")
    return FinalSolution(params=p, omega=omega, kbar=kbar, km=km, VE=VE, VN=VN, gH=gH, gAE=gAE,
                         gS=gS, gAN=gAN, VR=VR, info=dict(seconds=time.time() - t0, iters=iters))
