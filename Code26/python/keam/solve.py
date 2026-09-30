"""Value function iteration: translation of SimplerMod_May17_splines.m (lines 133-393).

State arrays are shaped (nI, nT, nE, nY, nZ) exactly like the MATLAB arrays so
that the output can be compared element-by-element with the saved
``policies.mat`` / ``Vfuns.mat``.

For each (type i, age it) the MATLAB code loops over (iy, iz, ie) and runs a
bisection on hours (employed) or search (non-employed).  Here the bisection is
vectorised over the nE*nY*nZ states of the block; the per-state stopping rules
and bound-update rules are identical to the MATLAB ones.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import time
import numpy as np

from .params import Params, Options
from . import functions as F


@dataclass
class Solution:
    VE: np.ndarray
    VU: np.ndarray
    V: np.ndarray
    gH: np.ndarray
    gS: np.ndarray
    gQ: np.ndarray
    n_iter: np.ndarray            # VFI iterations per (i, it)
    last_err: np.ndarray          # last sup-norm change per (i, it)
    params: Params = None
    options: Options = None
    info: dict = field(default_factory=dict)


# ----------------------------------------------------------------------------
def _bracket_above(egrid, ee):
    """MATLAB: iee = min(nE-1, find(egrid > ee, 1, 'first')) with the two guards.

    Returns 0-based index k such that the segment [k, k+1] is used.  Note that
    egrid[k] >= ee unless ee is above the grid (then k = nE-2).
    """
    nE = egrid.size
    k = np.searchsorted(egrid, ee, side="right")   # number of grid points <= ee
    k = np.minimum(k, nE - 2)
    k = np.where(ee < egrid[0], 0, k)
    return k


def _bracket_contain(egrid, ee):
    """0-based k with egrid[k] <= ee < egrid[k+1] (clipped to [0, nE-2]) and the
    linear weight on k+1.  Extrapolates linearly outside the grid."""
    nE = egrid.size
    k = np.searchsorted(egrid, ee, side="right") - 1
    k = np.clip(k, 0, nE - 2)
    w = (ee - egrid[k]) / (egrid[k + 1] - egrid[k])
    return k, w


# ----------------------------------------------------------------------------
def solve(p: Params, opts: Options | None = None, verbose: bool = False) -> Solution:
    opts = opts or Options.faithful()
    nI, nT, nE, nY, nZ = p.nI, p.nT, p.nE, p.nY, p.nZ
    egrid = p.egrid
    ftW = p.ftW
    kap = p.kap_matrix(opts)
    piT = p.piT_used(opts)                     # length nT-1
    beta, nu, nu_h = p.beta, p.nu, p.nu_h
    lossW, findW, piz, lamH = p.lossW, p.findW, p.piz, p.lamH
    wageH, BCwage, ym = p.wageH, p.BCwage, p.ym
    vf_tol = p.VFtol if opts.loose_vf_tolerance else opts.vf_tol_corrected
    max_iter = 10 ** 9 if opts.loose_vf_tolerance else opts.max_vf_iter_corrected

    shape = (nI, nT, nE, nY, nZ)
    VE = np.zeros(shape); VU0 = np.zeros(shape); V = np.zeros(shape)
    gH = np.ones(shape); gS = np.ones(shape); gQ = np.zeros(shape, dtype=np.int8)
    n_iter = np.zeros((nI, nT), int); last_err = np.full((nI, nT), np.nan)

    # --- initial guess (lines 157-169) ------------------------------------
    IE, IY, IZ = np.meshgrid(np.arange(nE), np.arange(nY), np.arange(nZ), indexing="ij")
    for i in range(nI):
        hp = F.hprod(p, ftW[i])
        w = F.wage(p, ftW[i], egrid)
        for it in range(nT - 1):
            ve0 = (F.utilC(p, wageH[it] * BCwage[IZ] * ym[IY] + BCwage[IZ] * w[IE] * 0.5
                           + hp * (1 - 0.5) ** nu_h) + F.utilL(p, 0.5) - kap[it, i])
            vu0 = F.utilC(p, hp * 0.2 + ym[IY] * wageH[nT - 2] * BCwage[IZ]
                          * lamH[IZ, np.maximum(IY - 1, 0), 0])
            VE[i, it] = ve0; VU0[i, it] = vu0; V[i, it] = np.maximum(ve0, vu0)

    # expectation weights W[iy, iz, iyy, izz] = beta * piz[iz, izz] * lamH[izz, iy, iyy]
    W = beta * np.einsum("zw,wyx->yzxw", piz, lamH)        # (nY, nZ, nY, nZ)
    Ws = W[IY.ravel(), IZ.ravel()]                           # (S, nY, nZ)
    iy_s, iz_s, ie_s = IY.ravel(), IZ.ravel(), IE.ravel()
    S = iy_s.size
    lossWs = lossW[None, None, :]                            # (1,1,nZ)
    findWs = findW[None, None, :]

    t0 = time.time()
    for it in range(nT - 2, -1, -1):
        pT = piT[it]
        last_age = (it + 1 > nT - 2)                         # it+1 is retirement
        for i in range(nI):
            if i > 0:                                        # lines 179-181
                V[i, it] = V[i - 1, it]; VU0[i, it] = VU0[i - 1, it]; VE[i, it] = VE[i - 1, it]
            elif it < nT - 2 and opts.loose_vf_tolerance:
                # MATLAB lines 286/360/378 assign the FULL arrays V=V0, VU0=VU1 where
                # V0 and VU1 were initialised to zero, so the flow-utility initial
                # guess (lines 157-169) survives only for the very first block.
                V[i, it] = 0.0; VU0[i, it] = 0.0; VE[i, it] = 0.0
            hp = F.hprod(p, ftW[i])
            w_e = F.wage(p, ftW[i], egrid)
            ymH = wageH[it] * BCwage[IZ] * ym[IY]            # (nE,nY,nZ) husband income
            ymH_s = ymH.ravel()
            w_s = w_e[ie_s]; e_s = egrid[ie_s]
            kap_it = kap[it, i]

            for it_v in range(max_iter):
                Vcur, VUcur = V[i, it], VU0[i, it]           # (nE,nY,nZ)
                Vnext, VUnext = V[i, it + 1], VU0[i, it + 1]

                # ============ EMPLOYED: bisection on hours (lines 190-283) ======
                hlow = np.zeros(S); hhigh = np.ones(S); h = np.zeros(S)
                active = np.ones(S, bool)
                VE1 = np.zeros(S); gH_blk = np.zeros(S)
                for ith in range(1, p.maxHSiter):
                    a = active
                    h[a] = 0.5 * (hlow[a] + hhigh[a])
                    ha = h[a]
                    YY = ymH_s[a] + BCwage[iz_s[a]] * w_s[a] * ha + hp * (1 - ha) ** nu_h
                    dUCdh = F.du_dC(p, YY) * (BCwage[iz_s[a]] * w_s[a] - hp * nu_h * (1 - ha) ** (nu_h - 1))
                    dULdh = F.du_dL(p, ha)
                    ee = F.exp2(p, e_s[a], ha)
                    de_dh = F.de2_dh(p, e_s[a], ha)
                    de_dh_typo = ((p.alpha_h if opts.alpha_h_typo else p.alpha_e)
                                  * p.psi * ha ** (p.psi - 1) * e_s[a])
                    if opts.roundup_continuation:
                        k = _bracket_above(egrid, ee); wk = None
                    else:
                        k, wk = _bracket_contain(egrid, ee)
                    dg = (egrid[k + 1] - egrid[k])
                    if opts.derivative_at_current_state:
                        ya, za = iy_s[a], iz_s[a]
                        dV1 = (Vcur[k + 1, ya, za] - Vcur[k, ya, za]) / dg * de_dh
                        dVU1 = (VUcur[k + 1, ya, za] - VUcur[k, ya, za]) / dg * de_dh_typo
                        if last_age:
                            dV2 = np.zeros_like(dV1); dVU2 = np.zeros_like(dV1)
                        else:
                            dV2 = (Vnext[k + 1, ya, za] - Vnext[k, ya, za]) / dg * de_dh
                            dVU2 = (VUnext[k + 1, ya, za] - VUnext[k, ya, za]) / dg * de_dh
                        dV1 = dV1[:, None, None]; dVU1 = dVU1[:, None, None]
                        dV2 = dV2[:, None, None]; dVU2 = dVU2[:, None, None]
                    else:
                        dg3 = dg[:, None, None]
                        dV1 = (Vcur[k + 1] - Vcur[k]) / dg3 * de_dh[:, None, None]
                        dVU1 = (VUcur[k + 1] - VUcur[k]) / dg3 * de_dh_typo[:, None, None]
                        if last_age:
                            dV2 = np.zeros_like(dV1); dVU2 = np.zeros_like(dV1)
                        else:
                            dV2 = (Vnext[k + 1] - Vnext[k]) / dg3 * de_dh[:, None, None]
                            dVU2 = (VUnext[k + 1] - VUnext[k]) / dg3 * de_dh[:, None, None]
                    X = ((1 - pT) * ((1 - lossWs) * dV1 + lossWs * dVU1)
                         + pT * ((1 - lossWs) * dV2 + lossWs * dVU2))
                    TD = (Ws[a] * X).sum(axis=(1, 2)) + dUCdh + dULdh
                    done = (np.abs(TD) < p.bisectTol) | (ith > p.maxHSiter - 2) | (ha < 0.001)
                    idx = np.flatnonzero(a)[done]
                    if done.any():
                        gH_blk[idx] = ha[done]
                        kd = k[done]
                        if opts.roundup_continuation:
                            EVe, EVu = Vcur[kd], VUcur[kd]
                            EVe2, EVu2 = Vnext[kd], VUnext[kd]
                        else:
                            wd = wk[done][:, None, None]
                            EVe = (1 - wd) * Vcur[kd] + wd * Vcur[kd + 1]
                            EVu = (1 - wd) * VUcur[kd] + wd * VUcur[kd + 1]
                            EVe2 = (1 - wd) * Vnext[kd] + wd * Vnext[kd + 1]
                            EVu2 = (1 - wd) * VUnext[kd] + wd * VUnext[kd + 1]
                        Xv = ((1 - pT) * ((1 - lossWs) * EVe + lossWs * EVu)
                              + pT * ((1 - lossWs) * EVe2 + lossWs * EVu2))
                        V1 = (Ws[idx] * Xv).sum(axis=(1, 2))
                        VE1[idx] = V1 + F.utilC(p, YY[done]) + F.utilL(p, ha[done]) - kap_it
                    nd = ~done
                    idx_nd = np.flatnonzero(a)[nd]
                    up = TD[nd] > 0
                    hlow[idx_nd[up]] = 0.8 * h[idx_nd[up]] + 0.2 * hlow[idx_nd[up]]
                    hhigh[idx_nd[~up]] = 0.8 * h[idx_nd[~up]] + 0.2 * hhigh[idx_nd[~up]]
                    active[idx] = False
                    if not active.any():
                        break
                # MATLAB loop order (iy, iz, ie): the last state visited is (ie=nE, iy=nY, iz=nZ)
                h_stale = gH_blk[np.ravel_multi_index((nE - 1, nY - 1, nZ - 1), (nE, nY, nZ))]

                # ============ NON-EMPLOYED: bisection on search (lines 293-357) ====
                def solve_search(ku, wu=None):
                    """Vectorised search bisection given the e' bracket index ku (and
                    interpolation weight wu when not rounding up).  Returns (s, VU)."""
                    if wu is None:
                        EVe_u, EVu_u = Vcur[ku], VUcur[ku]                # (S,nY,nZ)
                        EVe2_u, EVu2_u = Vnext[ku], VUnext[ku]
                    else:
                        wu3 = wu[:, None, None]
                        EVe_u = (1 - wu3) * Vcur[ku] + wu3 * Vcur[ku + 1]
                        EVu_u = (1 - wu3) * VUcur[ku] + wu3 * VUcur[ku + 1]
                        EVe2_u = (1 - wu3) * Vnext[ku] + wu3 * Vnext[ku + 1]
                        EVu2_u = (1 - wu3) * VUnext[ku] + wu3 * VUnext[ku + 1]
                    gain1 = np.maximum(EVe_u - EVu_u, 0.0)
                    gain2 = np.maximum(EVe2_u - EVu2_u, 0.0)
                    Wgain = (Ws * ((1 - pT) * findWs * gain1 + pT * findWs * gain2)).sum(axis=(1, 2))
                    hlow = np.zeros(S); hhigh = np.ones(S); s = np.zeros(S)
                    active = np.ones(S, bool)
                    VU1 = np.zeros(S); gS_blk = np.zeros(S)
                    for ith in range(1, p.maxHSiter):
                        a = active
                        s[a] = 0.5 * (hlow[a] + hhigh[a])
                        sa = s[a]
                        YY = ymH_s[a] + hp * (1 - sa) ** nu_h
                        dUCdh = -F.du_dC(p, YY) * (hp * nu_h * (1 - sa) ** (nu_h - 1))
                        TD = Wgain[a] * nu * sa ** (nu - 1) + dUCdh
                        done = (np.abs(TD) < p.bisectTol) | (sa < 0.0001) | (ith > p.maxHSiter - 2)
                        idx = np.flatnonzero(a)[done]
                        if done.any():
                            sd = sa[done]
                            gS_blk[idx] = sd
                            pf = (findWs * sd[:, None, None] ** nu)            # (nd,1,nZ)
                            Xv = ((1 - pT) * ((1 - pf) * EVu_u[idx] + pf * EVe_u[idx])
                                  + pT * ((1 - pf) * EVu2_u[idx] + pf * EVe2_u[idx]))
                            V1 = (Ws[idx] * Xv).sum(axis=(1, 2))
                            VU1[idx] = V1 + F.utilC(p, YY[done]) + F.utilL(p, p.lhome)
                        nd = ~done
                        idx_nd = np.flatnonzero(a)[nd]
                        up = TD[nd] > 0
                        hlow[idx_nd[up]] = 0.8 * s[idx_nd[up]] + 0.2 * hlow[idx_nd[up]]
                        hhigh[idx_nd[~up]] = 0.8 * s[idx_nd[~up]] + 0.2 * hhigh[idx_nd[~up]]
                        active[idx] = False
                        if not active.any():
                            break
                    return gS_blk, VU1

                if not opts.stale_h_nonemployed:
                    ee_u = F.exp2(p, e_s, 0.0)
                    if opts.roundup_continuation:
                        gS_blk, VU1 = solve_search(_bracket_above(egrid, ee_u))
                    else:
                        ku, wu = _bracket_contain(egrid, ee_u)
                        gS_blk, VU1 = solve_search(ku, wu)
                else:
                    # MATLAB re-uses the scalar `h` as the search variable, so the
                    # experience update of each state uses the converged search of the
                    # PREVIOUS state in loop order (iy, iz, ie); the first state uses
                    # the last employed hours.  Solve for the candidate brackets
                    # {ie-1, ie, ie+1}, then walk the states in MATLAB order.
                    cands = {}
                    for off in (-1, 0, 1):
                        kc = np.clip(ie_s + off, 0, nE - 2)
                        cands[off] = solve_search(kc)
                    gS_blk = np.zeros(S); VU1 = np.zeros(S)
                    h_prev = h_stale
                    for iy_m in range(nY):
                        for iz_m in range(nZ):
                            for ie_m in range(nE):
                                sidx = np.ravel_multi_index((ie_m, iy_m, iz_m), (nE, nY, nZ))
                                k_act = int(_bracket_above(egrid, F.exp2(p, egrid[ie_m], h_prev)))
                                off = k_act - int(np.clip(ie_m, 0, nE - 2))
                                if off in cands and int(np.clip(ie_m + off, 0, nE - 2)) == k_act:
                                    gS_blk[sidx], VU1[sidx] = cands[off][0][sidx], cands[off][1][sidx]
                                else:  # rare: solve this single state on the fly
                                    kk = np.full(S, k_act)
                                    s1, v1 = solve_search(kk)
                                    gS_blk[sidx], VU1[sidx] = s1[sidx], v1[sidx]
                                h_prev = gS_blk[sidx]

                # ============ update (lines 360-381) ==============================
                VE_new = VE1.reshape(nE, nY, nZ)
                VU_new = VU1.reshape(nE, nY, nZ)
                V0 = np.maximum(VU_new, VE_new)
                if opts.loose_vf_tolerance:
                    err = np.abs(np.max(V0 - Vcur))      # MATLAB line 377: abs(max(...)), not a sup norm
                else:
                    err = np.max(np.abs(V0 - Vcur))
                VE[i, it] = VE_new; VU0[i, it] = VU_new; V[i, it] = V0
                gH[i, it] = gH_blk.reshape(nE, nY, nZ)
                gS[i, it] = gS_blk.reshape(nE, nY, nZ)
                gQ[i, it] = (VU_new > VE_new)
                if err < vf_tol:
                    break
            n_iter[i, it] = it_v + 1; last_err[i, it] = err
            if verbose:
                print(f"  age {it} type {i:2d}: {it_v + 1:4d} iterations, last change {err:.3e}, "
                      f"elapsed {time.time() - t0:6.1f}s")

    return Solution(VE=VE, VU=VU0, V=V, gH=gH, gS=gS, gQ=gQ, n_iter=n_iter,
                    last_err=last_err, params=p, options=opts,
                    info={"seconds": time.time() - t0})
