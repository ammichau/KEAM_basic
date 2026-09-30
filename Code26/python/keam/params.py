"""Model parameters and behavioural options.

The default values reproduce ``Code26/SimplerMod_May17_splines.m`` (and the
saved ``Code26/Solution/Baseline/paras.mat``).  Comments give the MATLAB line
number of each assignment so that departures can be traced.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields, replace
import numpy as np


@dataclass
class Options:
    """Switches controlling departures between the MATLAB code and the model text.

    Each flag defaults to the MATLAB behaviour (``faithful``).  Setting a flag to
    ``False`` replaces the MATLAB behaviour with the behaviour described in the
    paper/slides.  The DEPARTURES.md entry number is given for each flag.
    """

    # ---- solver (SimplerMod_May17_splines.m) ------------------------------
    #: D-S1  piT(it) linear index into the 3x4 ageing matrix (young age w.p.
    #: 0.982 per quarter, middle/old never age).  False: use piT(it, it+1).
    aging_linear_index: bool = True
    #: D-S2  e' for the non-employed uses the stale hours value left over from
    #: the employed loop instead of h=0 (experience grows while non-employed).
    stale_h_nonemployed: bool = True
    #: D-S3  continuation values evaluated at the first grid point ABOVE e'
    #: (spline at a knot) instead of interpolated at e'.
    roundup_continuation: bool = True
    #: D-S4  marginal value of experience evaluated at today's (y, z) instead
    #: of tomorrow's (y', z') inside the expectation.
    derivative_at_current_state: bool = True
    #: D-S5  job-loss branch of the experience derivative scaled by alpha_h
    #: (home production curvature) instead of alpha_e.
    alpha_h_typo: bool = True
    #: D-S6  fixed-cost offset kapbar applied by WAGE type over columns 1:5
    #: instead of by fixed-kappa type over columns 1:4.
    kapbar_indexing_bug: bool = True
    #: D-S7  young-only cost multiplier applied at every age.
    young_cost_all_ages: bool = True
    #: D-S8  loose VFI tolerance (0.1 in value units) and previous type's
    #: value as the initial guess.  False: tolerance `vf_tol_corrected`.
    loose_vf_tolerance: bool = True
    vf_tol_corrected: float = 1e-6
    max_vf_iter_corrected: int = 3000

    # ---- simulator (SimplerMod_May17_sim.m) --------------------------------
    #: D-M1 interpolation weights in experience are reversed (weight on the
    #: lower grid point is (e - e_lo)/de).
    reversed_e_weights: bool = True
    #: D-M2 quit policy interpolation across types uses reversed weights and
    #: the interpolated 0/1 rule is applied as a quit probability.
    reversed_type_weights_gQ: bool = True
    #: D-M3 same seed (222) used for fixed and young kappa draws.
    same_seed_kappa_draws: bool = True
    #: D-M4 wife's simulated wage omits the recession penalty BCwage.
    no_bcwage_in_sim: bool = True
    #: D-M5 everyone starts employed at birth (Nw_init draw overridden).
    all_employed_at_birth: bool = True
    #: D-M6 husband transitions and wife job finding use last period's z.
    lagged_z_in_sim: bool = True
    #: D-M7 quit checked before job loss (P(loss) = (1-Q) lambda_u).
    quit_before_loss: bool = True

    # ---- statistics ---------------------------------------------------------
    #: D-R1 "expansion" indicator equals one in every period.
    expansion_includes_recessions: bool = True
    #: D-R2 women outside the four career definitions dropped from the shares.
    drop_undefined_careers: bool = True
    #: D-R3 the last (never simulated) birth cohort keeps its initial employment
    #: draws in the arrays and is counted in the cycle statistics' numerators.
    phantom_last_cohort: bool = True

    @classmethod
    def faithful(cls) -> "Options":
        return cls()

    @classmethod
    def corrected(cls, **overrides) -> "Options":
        kw = {f.name: False for f in fields(cls) if f.type == "bool"}
        kw.update(overrides)
        return cls(**kw)

    def with_(self, **overrides) -> "Options":
        return replace(self, **overrides)


@dataclass
class Params:
    # --- preferences (lines 30-34, 50) -------------------------------------
    beta: float = 0.97      # discount factor, quarterly (line 30)
    crra: float = 2.0       # curvature of consumption utility (31)
    eta: float = 1.4        # curvature of hours disutility (32)
    mu: float = 0.5         # weight on hours disutility (33)
    phi_c: float = 0.5      # scale on consumption utility (50)
    r: float = 0.0          # interest rate, unused: no savings (34)
    # --- experience technology (35-38) --------------------------------------
    alpha_e: float = 0.026  # theta_e in the slides (35)
    delta_e: float = 0.005  # depreciation (36)
    psi: float = 0.8        # psi_e curvature on hours (37)
    xi: float = 0.85        # xi_e curvature of experience in wage (38)
    # --- home production and search (40-45) ----------------------------------
    lhome: float = 0.0      # hours disutility argument when non-employed (40)
    z_h: float = 0.45       # (41)
    ybar_h: float = 0.1     # (42)
    alpha_h: float = 0.2    # (43)
    nu: float = 0.4         # search efficiency exponent (44)
    nu_h: float = 0.5       # home production curvature in hours (45)
    # --- wages (47-49) ---------------------------------------------------------
    ftWbase: float = 1.0    # high fixed wage type (47)
    ftWlow_scale: float = 0.5  # low type = ftWbase * 0.5 (line 90)
    gam_e: float = 0.5      # weight on experience in wage (48)
    tau_wf: float = 0.8     # gender wage gap (49)
    # --- fixed cost of work (51-52, 92-107) -----------------------------------
    kapbase: float = 0.075  # (51)
    kaplow: float = -0.04   # (52)
    kap_mid_scale: float = 1.0   # kap(2,:) multiplier (95)
    kap_old_scale: float = 1.0   # kap(3,:) multiplier (96)
    kayscale: tuple = (1.0, 0.9, 0.85, 0.5)  # (97)
    # --- experiment scalings (shell) -------------------------------------------
    kapscale: float = 1.0
    wagegapscale: float = 1.0
    rtoexpscale: float = 1.0
    # --- grid sizes and tolerances (53-57) --------------------------------------
    nI: int = 16
    nT: int = 4
    nY: int = 3
    nE: int = 30
    nZ: int = 2
    bisectTol: float = 1e-5
    VFtol: float = 0.1
    maxViter: int = 100
    maxHSiter: int = 100
    # --- stochastic environment (66-82) ------------------------------------------
    piz: np.ndarray = field(default_factory=lambda: np.array([[0.9, 0.1], [0.3, 0.7]]))
    lossW: np.ndarray = field(default_factory=lambda: np.array([0.04, 0.07]))
    lamH: np.ndarray = field(default_factory=lambda: np.array([
        [[0.96, 0.0, 0.04], [0.3, 0.60, 0.1], [0.0, 0.95, 0.05]],
        [[0.93, 0.0, 0.07], [0.25, 0.60, 0.15], [0.0, 0.88, 0.12]]]))
    findW: np.ndarray = field(default_factory=lambda: np.array([0.8, 0.5]))
    wageH: np.ndarray = field(default_factory=lambda: np.array([0.8, 0.8, 0.78]))
    BCwage: np.ndarray = field(default_factory=lambda: np.array([1.0, 0.85]))
    ym: np.ndarray = field(default_factory=lambda: np.array([1.0, 0.75, 0.0]))
    tsize: int = 4
    age_years: tuple = (14, 14, 8)   # expected years in each age group (line 82)
    egrid_max: float = 2.0           # line 112
    egrid_curv: float = 1.5          # line 112

    # ------------------------------------------------------------------ derived
    def __post_init__(self):
        self.piz = np.asarray(self.piz, float)
        self.lossW = np.asarray(self.lossW, float)
        self.lamH = np.asarray(self.lamH, float)
        self.findW = np.asarray(self.findW, float)
        self.wageH = np.asarray(self.wageH, float)
        self.BCwage = np.asarray(self.BCwage, float)
        self.ym = np.asarray(self.ym, float)

    # experiment scalings exactly as in lines 60-61 of the solver
    @property
    def tau_wf_eff(self) -> float:
        return 1.0 - self.wagegapscale * (1.0 - self.tau_wf)

    @property
    def gam_e_eff(self) -> float:
        return self.gam_e * self.rtoexpscale

    @property
    def egrid(self) -> np.ndarray:
        g = np.linspace(0.01, 1.0, self.nE)
        return self.egrid_max * g ** self.egrid_curv

    @property
    def ftW(self) -> np.ndarray:
        f = np.full(self.nI, self.ftWbase, float)
        f[: self.nI // 2] *= self.ftWlow_scale
        return f

    def piT_matrix(self) -> np.ndarray:
        """3x4 ageing transition matrix as written on line 82."""
        y1, y2, y3 = self.age_years
        t = self.tsize
        return np.array([
            [1 - 1 / (y1 * t), 1 / (y1 * t), 0, 0],
            [0, 1 - 1 / (y2 * t), 1 / (y2 * t), 0],
            [0, 0, 1 - 1 / (y3 * t), 1 / (y3 * t)]])

    def piT_used(self, opts: Options) -> np.ndarray:
        """Probability of ageing out of group it (it = 0,1,2) as used in the solver.

        MATLAB indexes the 3x4 matrix linearly: piT(1)=piT(1,1), piT(2)=piT(2,1),
        piT(3)=piT(3,1) -> [0.982, 0, 0].  The intended value is piT(it, it+1).
        """
        P = self.piT_matrix()
        if opts.aging_linear_index:
            return P.flatten(order="F")[: self.nT - 1]
        return np.array([P[i, i + 1] for i in range(self.nT - 1)])

    def kap_matrix(self, opts: Options) -> np.ndarray:
        """Fixed cost of work by (age, type), shape (nT, nI). Lines 92-107."""
        nT, nI = self.nT, self.nI
        kapbar = np.array([self.kaplow, 0.0])
        kap = np.full((nT, nI), self.kapbase, float)
        kap[1, :] *= self.kap_mid_scale
        kap[2, :] *= self.kap_old_scale
        kayscale = np.asarray(self.kayscale, float)
        for iw in range(2):
            for ik in range(2):
                for j in range(4):
                    c = iw * 8 + ik * 4 + j
                    if opts.young_cost_all_ages:
                        # MATLAB line 102-103: kap(:, c) = kapscale*kap(:, c)/kayscale(j)
                        kap[:, c] = self.kapscale * kap[:, c] / kayscale[j]
                    else:
                        kap[0, c] = self.kapscale * kap[0, c] / kayscale[j]
                        kap[1:, c] = self.kapscale * kap[1:, c]
            if opts.kapbar_indexing_bug:
                # MATLAB line 106: kap(:, (iw-1)*8+1 : (iw-1)*8+5) += kapbar(iw)
                kap[:, iw * 8: iw * 8 + 5] += kapbar[iw]
            else:
                for ik in range(2):
                    kap[:, iw * 8 + ik * 4: iw * 8 + ik * 4 + 4] += kapbar[ik]
        return kap

    # --------------------------------------------------------------- loaders
    @classmethod
    def from_paras_mat(cls, path: str) -> "Params":
        """Build a Params object from a saved MATLAB paras.mat."""
        import scipy.io as sio
        m = sio.loadmat(path, squeeze_me=True)
        kw = {}
        for name in ["beta", "crra", "eta", "mu", "phi_c", "r", "alpha_e", "delta_e",
                     "psi", "xi", "lhome", "z_h", "ybar_h", "alpha_h", "nu", "nu_h",
                     "ftWbase", "kapbase", "kaplow", "kapscale", "wagegapscale",
                     "rtoexpscale", "bisectTol", "VFtol", "tsize"]:
            if name in m:
                kw[name] = float(m[name])
        for name in ["nI", "nT", "nY", "nE", "nZ", "maxViter", "maxHSiter", "tsize"]:
            if name in m:
                kw[name] = int(m[name])
        for name in ["piz", "lossW", "lamH", "findW", "wageH", "BCwage", "ym"]:
            kw[name] = np.asarray(m[name], float)
        kw["kayscale"] = tuple(np.asarray(m["kayscale"], float))
        # paras.mat stores the *scaled* tau_wf and gam_e; undo the scaling so
        # that tau_wf_eff / gam_e_eff reproduce the stored values.
        wgs = kw.get("wagegapscale", 1.0)
        rts = kw.get("rtoexpscale", 1.0)
        kw["tau_wf"] = 1.0 - (1.0 - float(m["tau_wf"])) / wgs
        kw["gam_e"] = float(m["gam_e"]) / rts
        p = cls(**kw)
        egrid = np.asarray(m["egrid"], float)
        assert np.allclose(p.egrid, egrid), "egrid formula differs from saved grid"
        ftW = np.asarray(m["ftW"], float)
        assert np.allclose(p.ftW, ftW), "ftW differs from saved"
        return p

    def replace(self, **kw) -> "Params":
        return replace(self, **kw)
