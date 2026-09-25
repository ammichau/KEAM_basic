# Departures between the Code26 computation and the model in the slides/paper

Reference documents: slides `KEAM_Klein.pdf` (v. 4.28.2025, "the winner") and
paper draft `EM_DemogBCtrends.pdf` (Feb 2023). Code: `Code26/SimplerMod_May17_splines.m`
(solver), `Code26/SimplerMod_May17_sim.m` (simulator), `Code26/ModelSim_22.do` (Stata).
Line numbers refer to those files. Each item has an ID; solver items S*, simulator
items M*, reporting items R*, structural items A*, parameter items P*. Where a
Python flag exists it is named (`keam.params.Options`); `True` = MATLAB behaviour.

**Validation status.** The Python translation (`Code26/python`) in faithful mode
reproduces all eight saved solutions (`Code26/Solution/*`) to machine precision and
the stored `Output/Baseline/SimulStats.xls` and `Code26/CrossCohort.xlsx` career shares
to 1e-6. So the list below describes what the code that produced the slides' numbers
actually does. Note that the committed solver has a typo on line 82 (`14se*tsize`) that
MATLAB cannot parse; the saved solutions were produced before it was introduced.

The quantitative impact of switching each item off individually is in the table at
the end (`Code26/python/output/departure_impacts.csv`).

---

## A. Structure: what is computed versus what is written

| ID | Text (slides unless noted) | Code | Where |
|---|---|---|---|
| A1 | Households hold assets `a` with `R = 0` and a borrowing limit (slides p.13, 15, 47). | No savings at all; consumption equals current income every period. README calls this the "stripped down model". | solver 199, 312 |
| A2 | Wage `w = φ(Z) τ_w ω (1 + γ_e e^ξ)`: experience multiplies the fixed type (slides p.25, 28). | `w = τ_w (ω + γ_e e^ξ)`: additive. For the low type (ω = 0.5) the return to experience relative to the fixed component is twice the slides' form. The paper's mechanism section (p.32-33) uses the additive form. | `wage.m` |
| A3 | Husband unemployed receives UI, `φ_m(u) < φ_m(r) < 1` (paper p.13-14). | Husband income is **zero** when unemployed: `ym = [1, 0.75, 0]`. The "recently unemployed" state is the re-employed-with-scar state (60% persistence per quarter, i.e. 2.5 quarters, vs "3.3-4 years" in the paper). | solver 69-70, 76 |
| A4 | Extra cost of work `κ_m` applies during child-raising years only; "Life-cycle women" quit young and return (slides p.24-25, 29). | The young-age multiplier `1/kayscale` is applied to **all ages** (`kap(:, …)`, and `kap(2,:)*1`, `kap(3,:)*1`). The cost matrix has identical rows: there is no life-cycle variation in the cost of work in the solved model. The simulator's panel export (line 1163-1167) and the Stata code (`Kap_y`) assume the young-only structure. Flag `young_cost_all_ages`. | solver 94-104 |
| A5 | Model period is a month (paper p.21, β = 0.99). Slides silent. | Quarterly (`tsize = 4`, β = 0.97). Annual discounting is consistent (0.886 vs 0.885), but all rates (job loss 4/7%, husband transitions, search efficiency) are quarterly. | solver 30, 81 |
| A6 | Retiree solves a consumption-savings problem (paper p.16). | Retirement value is 0 for every state and, because of S1, never entered in expectation. | solver 134 |
| A7 | Age groups 25-39, 40-54, 55-64 (15/15/10 years). | Solver expected durations 14/14/8 years; simulator deterministic 14/14/11 years (ages 25-63). | solver 82; sim 54, 75-83 |
| A8 | Aggregate state Markov. | Solver beliefs `piz = [.9 .1; .3 .7]` (25% of time in recession, mean spell 3.3 quarters); simulation feeds NBER dates 1955-2019 (15.7% of quarters in the 1973-2019 window). Beliefs and realised process differ. | solver 66; sim 230-267 |
| A9 | Job offer leads to `V^E` (slides p.47) or to `V = max(V^E, V^N)` (paper p.16). | `max`: an offer can be declined at no cost. | solver 320, 336 |
| A10 | Permanent cost `κ̄ ~` truncated normal on [0, 0.075], `κ_m ~ U[1, 2.27]` (slides p.29). | Two fixed offsets {-0.04, 0} (negative costs possible, minimum total cost 0.035) and four young multipliers {1, 1.11, 1.18, 2.0}. The simulator draws a truncated normal on [-0.04, 0] and **linearly interpolates the two solved policies**. | solver 51-52, 92-107; sim 107-128, 173-182 |
| A11 | Permanent productivity ω drawn from a distribution (slides p.13, 25). | Two solved wage types {0.5, 1}; the simulator draws a truncated normal on [0.5, 1] (10-point grid) and linearly interpolates policies between the two solved types. The 0/1 quit rule becomes a quit probability equal to the interpolation weight (a quit lottery). No simulated woman's policy is the solution of her own problem. | sim 86-103, 173-182 |
| A12 | Unemployed vs out of the labour force. | A non-employed woman is "unemployed" if search intensity > 0.2, otherwise NiLF. In the baseline 0.03% of the population is unemployed; the slides' U-rates (p.36: 7.2% … 4.9%) cannot come from this code. | sim 60, 453-457 |
| A13 | Experience `e ∈ [1, ē]`, initial e "one-tenth to one-half of maximum" (paper p.27 fn). | Grid on [0.002, 2]; initial e uniform on grid points 12-18 (0.48-0.91). 3-7% of simulated person-quarters have e above the grid maximum (flat policy extrapolation, wage keeps growing). | solver 111-112; sim 273 |

## P. Parameter values (slides win)

| Parameter | Slides | Paper Table 2/3 | Code (Baseline paras.mat) |
|---|---|---|---|
| ψ_e (hours curvature in experience) | 0.66 | 0.8 | **0.8** |
| θ_e | 0.025 | 0.026 | **0.026** |
| ξ_e | 0.8 | 0.8 | **0.85** |
| δ_e | 0.005 | 0.005 | 0.005 |
| γ_e 1940 (cohorts 1950-80) | 0.50 (0.55, 0.58, 0.68, 0.69) | 0.3 | 0.5 (0.5475, 0.575, 0.675, 0.685) |
| φ(recession) | 0.88 | 0.88 (text p.22: 0.85) | **0.85** (`BCwage`) |
| ν (search) | 0.5 | 0.5 | **0.4** |
| ν_h (home production hours) | 0.65 | 0.65 | **0.5** |
| η | 1.4 | 1.4 | 1.4 |
| γ (CRRA) | 2 | 2 | 2 |
| µ_h | — | 1.0 | 0.5 with `phi_c = 0.5` on consumption utility. Dividing through by `phi_c`: µ = 1.0 but every κ is effectively doubled (0.07-0.30). |
| β | — | 0.99 monthly | 0.97 quarterly |
| ȳ_h, z_h, α_h | form only | 0.11, 0.45, 0.21 | 0.1, 0.45, 0.2 |
| κ̄ range | [0, 0.075] | κ̄ = 0.075 | offsets [-0.04, 0] plus 0.075/kayscale |
| κ_m range | [1, 2.27] | 0.17 (sic) | {1, 1.11, 1.18, 2.0}, all ages (A4) |
| job finding efficiency in recession | 15% lower | — | `findW = [0.8, 0.5]`: **37.5% lower** |
| exogenous job loss, women | 1.7x in recession | 3%/5% (men 4%/7%) | 4%/7% (1.75x) |
| husband scar | — | 0.78, 3.4 yrs, 2.5x loss | 0.75, 2.5 quarters, 2.5x |
| husband income by wife's age | — | young 0.89, old 0.94 | [0.8, 0.8, 0.78] |
| τ_w 1940 | (moment 0.71) | τ_1940 = 1.0 | 0.8 |
| τ_w by cohort | gap **closes** (0.71→0.77) | — | τ_w **falls**: 0.80, 0.786, 0.762, 0.720, 0.718 (`wagegapscale` 1.07-1.41 > 1). The measured wage-gap moment closes only through returns to experience. Only the `Wgap_dcr` experiment (τ_w = 0.84) moves the way the text describes. |
| kapscale by cohort | 1.0, 1.01, 0.88, 1.0, 1.04 | — | 1.0, 1.01, 0.84, 1.0, 1.037 (`SimplerMod_Cohorts_Experiment.m` has 0.80/1.2/1.095 for 1950 as the live values and the 1960 file says 0.84 whereas the slides say 0.88) |

## S. Solver implementation (`SimplerMod_May17_splines.m`)

| ID | Flag | Where | What happens | Consequence |
|---|---|---|---|---|
| S1 | `aging_linear_index` | 82, 243-246, 259-262, 319-321, 335-338 | `piT` is a 3x4 matrix but is indexed `piT(it)` (linear, column-major): `piT(1) = 0.982`, `piT(2) = piT(3) = 0`. | Young women age into the middle group with probability **0.982 per quarter**; middle-aged and old women **never** age or retire (infinite horizon). The intended values are 1/56, 1/56, 1/32. All life-cycle incentives (returns to experience near retirement, young-age quits) are distorted. |
| S2 | `stale_h_nonemployed` | 297 | In the non-employed block `ee = exp2(egrid(ie), h)` uses whatever `h` currently holds: the last employed hours for the first state, then the previous state's converged **search** intensity (the same variable `h` is reused). | Experience of the non-employed does not depreciate in the continuation value (and rises when the previous state's search exceeds 0.13). The simulator, by contrast, depreciates it (`exp2(e, 0)`). The solver's quit decision therefore does not see the cost of losing experience; the mechanism in slides p.20/32 is switched off. |
| S3 | `roundup_continuation` | 205-212, 255-258, 298-305, 331-334 | `iee = find(egrid > ee, 1, 'first')` and `spline(egrid, V, egrid(iee))` evaluate V **at the grid point above** e' (a spline at a knot is the knot value). The derivative uses the segment above e'. | Continuation experience is rounded up every period; no interpolation despite the file name. Employed and non-employed continuations are rounded differently (S2), so the employed continuation can sit at a lower experience than the non-employed one. |
| S4 | `derivative_at_current_state` | 218-235 | `dVedh(1,izz,iyy) = (V(i,it,iee+1,iy,iz) - V(i,it,iee,iy,iz))/…`: the marginal value of experience is taken at today's `(iy, iz)`, not at `(iyy, izz)` inside the expectation. | The hours FOC ignores that the value of experience differs across future husband/aggregate states: the precautionary-hours channel (slides p.17, 43) is absent from the intensive margin. |
| S5 | `alpha_h_typo` | 229 | `dVUedh(1) *= alpha_h*psi*h^(psi-1)*egrid(ie)` uses `alpha_h = 0.2` (home-production curvature) instead of `alpha_e = 0.026`. | The experience gain in the job-loss branch is overstated 7.7x, inflating the dynamic return to hours by roughly 30%. |
| S6 | `kapbar_indexing_bug` | 106 | `kap(:, (iw-1)*8+1 : (iw-1)*8+5) += kapbar(iw)`: the fixed-cost offset is added by **wage type** over five columns instead of by fixed-cost type over four. | High-wage types (9-16) all get offset 0, low-wage types 1-5 get -0.04, 6-8 get 0. The simulator's interpolation across the "fixed cost" dimension (which assumes types 1-4/9-12 low, 5-8/13-16 high) is fed policies with the wrong costs. |
| S7 | `young_cost_all_ages` | 95-96, 102-103 | See A4. | No life-cycle cost variation. |
| S8 | `loose_vf_tolerance` | 134, 175-181, 286, 360, 377-381 | (i) `VFtol = 0.1` in value units (V ≈ -13); (ii) the criterion is `abs(max(V0 - V))`, the least-negative change, not the sup norm; (iii) `V = V0` and `VU0 = VU1` assign whole arrays that were initialised to zero, so the initial guess is wiped for type 1 at the younger ages; (iv) each type starts from the previous type's value function. | Type 1 stops after 46-47 Bellman updates; **types 2-16 perform a single Bellman update** each (see `n_iter` in the verification script). The 16 types' policies are barely differentiated by their own parameters. Converged value functions differ from the stored ones by up to 0.43 (3% of V). |
| S9 | — | 82 | `1/(14se*tsize)` is a syntax error. | The committed solver does not run; the `.asv` autosave is identical. Fixed in the Python translation. |

## M. Simulation (`SimplerMod_May17_sim.m`)

| ID | Flag | Where | What happens | Consequence |
|---|---|---|---|---|
| M1 | `reversed_e_weights` | 421, 433, 449, 496, 515, 531 | Linear interpolation in experience puts weight `(e - e_lo)/Δ` on the **lower** grid point. | Policies are mirrored inside each grid cell; the quit threshold shifts by up to one cell. |
| M2 | `reversed_type_weights_gQ` | 181-182 | The quit rule is interpolated across the fixed-cost and wage types with reversed weights (a woman at the low end gets the high type's rule) and the interpolated 0/1 rule is used as a quit probability. | Low-wage women mostly use high-wage quit rules and vice versa; quits become lotteries. `gH`, `gS` use `interp1` and are correct. |
| M3 | `same_seed_kappa_draws` | 117, 137 | `rng(222)` is used for both the fixed-cost draw and the young-cost draw. | The two cost components are perfectly rank-correlated instead of independent. |
| M4 | `no_bcwage_in_sim` | 414, 476 | The wife's simulated wage is `wage(ftW_i, e)` without `BCwage(iz)`, while the solver's budget includes it. | Wife's income, household income and the "wife share" statistics ignore the 15% recession wage cut that her decisions were based on; the counter-cyclicality of her income share is overstated. |
| M5 | `all_employed_at_birth` | 417 | `EmpI_i(ii,tt) = 1` overrides the 8% `Nw_init` draw; the 8% keep `UnempI = 1`. | Those women are counted as both employed and unemployed in the birth quarter and generate a spurious U→E flow the next quarter. |
| M6 | `lagged_z_in_sim` | 479-490 | Husband transitions are drawn with last period's `(iy, iz)`; job finding uses `findW(z_t)`. The solver uses `λ_H(z')` and `findW(z')`. | Timing of the aggregate state in the husband and job-finding processes differs between beliefs and simulation by one quarter. |
| M7 | `quit_before_loss` | 500-505 | Quit is checked first, job loss only `elseif`. Solver: loss with probability λ_u, quit in the `1-λ_u` branch. | P(quit) = Q instead of (1-λ_u)Q; measured quits ~5% too high relative to the model. |
| M8 | (always replicated) | 494, 449/531, 547 | `ee` is recomputed only when the woman enters the period employed. The search policy of the non-employed is looked up at the experience **index of her last employed quarter** but weighted with her current (depreciated) `e`, which extrapolates and can turn negative. MATLAB then computes `(negative)^0.4` as a complex number and `>` compares real parts (87 person-quarters in the baseline). | Search intensity of the long-term non-employed is evaluated at the wrong experience. |
| M9 | — | 274-348 | Types (100 women) and all idiosyncratic shock paths are drawn once and copied to every birth cohort. | 64 cohorts share 100 life histories; effective sample size is 100, and cohort differences are pure calendar effects. |
| M10 | — | 365-371, 1181 | The value-at-birth and panel `ValueFn` columns use types 1 and 2 (both low-wage) with reversed weights. | The welfare/"Vfun" columns are not the model's value function. |
| M11 | — | 1184 | Panel consumption utility adds `Wage*Hours` to `Inc_hh`, which already contains it. | Wife's earnings double-counted in the exported `UtilC`. |
| M12 | — | 474 | `AgeY = age0 + ceil(it/4)`: quarters 2-4 of the first year are age 26. | Cohort/age bins in Stata off by one in the first year. |

## R. Statistics and reporting

| ID | Flag | Where | What happens | Consequence |
|---|---|---|---|---|
| R1 | `expansion_includes_recessions` | 957-962 | `ExpI(1,t) = 1` is set in recessions too (should be 0). | Every "expansion" column (`Exp_means`, `CrossCohort.xlsx`, slides p.30/37) is an **unconditional** mean. Recession-minus-expansion differences are understated by the factor (1-p) ≈ 0.84. |
| R2 | `drop_undefined_careers` | 908-933 | Women with employment rate between the NiLF cut-off and 0.8 who are not life-cycle are "undefined" and excluded from the denominator. | 31% of baseline women (21-26% in the cohorts) are dropped from the reported career shares (slides p.30, 33, 36). |
| R3 | `phantom_last_cohort` | 273-311 vs 402 | Initial draws are written for all 64 cohorts but only 63 are simulated. | Untouched birth-quarter entries of the 64th cohort inflate employment counts in one quarter (1e-4 on the means). |
| R4 | — | 1103-1112 | Header labels and values are misaligned in `xlswrite`: the 'EU' row holds the E→N rate and 'EN' holds E→U; the 'Search' … 'Wife Share' columns are shifted by the two `BirthV` entries. | Spreadsheet columns must be re-labelled before use. |
| R5 | — | `CrossCohort.xlsx` Cycle sheet | Quarterly flow rates are divided by 4 and compared to monthly CPS rates (slides p.30: "quit E-NonE expansion 3.5%" is 14.06% per quarter ÷ 4). | The correct monthly equivalent of 14.06% per quarter is 4.9%; the E→nonE rate 18.4%/quarter is 6.5% monthly, not 4.6%. |
| R6 | `SimConfig.NiLFdef` | 59 | Committed value 0.21; the stored baseline output and all `CrossCohort.xlsx` career shares reproduce only with the commented alternative **0.3** (`cyclefactor = 0.5`). | The slides' career shares use NiLF = employment rate < 0.3. |
| R7 | — | 892-920 | Career taxonomy: employment rate > 0.8 → career (share of quarters with hours > 0.39 above 0.7) or part-time; young rate < 0.5 × middle rate → life-cycle; else rate < NiLFdef → NiLF; quarters 12-156 (ages 28-63). Slides p.24: annual hours ≥ 1500 / 400-1500 / < 400 over ages 25-54, life-cycle = NiLF → career. | Definitions differ from the data definitions they are compared with. |
| R8 | — | 935 | `wgap = mWage_E(2)/wageH(1)` prints 1.25 in the baseline. The slides' 0.71 (p.30, 35) and the employment rates 62%/67% (p.35, 37; the cycle sheet says 51%, cross-section 47-54%) are not produced by any formula in the code. | Targets in the slides cannot be traced to the code. |
| R9 | — | slides p.37 | The 1940 row's quit change (-33.0%) is the 1980 value in `CrossCohort.xlsx` (1940: -26.7%); wife-share change +2.5 vs +2.1 in the sheet. | Transcription. |
| R10 | — | `ModelSim_22.do` 161-185 | `hh_hours = hours + h_emp*40` mixes a time share with weekly hours; `collapse (sum)` by recession then `dbc = sum_rec - sum_exp` is dominated by the number of quarters. `Agg_relto_Base.csv` shows the **baseline's change relative to itself as -84%**. | Slides p.38 ("-115.2%", "+114.8%") are artifacts of this block. |
| R11 | — | `Output/{RoE_incr,Wgap_dcr,Kap_dcr}/SimulStats2.xls` | These spreadsheets do not correspond to the saved `Code26/Solution/*` solutions (the shell references a `_betteronBC` solver that is not in the repository). Re-simulating the saved comparative-static solutions gives different career shifts from slides p.33-34 (e.g. RoE↑: PT -8.7, LC -13.3, Career +31.3, NiLF -9.4 points vs -5, -10, +33, -18). | The comparative statics in the slides are not reproducible from the repository. |

## Q. Questions for you (answers change the "final model")

1. Assets: keep the no-savings model, or add the asset state from the slides?
2. Cost of work: should `κ_m` apply only to ages 25-39 (text) — this changes the life-cycle mechanism and the meaning of the 16 types.
3. Wage form: additive `τ(ω + γ e^ξ)` (code, paper mechanism section) or multiplicative `τ ω (1 + γ e^ξ)` (slides)?
4. Cohort calibration: is the fall in τ_w across cohorts intended (the *raw* gender wage penalty widening while the measured gap closes)?
5. Period length: quarter or month? All transition rates need re-mapping if month.
6. Types: replace the 16-type solve plus interpolation with solving each simulated woman's own type (the corrected Python code can loop over types cheaply)?
7. Career taxonomy: keep the employment-rate version, or implement the hours-based definitions of slides p.24?
8. Husband's UI replacement rate when unemployed (currently 0).

## Quantitative impact of each item (baseline parameters, `NiLFdef = 0.3`)

Each row switches **one** flag to the textbook behaviour with everything else as in
MATLAB; the last rows switch all solver flags, then everything. Generated by
`Code26/python/scripts/departure_impacts.py`.

<!-- IMPACT_TABLE -->
