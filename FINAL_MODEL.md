# Final model specification (`Code26/python/keam/final`)

This is the model the paper describes, with the decisions taken on 2026-09-25: monthly
period, assets with a borrowing constraint, a 30% replacement rate for the husband's
unemployment income, and a compensated wage-gap experiment. Everything below is implemented;
the table at the end lists what is calibrated.

## Environment

* **Period:** one month. Women enter at 25, pass through age groups 25-39, 40-54, 55-64
  (stochastic ageing with the correct expected durations in the solver, deterministic ages in
  the simulation) and retire at 65.
* **Type** θ = (ω, κ̄, κ_m), fixed at entry and solved for explicitly: ω is log-normal with
  sd 0.37 (the wage fixed-effect dispersion) on 5 points normalised to mean 1; κ̄ is a
  truncated normal on [0, κ̄_max] on 5 points; κ_m is uniform on [1, κ_m,max] on 4 points.
  100 types with equal weights. **Option (2026-09-27):** a permanent home-productivity type, `n_zh`
  equiprobable multipliers 1 ± `zh_spread` on z_h (types × n_zh; `make_types4`), motivated by the
  never-working share being one wage-type cell in every calibration (RESULTS.md 1b). On the coarse grid
  at the version-7b parameters a spread of 0.3 moves the never-working share from 0.30 to 0.22 and the
  career shares to 0.28 / 0.29 / 0.21 / 0.22 (data 0.31 / 0.28 / 0.19 / 0.22) without recalibration;
  On the 100-type grid at the corrected KPR calibration (v7c) the type does the opposite: a spread of
  0.15 / 0.30 / 0.45 raises the never-working share from +12% to +15% / +22% / +25% and polarises the
  shares (career up, life-cycle and part-time down; `output/eval_v7c_zh*.out`), and the least-squares
  recalibration with the type (v8c) did not improve on v7c. Not adopted; the coarse-grid result was a
  27-type artefact.
* **Added-worker effect (2026-09-27).** The model's wives insure against the husband's risk by not quitting,
  not by entering: the monthly labor-force entry ratio (husband E→U versus stays E; data 1.60, Guner,
  Kulikova and Valladares-Esteban 2025) is 1.01 (v7c) and 1.11 (v4e), and neither removing assets, lengthening
  the husband's spells nor a persistent cost shock (ratios 0.66-0.91) raises it (`output/awe_levers_v7c.md`,
  `output/awe_rho_scan_v7c.json`). Entry is decided by the wage type, not by the household's state.
* **States:** experience e ∈ [0, 2], assets a ≥ 0, husband x_m ∈ {E, R, U}, aggregate
  Z ∈ {expansion, recession}, employment status, and a cost-of-work shock κ_T (discrete
  normal with sd σ_κ) realised at the start of each month, before the quit decision. In the
  iid version (`rho_kT = 0`, 5 nodes) the shock is redrawn every month; in the persistent
  version (`rho_kT = ρ_κ > 0`, 3 nodes) it keeps its value with probability ρ_κ and is
  otherwise redrawn, so a bad spell lasts 1/(1−ρ_κ) months on average and the current node
  is a state variable of both value functions (`keam/final/solve.py`). The persistent version
  is under calibration (`scripts/calibrate_ls.py`, `output/final_calib_rho_*`); the results in
  `RESULTS.md` use the iid version unless stated.
* **Preferences:** u(c) − µ h^{1+η}/(1+η) − κ_τ − κ_T when employed, u(c) when not;
  u = c^{1−γ}/(1−γ), γ = 2, η = 1.4, β = 0.99. κ_τ = κ̄ κ_m at ages 25-39 and κ̄ afterwards.
* **Wage:** w = φ(Z) τ_w ω (1 + γ_e e^ξ) with γ_e = 0.5, ξ = 0.8, φ(rec) = 0.88.
* **Experience:** e' = min(2, (1−δ)e + θ e h^ψ) while employed, e' = (1−δ)e otherwise;
  δ = 0.005, θ = 0.025, ψ = 0.66 (slides p.28, monthly). Initial e uniform on [0.2, 1.0].
* **Home production:** f(ω)(1−h)^{ν_h} with f = ȳ_h + z_h ω^{α_h}; the non-employed use
  f(ω)(1−s)^{ν_h}. **Child care:** at ages 25-39 home productivity is multiplied by
  m_c ≥ 1 (`home_young_mult`), the paper's "opportunity cost of home production around
  child bearing" (p.12); together with the utility-cost multiplier κ_m this generates the
  life-cycle women. ȳ_h, z_h, α_h, ν_h and m_c are calibrated (see below).
* **Search and job loss:** job finding s^ν λ_f(Z), ν = 0.5, λ_f(rec) = 0.85 λ_f(exp);
  exogenous job loss λ_u(Z). A non-employed woman counts as unemployed if s ≥ s̄.
* **Husband:** states E, R (re-employed with scar), U. Monthly transitions: E→U 1.35%
  (2.4% in recessions); U→R 35% (28%); R→E 2.5% (scar lasts 3.3 years on average); R→U
  2.5 times the E→U rate. Income y_H(τ) × (1, 0.85, 0.30) by state, y_H = (0.89, 1.0, 0.94)
  by the wife's age group, times φ_H(Z) with φ_H(rec) = 0.88 (assumed equal to the wife's).
  Variant under evaluation (author request, 2026-09-26): `ui_rec_mult` multiplies the U-state share
  in recessions (0.5: replacement 15% instead of 30%, standing in for longer unemployment spells).
  At the adopted parameters this alone lowers the recession quit rate from 2.50% to 2.27% and the
  recession employment drop from 1.69 to 1.11 points (`output/eval_ui_direct.out`); the
  recalibrated version is `output/final_calib_ui_full.*`. **Version 3** (`output/final_calib_v3_full.*`,
  objective 0.113, the best fit so far) combines the recession UI cut with a 10% fall in the wife's
  job-finding efficiency in recessions (λ_f(rec) = 0.90 λ_f(exp) instead of 0.85) and a calibrated
  recession rise in her own job-loss rate; in it precautionary labor supply accounts for 29% of the
  recession quit drop and job hoarding for 25% (`output/channels_v3.md`).
* **Job-finding cyclicality target (author's data, 2026-09-26).** The standard deviation of the
  cyclical component of the UE rate is 0.0765 for men and 0.0686 for women. Under the two-state
  aggregate process the model counterpart is |log(UE_rec / UE_exp)| √(π_exp π_rec), computed from the
  simulated UE rate of unemployed searchers (`keam/final/moments.py`, "sd log UE (women)"; the
  husband's analogue reproduces the men's number with the assumed 0.35 / 0.28 finding rates). The
  measured rate falls less than the efficiency λ_f because search rises in recessions, so matching
  0.0686 needs an efficiency fall of roughly 25%. Versions: **v4** (target added, fall calibrated;
  the least-squares run left it at 15%), **v4b** (fall fixed at 25%), **v5** (v4b's assumptions with
  log utility, γ = 1, for a balanced growth path; `u(c) = log c` in `keam/final/solve.py`).
* **Epstein-Zin option (2026-09-27).** `ez_rra` > 1 with γ = 1 gives Epstein-Zin preferences with unit
  intertemporal elasticity and relative risk aversion `ez_rra`: next-period risks (husband's state,
  aggregate state, cost shock, own job loss, job finding, ageing) enter through the certainty
  equivalent −(1/θ) log E exp(−θ V) with θ = (ez_rra − 1)(1 − β), the log form of
  V = c^(1−β) [E V'^(1−α)]^(β/(1−α)); mortality acts as discounting. Nests expected utility at
  `ez_rra` = 1 (checked to 1e-9). Finding: at the log-utility calibration (v5b) the precautionary
  share of the recession quit drop stays at 7-9% for risk aversion 1, 2, 5 and 10
  (`output/channels_v5b_rra{2,5,10}.md`), so the channel is governed by the intertemporal
  elasticity (period-utility curvature, γ = 2 gives 0.5), not by risk aversion. Balanced growth
  with a precautionary labor-supply channel therefore needs non-separable consumption-hours
  preferences of the King-Plosser-Rebelo class.
* **King-Plosser-Rebelo option (2026-09-27).** `kpr = True`: the flow utility is U(x) with the
  composite x = log c − µ h^{1+η}/(1+η) − κ_τ and U(x) = exp((1−γ) x)/(1−γ) (U(x) = x at γ = 1), i.e.
  u = [c·exp(−v(h) − κ)]^{1−γ}/(1−γ): consistent with balanced growth for any γ, nesting the separable
  log model exactly at γ = 1 (checked). The costs µ, κ̄, κ_m are then in log-consumption units, as
  in the log model; the transitory shock κ_T remains an additive shock to the value at the quit
  decision. At the log-utility calibration (v5b) with γ = 2 imposed, the precautionary share rises
  from 7% to 22% (hoarding 50%; `output/channels_v5b_kpr2.md`); version 7 recalibrates this
  specification (`output/final_calib_v7_full.*`).
  With `kT_mult = True` (2026-09-27) the transitory shock is treated like the fixed costs: employment in a
  month with shock κ_T is taxed proportionally, u = [c·ψ(h)·exp(−(κ + κ_T))]^{1−γ}/(1−γ), so that the
  participation comparison V^E − V^N is scale-invariant (balanced growth). The shock is drawn at the start
  of the month, before the quit/accept decision, and the employed woman's value carries today's node; with
  an iid shock the continuation and V^N do not depend on it. At γ = 1 this is identical to the additive
  shock (`output/verify_kT_mult.md`). Version 9 (`output/final_calib_v9_full.*`) is calibrated with it.
* **Assets:** c + a' = income + a, a' ≥ 0, gross return 1 (net rate zero). The asset grid has
  20 points on [0, 15] (about ten months of household income) with more points near zero;
  a' is chosen on the grid, policies are interpolated bilinearly in (e, a) in the simulation.
* **Retirement:** absorbing, income 0.5 × middle-age husband income, monthly death hazard
  1/240; V_R(a) solved by value iteration. Wife's own pension is ignored.
* **Aggregate state:** monthly Markov chain with P(exp→rec) = 0.015, P(rec→exp) = 0.09
  (14% of months in recession, 11-month spells). The simulation can use either a drawn path
  (calibration) or the NBER dates 1955-2019 (cohort narrative).

## Bellman equations (per type)

V^E(e,a,y,Z) = max_{h,a'} u(c) − v(h) − κ_τ + β E[(1−λ_u(Z)) V(e',a',y',Z') + λ_u(Z) V^N(e',a',y',Z')]
with c = w h + f(1−h)^{ν_h} + y_m + a − a', e' = g(e,h);

V^N(e,a,y,Z) = max_{s,a'} u(c) + β E[(1 − s^ν λ_f(Z)) V^N(e',a',y',Z') + s^ν λ_f(Z) V(e',a',y',Z')]
with c = f(1−s)^{ν_h} + y_m + a − a', e' = (1−δ)e;

V = E_{κ_T} max{V^E − κ_T, V^N}; the expectation E also covers ageing (and retirement from
the last group). Hours are on a 20-point grid, search on 21 points, savings on the 20-point
asset grid; continuation values are linearly interpolated in e'. Modified policy iteration
(one maximisation, 25 evaluation steps) to a sup-norm tolerance of 1e-5.

## Simulation and moments

Annual entry cohorts of 60 women with independent type and shock draws, 90 cohorts on a
drawn aggregate path; moments over the fully populated calendar months. Flows are monthly:
a quit is a woman working last month who chooses non-employment this month; E→nonE adds
exogenous losses. Careers use annual hours (4,000-hour endowment) over ages 25-54 exactly as
on slide 24: career ≥ 1,500, part-time 400-1,500, NiLF < 400, life-cycle = ≥ 1,500 at 40-54
and < 600 at 25-39 (supersedes). The wage gap is the ratio of the employed wife's full-time
equivalent monthly earnings to the employed husband's monthly income.

## Experiments

* Returns to experience: γ_e scaled.
* Compensated wage gap: τ_w scaled by g and y_H scaled by 1 − s_w (g − 1)/(1 − s_w), where
  s_w is the baseline wife share of household income, so that expected household income at
  baseline behaviour is unchanged.
* Cost of work: κ̄_max scaled.
* Each is sized by bisection to reproduce the 1970s cohort employment rate; cohort accounting
  uses the observed τ_w and γ_e paths and the residual cost.
* Mechanism counterfactuals: acyclical husband transitions, acyclical job finding, no recession
  wage cut, acyclical own job loss.

## Calibrated parameters and targets

| parameter | role | target |
|---|---|---|
| µ | hours disutility | hours of employed 0.40 |
| κ̄_max, κ_m,max | cost of work levels | career shares (LC 31, PT 28, career 19, NiLF 22) |
| σ_κ | cost-of-work shock | monthly quit rates 3.4% / 2.8% |
| ρ_κ (persistent version) | persistence of the cost shock | employment among ever-working women (E/pop jointly with the never-working share) |
| τ_w | wage penalty | wage gap 0.71 |
| λ_f | job-finding efficiency | employment 0.62 |
| λ_u(exp), λ_u(rec) | exogenous loss | E→nonE 5.0% / 4.8% |
| ȳ_h, z_h, α_h, ν_h | home production level, its slope in ω, curvature in hours | NiLF and part-time shares, hours |
| m_c (child-care multiplier) | home productivity at 25-39 | life-cycle share |

Externally set: β, γ, η, γ_e, ξ, φ, δ, θ, ψ, ν, s̄ = 0.25, husband process, ageing, pension.
Not targeted, reported: unemployment rate, wife's share of income, consumption response to the
husband's job loss, assets.
