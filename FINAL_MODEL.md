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
  100 types with equal weights.
* **States:** experience e ∈ [0, 2], assets a ≥ 0, husband x_m ∈ {E, R, U}, aggregate
  Z ∈ {expansion, recession}, employment status, and an iid transitory cost shock κ_T
  (5-point discrete normal with sd σ_κ) drawn at the start of each month.
* **Preferences:** u(c) − µ h^{1+η}/(1+η) − κ_τ − κ_T when employed, u(c) when not;
  u = c^{1−γ}/(1−γ), γ = 2, η = 1.4, β = 0.99. κ_τ = κ̄ κ_m at ages 25-39 and κ̄ afterwards.
* **Wage:** w = φ(Z) τ_w ω (1 + γ_e e^ξ) with γ_e = 0.5, ξ = 0.8, φ(rec) = 0.88.
* **Experience:** e' = min(2, (1−δ)e + θ e h^ψ) while employed, e' = (1−δ)e otherwise;
  δ = 0.005, θ = 0.025, ψ = 0.66 (slides p.28, monthly). Initial e uniform on [0.2, 1.0].
* **Home production:** f(ω)(1−h)^{ν_h} with f = ȳ_h + z_h ω^{α_h}, z_h = 0.45, α_h = 0.21,
  ν_h = 0.65; the non-employed use f(ω)(1−s)^{ν_h}.
* **Search and job loss:** job finding s^ν λ_f(Z), ν = 0.5, λ_f(rec) = 0.85 λ_f(exp);
  exogenous job loss λ_u(Z). A non-employed woman counts as unemployed if s ≥ s̄.
* **Husband:** states E, R (re-employed with scar), U. Monthly transitions: E→U 1.35%
  (2.4% in recessions); U→R 35% (28%); R→E 2.5% (scar lasts 3.3 years on average); R→U
  2.5 times the E→U rate. Income y_H(τ) × (1, 0.85, 0.30) by state, y_H = (0.89, 1.0, 0.94)
  by the wife's age group, times φ_H(Z) with φ_H(rec) = 0.88 (assumed equal to the wife's).
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
| σ_κ | transitory cost shock | monthly quit rates 3.4% / 2.8% |
| τ_w | wage penalty | wage gap 0.71 |
| λ_f | job-finding efficiency | employment 0.62, U rate 5% |
| λ_u(exp), λ_u(rec) | exogenous loss | E→nonE 5.0% / 4.8% |
| ȳ_h | home production level | wife share of income 22.5% |
| s̄ | unemployment definition | U rate 5% |

Externally set: β, γ, η, γ_e, ξ, φ, δ, θ, ψ, z_h, α_h, ν_h, ν, husband process, ageing, pension.
