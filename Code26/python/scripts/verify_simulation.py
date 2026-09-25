"""Simulate every saved MATLAB solution in faithful mode and compare with the stored
spreadsheets (Output/<name>/SimulStats*.xls and Code26/CrossCohort.xlsx).

NOTE: the stored outputs were produced with NiLFdef = 0.3 (the commented-out
alternative on line 59 of SimplerMod_May17_sim.m), not the committed 0.21.
"""
import sys, os, time, warnings
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import numpy as np, pandas as pd
from keam import Options, simulate, stats, SimConfig
from keam.matio import load_matlab_solution
warnings.simplefilter("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
CODE26 = os.path.join(HERE, "..", "..")
ROOT = os.path.join(CODE26, "..")
nilfdef = float(sys.argv[1]) if len(sys.argv) > 1 else 0.3

cc = pd.read_excel(os.path.join(CODE26, "CrossCohort.xlsx"), sheet_name="Careers", header=None)
cc_careers = {1940: cc.iloc[2:6, 1].values, 1950: cc.iloc[2:6, 3].values, 1960: cc.iloc[2:6, 5].values,
              1970: cc.iloc[2:6, 7].values, 1980: cc.iloc[2:6, 9].values}   # PT, Cycle, Career, NiLF
xls = {"Baseline": os.path.join(ROOT, "Output", "Baseline", "SimulStats.xls"),
       "RoE_incr": os.path.join(ROOT, "Output", "RoE_incr", "SimulStats2.xls"),
       "Wgap_dcr": os.path.join(ROOT, "Output", "Wgap_dcr", "SimulStats2.xls"),
       "Kap_dcr": os.path.join(ROOT, "Output", "Kap_dcr", "SimulStats2.xls")}
cohort_of = {"Baseline": 1940, "Cohort1950": 1950, "Cohort1960": 1960, "Cohort1970": 1970, "Cohort1980": 1980}

for sub in ["Baseline", "RoE_incr", "Wgap_dcr", "Kap_dcr", "Cohort1950", "Cohort1960", "Cohort1970", "Cohort1980"]:
    sol = load_matlab_solution(os.path.join(CODE26, "Solution", sub))
    cfg = SimConfig(NiLFdef=nilfdef)
    r = simulate(sol.params, sol, Options.faithful(), cfg)
    ca = stats.careers(r); cy = stats.cycle(r); cs = stats.cross_section(r)
    print(f"===== {sub}  (tau_wf={sol.params.tau_wf_eff:.3f}, gam_e={sol.params.gam_e_eff:.4f}, kapscale={sol.params.kapscale})")
    mine = np.array([ca["PT"], ca["Cycle"], ca["Career"], ca["NiLF"]])
    print(f"  careers python : PT={mine[0]:.6f} Cycle={mine[1]:.6f} Career={mine[2]:.6f} NiLF={mine[3]:.6f}  (undefined dropped: {ca['Undefined']:.3f})")
    if sub in cohort_of:
        t = cc_careers[cohort_of[sub]].astype(float)
        print(f"  CrossCohort.xlsx: PT={t[0]:.6f} Cycle={t[1]:.6f} Career={t[2]:.6f} NiLF={t[3]:.6f}   max|diff|={np.abs(mine-t).max():.2e}")
    if sub in xls:
        x = pd.ExcelFile(xls[sub])
        car = x.parse("Careers", header=None).iloc[1, :4].values.astype(float)
        print(f"  SimulStats.xls  : PT={car[0]:.6f} Cycle={car[1]:.6f} Career={car[2]:.6f} NiLF={car[3]:.6f}   max|diff|={np.abs(mine-car).max():.2e}")
        cyc = x.parse("Cycle", header=None)
        rows = {str(cyc.iloc[k, 0]): (cyc.iloc[k, 1], cyc.iloc[k, 2]) for k in range(1, 18)}
        # xlswrite labels: 'EU' row holds m_EN and 'EN' row holds m_EU (see DEPARTURES D-R4)
        pairs = [("H_Emp", "H_Emp"), ("Emp", "Emp"), ("EE", "EE"), ("EU", "EN"), ("EN", "EU"), ("EnonE", "EnonE"),
                 ("NonEE", "NonEE"), ("Hours", "Hours"), ("Quit", "Quit"), ("Wage", "Wage"),
                 ("Household Income", "IncHH"), ("Wife Income", "IncW"), ("Wife Share of Inc", "Wshare")]
        worst = 0.0
        for lab, key in pairs:
            m_rec, m_exp = float(rows[lab][0]), float(rows[lab][1])
            d = max(abs(cy[key]["rec"] - m_rec), abs(cy[key]["exp"] - m_exp))
            worst = max(worst, d)
            print(f"    {lab:18s} rec {cy[key]['rec']:.6f} vs {m_rec:.6f} | exp {cy[key]['exp']:.6f} vs {m_exp:.6f}")
        print(f"  cycle stats max|diff| = {worst:.2e}")
        xs = x.parse("CrossSection", header=None)
        print(f"    hours young/mid/old python {cs['mHours'].round(6)} vs xls {xs.iloc[2,1:4].values}")
        print(f"    quit rate by age python {cs['Qrate'].round(6)} vs xls {xs.iloc[11,1:4].values}")
