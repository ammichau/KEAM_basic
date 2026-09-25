"""Render output/departure_impacts.csv as a markdown table into DEPARTURES.md."""
import os, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
csv = os.path.join(HERE, "..", "output", "departure_impacts.csv")
doc = os.path.join(HERE, "..", "..", "..", "DEPARTURES.md")
df = pd.read_csv(csv, index_col=0)
keep = ["E/pop (all ages)", "hours | E", "quit rate/qtr (exp)", "quit rate/qtr (rec)",
        "E->nonE/qtr (exp)", "dE/pop rec-exp (pts)", "dHours rec-exp (%)", "wife inc share (exp)",
        "share PT", "share Lifecycle", "share Career", "share NiLF", "share undefined (of all)",
        "share person-qtrs e > grid max", "mean gQ (working ages)", "max|dV| vs faithful"]
short = {"E/pop (all ages)": "E/pop", "hours | E": "hours|E", "quit rate/qtr (exp)": "quit/q exp",
         "quit rate/qtr (rec)": "quit/q rec", "E->nonE/qtr (exp)": "E→nonE/q", "dE/pop rec-exp (pts)": "ΔE rec-exp (pts)",
         "dHours rec-exp (%)": "Δhours rec (%)", "wife inc share (exp)": "wife share", "share PT": "PT",
         "share Lifecycle": "Lifecycle", "share Career": "Career", "share NiLF": "NiLF",
         "share undefined (of all)": "undefined", "share person-qtrs e > grid max": "e>grid",
         "mean gQ (working ages)": "mean gQ", "max|dV| vs faithful": "max|ΔV|"}
t = df.loc[keep].T
lines = ["| scenario | " + " | ".join(short[k] for k in keep) + " |",
         "|---|" + "---|" * len(keep)]
for name, row in t.iterrows():
    cells = []
    for k in keep:
        v = row[k]
        cells.append("" if pd.isna(v) else (f"{float(v):.3f}" if abs(float(v)) < 100 else f"{float(v):.1f}"))
    lines.append(f"| {name} | " + " | ".join(cells) + " |")
table = "\n".join(lines)
table += ("\n\nColumns: E/pop = employment rate over all simulated person-quarters; quit/q = quarterly quit "
          "rate (true expansion / recession quarters, 1973-2019 window); E→nonE/q = quarterly employment-to-"
          "non-employment rate in expansions; ΔE = employment-rate change in recessions in percentage points; "
          "career shares with the undefined group dropped (R2) and its share of all women; e>grid = share of "
          "person-quarters with experience above the solver grid; mean gQ = share of solved states where "
          "quitting is optimal; max|ΔV| = largest change in the value function relative to the MATLAB solution.\n")
src = open(doc).read()
if "<!-- IMPACT_TABLE -->" in src:
    src = src.replace("<!-- IMPACT_TABLE -->", "<!-- IMPACT_TABLE -->\n" + table)
else:
    src = src[: src.index("| scenario |")] + table
open(doc, "w").write(src)
print(table)
