from pathlib import Path
import re
import pandas as pd

ROOT = Path(r"C:\Users\hsollim\behavior_task\TRAP_analysis_sync")
wb = pd.read_excel(ROOT / "Hansol Lim 561 cell counts + densitynew.xlsx", sheet_name="All Samples")
wb["base"] = wb["acronym"].astype(str).str.replace(r"-[LR]$", "", regex=True)
man = pd.read_csv(ROOT / "TRAP_sample_manifest.csv")
man = man[man["include"] == 1]
meta = man.drop_duplicates("mouse_id").set_index("mouse_id")[["delivery", "phase"]].to_dict("index")
dens = [c for c in wb.columns if c.endswith("density (cells/mm^3)")]
count = [c for c in wb.columns if c.endswith("count")]

for reg in ["ORBm", "BMAp", "COAa"]:
    rsub = wb[wb["base"] == reg]
    rows = []
    for dcol, ccol in zip(dens, count):
        m = re.match(r"(HaLi_\d+_\d+)", dcol)
        if not m or m.group(1) not in meta:
            continue
        mid = m.group(1)
        d = float(pd.to_numeric(rsub[dcol], errors="coerce").mean())
        c = float(pd.to_numeric(rsub[ccol], errors="coerce").sum())
        if pd.isna(d):
            continue
        md = meta[mid]
        cohort = re.search(r"HaLi_(\d+)_", mid).group(1)
        rows.append({"mouse_id": mid, "cohort": cohort, "group": md["delivery"],
                     "phase": md["phase"], "density": d, "count": c})
    df = pd.DataFrame(rows).sort_values("density", ascending=False)
    hi, lo = df.iloc[0], df.iloc[-1]
    print("=" * 64)
    print(reg)
    print(f"  HIGHEST: {hi.mouse_id} | {hi.group} {hi.phase} | dens={hi.density:.1f} | count={hi['count']:.0f}")
    print(f"  LOWEST:  {lo.mouse_id} | {lo.group} {lo.phase} | dens={lo.density:.1f} | count={lo['count']:.0f}")
    df102 = df[df.cohort == "102125"]
    print(f"  102125 hi: {df102.iloc[0].mouse_id} | {df102.iloc[0].group} {df102.iloc[0].phase} | {df102.iloc[0].density:.1f}")
    print(f"  102125 lo: {df102.iloc[-1].mouse_id} | {df102.iloc[-1].group} {df102.iloc[-1].phase} | {df102.iloc[-1].density:.1f}")
    # Passive hi for COAa especially
    p = df[df.group == "Passive"]
    a = df[df.group == "Active"]
    print(f"  Active hi:  {a.iloc[0].mouse_id} | {a.iloc[0].phase} | {a.iloc[0].density:.1f}")
    print(f"  Passive hi: {p.iloc[0].mouse_id} | {p.iloc[0].phase} | {p.iloc[0].density:.1f}")
