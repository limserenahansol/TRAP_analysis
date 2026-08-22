"""
Fig · Cluster 3 region-selection heatmap — same style as Fig 4 (cluster 4).

Rows ranked by Post+Reinstatement Active−Passive (within-phase z).
Highlights:
  - solid black: interest / discussed (COAa, AVP, SUM)
  - dashed green: behaviour-match (A>P Post/Rein, A<P Withdrawal)

Out: region_selection_PPT/Fig_cluster3_region_selection.png
     forMark/Fig_cluster3_region_selection.png
"""
from pathlib import Path
import shutil
import trap_region_selection_ppt as T

OUT = T.OUT
FM = Path(r"C:\Users\hsollim\behavior_task\TRAP_analysis_sync") / \
    "TRAP_OUTPUT_calculated_mm3" / "forMark"

# Cluster-3 regions of current interest (discussed recently)
SEL_C3 = ["COAa", "AVP", "SUM"]

C3_NOTE = (
    "Cluster 3 (n=30, hypothalamic / midline-thalamic family)\n"
    "rows ranked by Post+Reinstatement Active−Passive · "
    "boxed = interest (COAa, AVP, SUM)\n"
    "• COAa: Passive > Active at both Post & Withdrawal (raw density).\n"
    "• AVP, SUM: behaviour-match pattern (Active Post peak + Rein rebound; Passive WD peak).\n"
    "• solid black = interest · dashed green = behaviour-match (A>P craving, A<P Withdrawal)."
)

ap = T.load_ap_long()
mat = T.ap_delta_matrix(ap)
bm = T.behavior_match_set(mat, 3)

n_c3 = (mat["cluster"] == 3).sum()
out_png = OUT / "Fig_cluster3_region_selection.png"
d = T._selection_heatmap(
    mat, 3, set(SEL_C3), out_png,
    f"Fig · Cluster 3 (n={n_c3}, hypothalamic / midline-thalamic family): region selection\n"
    "rows ranked by Post+Reinstatement Active−Passive · boxed = interest (COAa, AVP, SUM)",
    criteria_note=C3_NOTE,
    behavior_match=bm,
)
csv_path = OUT / "Fig_cluster3_ranked.csv"
d.to_csv(csv_path, index=False)

# also copy to forMark
shutil.copy(out_png, FM / out_png.name)
shutil.copy(csv_path, FM / csv_path.name)

print("Saved:", out_png)
print("Also:", FM / out_png.name)
print(f"Cluster 3 n={n_c3}")
print(f"Interest boxed: {SEL_C3}")
print(f"Behaviour-match (dashed): {sorted(bm)}")
print("\nRanked (top 15 by Post+Rein score):")
print(d[["acronym", "During", "Post", "Withdrawal", "Reinstatement", "score"]].head(15).to_string(index=False))
