# Cross-project: opioid behavior ↔ TRAP ↔ Genelist (WMB)

**Status:** Updated 2026-09-09 with the ORBm/BMAp Xenium add-on and collaborator gene lists.

## Shared ideas to document here (not in three places)

- **Cohort / mouse IDs:** which animals appear in more than one pipeline; exclusions; Active–Passive pairing rules.
- **Phase / day alignment:** what “day” means in behavior vs TRAP vs gene lists.
- **Outputs that bridge projects:** e.g. region lists, phenotype scores, PR-related summaries — file names + repo links.

## Repositories

- Behavior / opioid MATLAB: [opioidaddiction-matlab](https://github.com/limserenahansol/opioidaddiction-matlab) — see [../opioidaddiction-matlab.md](../opioidaddiction-matlab.md)
- TRAP: [TRAP_analysis](https://github.com/limserenahansol/TRAP_analysis) — see [../TRAP_analysis.md](../TRAP_analysis.md)
- Genelist WMB: [Genelist_analysis_WMB](https://github.com/limserenahansol/Genelist_analysis_WMB) — see [../Genelist_analysis_WMB.md](../Genelist_analysis_WMB.md)

## Outputs that bridge projects (2026-09-09)

- Xenium shared panel for ORBm + BMAp (TRAP2 × Ai14; Active vs yoked Passive): [MSGS111 workbook](https://github.com/limserenahansol/Genelist_analysis_WMB/blob/main/v3/outputs/FINAL_Xenium_panel_ORBm_BMAp_for_MSGS111.xlsx).
- Jesse PL-ILA-ORB morphine DEGs vs that panel: [Jesse_ORB_vs_Xenium_panel.xlsx](https://github.com/limserenahansol/Genelist_analysis_WMB/blob/main/v3/outputs/Jesse_ORB_vs_Xenium_panel.xlsx).
- Dan GSE283418 98-gene spatial panel vs BMAp: [GSE283418_vs_BMAp_panel.xlsx](https://github.com/limserenahansol/Genelist_analysis_WMB/blob/main/v3/outputs/GSE283418_vs_BMAp_panel.xlsx).
- PI deck: [ORBm_BMAp_Jesse_Dan_panel_decision.pptx](https://github.com/limserenahansol/Genelist_analysis_WMB/blob/main/v3/outputs/ORBm_BMAp_Jesse_Dan_panel_decision.pptx) (also copied under `docs/xenium_ORBm_BMAp_2026-09/` in this repo).

## Open questions

- Keep all 14 Dan extras (114 custom) or trim to `Lamb3` + 012 VGLUT1 border genes to stay near the 100-slot cap?
- Add Jesse morphine-state genes (`Per2`, `Pcsk1`, `Per1`, `Rxfp1`) if the Xenium readout should include 5-day morphine dependence, not only cell type + TRAP?
