# Genelist_analysis_WMB (summary)

**Canonical code:** https://github.com/limserenahansol/Genelist_analysis_WMB

## What belongs on this page

- Gene list / WMB analysis scope, inputs (atlas, cell sets), and outputs **as narrative** + links to scripts or notebooks in that repo.
- Relationship to anatomy or TRAP-derived region sets (if you join on region or ontology IDs, document the join rule here once).

## Current Xenium add-on (ORBm + BMAp) — 2026-09-09

Shared order workbook: [`FINAL_Xenium_panel_ORBm_BMAp_for_MSGS111.xlsx`](https://github.com/limserenahansol/Genelist_analysis_WMB/blob/main/v3/outputs/FINAL_Xenium_panel_ORBm_BMAp_for_MSGS111.xlsx)

- **144** curated genes (44 free on Xenium Mouse Brain v1 + **100** custom). All 20 cell types (12 ORBm + 8 BMAp) remain separable.
- **+14 Dan / GSE283418 genes** appended (block `12_GSE283418_added`) → **158** shared / **114** custom if all are kept.
- Jesse morphine DEGs were compared, **not** appended. After Allen ORBm re-score (106,122 cells): add `Rxfp1` (free), `Per2`, `Pcsk1`, `Per1`, `Chrm1`, `Grm8`. Optional: `Gpr26`, `Mas1`, `Camk2g`. Skip `Mchr1`.
- 4-slide deck: [`Add_genes_Jesse_Allen_Dan_4slides.pptx`](https://github.com/limserenahansol/Genelist_analysis_WMB/blob/main/v3/outputs/Add_genes_Jesse_Allen_Dan_4slides.pptx)

Local copies of the deck and figures: [`docs/xenium_ORBm_BMAp_2026-09/`](../../docs/xenium_ORBm_BMAp_2026-09/).

| Collaborator | Data | Decision |
|---|---|---|
| Daniel Berg / Scherrer | GSE283418 Resolve 98-gene smFISH | 38 already on panel; 14 added; rest CEA/glia/weak |
| Jesse Niehaus | PL-ILA-ORB morphine DEGs + ORB-enriched GPCRs | Core IEG/opioid already on panel; morphine-state layer still open |

Allen check: `Cckbr` and `Hcrtr2` are truly abundant in our ORBm WMB-10X tables. Jesse Glut #1 GPCRs (`Mas1`, `Gpr68`, …) were not in our 40-GPCR Allen pull.

## What does not belong here

- Full gene tables — link to repository artifacts or supplementary files.

## Cross-links

- [Cross-project bridge](cross/opioid_TRAP_genelist.md)
- [Index](index.md)
