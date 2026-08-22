# Colleague pack — TRAP final-PDF figures only

**Audience:** colleague building an AI analysis benchmark against Hansol’s TRAP story.  
**Scope:** only figures that appear in the final deck PDF (not the full Step-13 output dump).

| Item | Location |
|------|----------|
| **This pack** | `Research_Projects/02_TRAP_wholebrain/colleague_AI_benchmark_pack/` |
| **NEW — full explain PPT/PDF** | [`TRAP_colleague_explain_figure_code.pptx`](TRAP_colleague_explain_figure_code.pptx) · [`TRAP_colleague_explain_figure_code.pdf`](TRAP_colleague_explain_figure_code.pdf) |
| **Original final story PDF** | [`ppt_source/11TRAP_data_ORBm_BMAp_COA_for_next_step33.pdf`](ppt_source/11TRAP_data_ORBm_BMAp_COA_for_next_step33.pdf) |
| **How the PDF was built** | [`HOW_FINAL_PDF_WAS_MADE.md`](HOW_FINAL_PDF_WAS_MADE.md) |
| **Figure ↔ slide ↔ code** | [`FINAL_PDF_FIGURE_MAP.md`](FINAL_PDF_FIGURE_MAP.md) |
| **Biological aim + workflow** | [`01_BIOLOGICAL_AIM_WORKFLOW_FIGURES_CODE.md`](01_BIOLOGICAL_AIM_WORKFLOW_FIGURES_CODE.md) |
| **GitHub code** | https://github.com/limserenahansol/TRAP_analysis |

## Read order

1. **Open `TRAP_colleague_explain_figure_code.pdf`** (or `.pptx`) — figure + meaning + code on each slide  
2. `HOW_FINAL_PDF_WAS_MADE.md` — how the original story deck was assembled  
3. `FINAL_PDF_FIGURE_MAP.md` — file map for pack figures  
4. `figures/` + `code/` for offline assets  

Rebuild the explain deck anytime: `python build_colleague_explain_pptx.py`

## What is *not* in this pack

- Raw Step-13 universal-cluster PNGs that never entered the final deck  
- Large histology panels (original PDF slides 14–16 ORBm, 19–21 BMAp, 29–30 COAa) — open the original final PDF for those
