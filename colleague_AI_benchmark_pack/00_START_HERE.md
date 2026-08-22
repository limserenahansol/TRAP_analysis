# Colleague pack — TRAP whole-brain analysis (AI benchmarking)

**Audience:** colleague building an **AI analysis benchmark** against Hansol’s TRAP pipeline.  
**Purpose:** navigate biological aim → workflow → each figure → code → meaning, without digging the full MATLAB repo first.

| Item | Location |
|------|----------|
| **This pack** | `C:\Users\hsollim\Research_Projects\02_TRAP_wholebrain\colleague_AI_benchmark_pack\` |
| **Canonical MATLAB code** | https://github.com/limserenahansol/TRAP_analysis |
| **Local code sync** | `C:\Users\hsollim\behavior_task\TRAP_analysis_sync\` |
| **Final story PPT (PDF)** | [`ppt_source/11TRAP_data_ORBm_BMAp_COA_for_next_step33.pdf`](ppt_source/11TRAP_data_ORBm_BMAp_COA_for_next_step33.pdf) |
| **Master guide (read this)** | [`01_BIOLOGICAL_AIM_WORKFLOW_FIGURES_CODE.md`](01_BIOLOGICAL_AIM_WORKFLOW_FIGURES_CODE.md) |
| **Figure gallery** | [`figures/`](figures/) |
| **Code copies (offline)** | [`code/`](code/) |
| **Neuro-Agent task package** | `behavior_task\TRAP_analysis_sync\neuro_agent_trap_submission\` |

## Read order (10–20 min)

1. **Biological aim** → section 1 of the master guide  
2. **Funnel workflow** → section 2 (hundreds of regions → clusters → ORBm/BMAp/COAa)  
3. **Figure ↔ code ↔ meaning** → section 3 (matches PPT Fig 1, 3, 4, 5 + Finding B)  
4. **How to re-run** → section 4  
5. **What an AI agent should reproduce** → section 5 (rubric-style)

## Folder map

```text
colleague_AI_benchmark_pack/
  00_START_HERE.md                          ← you are here
  01_BIOLOGICAL_AIM_WORKFLOW_FIGURES_CODE.md ← full navigation doc
  02_CODE_INDEX_GITHUB.md                   ← GitHub paths only
  figures/                                  ← PNGs + CSVs from MATLAB Step 13
  code/                                     ← key .m files (mirrors GitHub layout)
  ppt_source/                               ← final PPT PDF
```

## Density / cohort used for these figures

- **Variant:** `calculated_mm3` (cells / sample volume in mm³)  
- **Mask:** `forebrain_no_bs` + `z_within_phase`  
- **Mice:** 20 included (both cohorts); 1 excluded (`HaLi_020326_11`)  
- **Output root on disk:**  
  `TRAP_analysis_sync\TRAP_OUTPUT_calculated_mm3\13_universal_cluster_PCA_density\forebrain_no_bs\z_within_phase\`
