# TRAP_analysis — Whole-brain TRAP density pipeline (MATLAB + Python)

Analysis pipeline for **TRAP whole-brain imaging**: density (cells/mm³) per Allen region, **Active vs Passive**, and **phase**. BRANCH-style stats, clustering, flip-direction analysis, **per-region statistics (Step 6)**, and follow-ups (Step 7).

| Component | Location |
|-----------|----------|
| **MATLAB** (Steps 1–7) | This repo root — `RUN_PIPELINE_ALL.m` |
| **Python** (Steps 1–7) | **`trap_pipeline_py/`** — `python run_pipeline.py --csv … --manifest …` |

**Run order (MATLAB):** [`WORKFLOW.md`](WORKFLOW.md) Steps 1–5 + **Step 6–7** via `RUN_PIPELINE_ALL.m` (or run Step 6/7 scripts individually).  
**Core question (Step 6):** For **each** brain region, do Active and Passive differ in TRAP density **within Reinstatement** and **within Withdrawal**? Default test: **Wilcoxon rank-sum** (`trap_config.m` → `phase_AP_test = 'ranksum'`).  
**Step 7:** Phase-delta screening + exploratory cross-phase folders (see `trap_pipeline_matlab/README.md` logic in root `README` section below).

**Where code vs outputs live:** [`FOLDERS_GUIDE.md`](FOLDERS_GUIDE.md)  
**Add cohorts / mice:** [`WHEN_YOU_ADD_MICE_EN_KR.md`](WHEN_YOU_ADD_MICE_EN_KR.md)

---

## Before you run (MATLAB)

1. **`init_TRAP_pipeline`** once per MATLAB session.
2. **`trap_config.m`** — paths, `phase_AP_test` (`ranksum` | `welch`), FDR vs raw p.
3. **`TRAP_cohort_CSVs.txt`** — one density CSV per line.
4. **`TRAP_sample_manifest.csv`** — column_name, delivery, phase, include.

## Python pipeline

```bash
cd trap_pipeline_py
pip install -r requirements.txt
python run_pipeline.py --csv /path/to/density.csv --manifest /path/to/TRAP_sample_manifest.csv --out TRAP_PYTHON_OUTPUT
```

See **`trap_pipeline_py/README.md`**.

---

## Folder layout (MATLAB)

| Folder | Step |
|--------|------|
| **`shared/`** | Helpers |
| **`Step_01_…` – `Step_05_…`** | BRANCH, correlations, clustering, flip, utilities |
| **`Step_06_phase_AP_contrasts/`** | Regionwise Active vs Passive + Step 7 helpers |
| **`trap_pipeline_py/`** | Python mirror (optional) |

**Outputs:** `TRAP_OUTPUT/06_regionwise_Active_vs_Passive/`, `TRAP_OUTPUT/07_followup/`, etc.

---

## Requirements

- MATLAB (Statistics Toolbox)
- Python 3.10+ recommended for `trap_pipeline_py` (see `trap_pipeline_py/requirements.txt`)
- Optional: UMAP (`run_umap`)
