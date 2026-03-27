# TRAP_analysis — Whole-brain TRAP density pipeline (MATLAB)

Analysis pipeline for **TRAP whole-brain imaging**: density (cells/mm³) per Allen region, **Active vs Passive**, and **phase**. BRANCH-style stats, PCA/UMAP, k-means, correlations, flip-direction analysis, regionwise statistics, and an optional **five-phase timeline** (Step 10).

**Quick links:** [Roadmap (schematic)](#roadmap-schematic) · [New data checklist](#when-you-have-new-data) · [Steps 1–10](#workflow-steps-1–10-matlab) · [What each figure type means](#what-the-figures-mean) · [Methods EN/KR](#methods--statistics-en--kr--추가-문서)

**Publication schematics (Steps 1–3, PNG):** [`docs/steps_1_2_3_publication_schematics/`](docs/steps_1_2_3_publication_schematics/) — raster figures plus [`PHASE_LABELS_AND_Z.txt`](docs/steps_1_2_3_publication_schematics/PHASE_LABELS_AND_Z.txt) (canonical phase labels at load vs optional within-phase z-scoring).

---

## Roadmap (schematic)

High-level order when you run **`init_TRAP_pipeline`** then **`RUN_PIPELINE_ALL`**:

```mermaid
flowchart TB
  subgraph IN["Inputs"]
    CFG["trap_config.m"]
    CSV["Density CSV(s)<br/>Allen id × samples"]
    CL["TRAP_cohort_CSVs.txt"]
    MAN["TRAP_sample_manifest.csv"]
  end
  subgraph S1["1–5 Exploration"]
    A1["Step 1 BRANCH<br/>global stats"]
    A2["Step 2a Clustering sweep<br/>K, silhouette, PCA"]
    A3["Step 3 v2<br/>regions × phases"]
    A4["Step 4 Flip<br/>Rein+With patterns"]
    A5["Step 5 Utilities<br/>name lists"]
  end
  subgraph S2["6–9 Active vs Passive & follow-ups"]
    B6["Step 6 Phase AP<br/>per-region tests"]
    B6b["Step 6b Phase-Δ screen"]
    B7["Step 7 Directional scenarios"]
    B8["Step 8 Within-group<br/>Rein vs With"]
    B9["Step 9 Same as 6–8<br/>forebrain mask"]
  end
  subgraph S3["10 Timeline"]
    B10["Step 10 Five-phase<br/>timeline AP + within"]
  end
  IN --> A1 --> A2 --> A3 --> A4 --> A5
  A5 --> B6 --> B6b --> B7 --> B8 --> B9 --> B10
```

How **rows and columns** flow into statistics (Steps 6–10 share this loader logic):

```mermaid
flowchart LR
  CSV["Cohort CSVs"] --> POOL["L+R mean<br/>per region"]
  MAN["Manifest"] --> POOL
  POOL --> NORM["Phase labels<br/>e.g. Reexposure→Rein"]
  NORM --> DROP["Drop Exclude<br/>if configured"]
  DROP --> MASK["Region mask<br/>match Step 3"]
  MASK --> Z["Optional within-phase<br/>z per region"]
  Z --> TEST["Wilcoxon ranksum<br/>where applicable"]
```

---

## When you have new data

Do these **before** or **instead of** a full rerun, depending on what changed.

### 1. Add or replace the density CSV(s)

- Copy each file into the **repo root** (next to `trap_config.m`) **or** put them anywhere and reference **absolute paths** in `TRAP_cohort_CSVs.txt`.
- Each CSV must have Allen columns such as **`id`**, **`name`**, **`acronym`**, **`parent_structure_id`**, **`depth`**, plus **one column per mouse/sample** (density).
- **Every cohort CSV must use the same `id` list** as cohort 1 (row order can differ; the code aligns by `id`).

### 2. Files you may need to change

| Situation | What to edit |
|-----------|----------------|
| New filename or location | **`trap_config.m`** → `C.csvPath` (single-file mode), **and/or** **`TRAP_cohort_CSVs.txt`** (one path per line). |
| Second (third, …) cohort file | Add a **new line** to **`TRAP_cohort_CSVs.txt`**. Line index = **`cohort_id`** in the manifest. |
| New mice / new columns | **`TRAP_sample_manifest.csv`** — one row per sample column you want in the analysis. |
| New phase labels (e.g. Baseline, During) | Manifest column **`phase`**. If the spelling is nonstandard, extend **`shared/trap_normalize_manifest_phase.m`**. |
| Five-phase timeline (Step 10) | Manifest must include samples for **`phase5_phases`** in **`trap_config.m`** (default: Baseline, During, Post, Withdrawal, Reinstatement) and a **baseline** phase for deltas. Details: **[`STEP10_NEW_DATA_FIVE_PHASE_WORKFLOW.md`](STEP10_NEW_DATA_FIVE_PHASE_WORKFLOW.md)**. |
| Same mice for Step 3 and Steps 6–9 | **`trap_config.m`** → `C.v2_sample_source = 'manifest'`. See **`STEP_CONSISTENCY_3_vs_6_8.md`**. |

### 3. Manifest columns (required)

For **each** sample column in the CSV(s):

- **`cohort_id`** — which line of `TRAP_cohort_CSVs.txt` (1, 2, …), or `1` if only one CSV.
- **`column_name`** — **exact** header string in that cohort’s CSV.
- **`delivery`** — e.g. `Active`, `Passive`.
- **`phase`** — e.g. `Withdrawal`, `Reinstatement`, `Baseline`, … (normalized; see `trap_normalize_manifest_phase.m`).
- **`include`** — `1` = use, `0` = skip.

Optional: **`mouse_id`** (for your notes; the pipeline matches data by **`column_name`** + **`cohort_id`**).

### 4. Run MATLAB

```matlab
cd('...\TRAP_analysis')   % folder that contains init_TRAP_pipeline.m
init_TRAP_pipeline
RUN_PIPELINE_ALL
```

- **Faster test:** in **`trap_config.m`**, set **`C.runMode = 'quick'`** (fewer permutations / bootstrap off).
- **Plots:** prefer the **MATLAB desktop** app; some **`matlab -batch`** environments skip PNG export.
- **Step 10 only** (after the rest has been run at least once with a valid manifest):

```matlab
trap_run_phase5_timeline_analysis
% Optional forebrain-heavy mask:
% trap_run_phase5_timeline_analysis(struct('phase_AP_row_filter_fn', @trap_AP_filter_forebrain_exclude_fiber_wm))
```

More detail: **[`WHEN_YOU_ADD_MICE_EN_KR.md`](WHEN_YOU_ADD_MICE_EN_KR.md)** (EN+KR), **[`DENSITY_CSV_SETUP.txt`](DENSITY_CSV_SETUP.txt)**.

---

## Workflow: Steps 1–10 (MATLAB)

Executed in order by **`RUN_PIPELINE_ALL.m`** (after **`init_TRAP_pipeline`**).

| Step | What it does | Main outputs under `TRAP_OUTPUT/` (or config paths) |
|------|----------------|-----------------------------------------------------|
| **1** | **BRANCH** — whole-brain summaries of Active vs Passive (global / tree / embedding views depending on script). | **`01_BRANCH_tables_and_figures/`** — CSVs + **`figures_described/`** |
| **2a** | **Clustering sweep** — scan **K**, **silhouette**, **stability**, **sample PCA**. | **`02_clustering_sweep/figures_described/`** |
| **3** | **Region clustering v2** — phase-aware clustering; builds **`TRAP_downstream_input.mat`** for Step 4. | **`03_region_clustering_v2/`** (+ RepRegions CSVs, figures) |
| **4** | **Flip / downstream** — regions whose **Reinstatement** and **Withdrawal** Active−Passive pattern meets **Conditions A/B/C** (joint across phases; different question from Step 6). | **`04_flip_downstream/figures_described/`** |
| **5** | **Utilities** — export region name lists (e.g. depth 5–6). | Next to v2 / config paths |
| **6** | **Phase-specific Active vs Passive** — per region, **within Reinstatement** and **within Withdrawal** separately; **Wilcoxon rank-sum**; FDR trees / volcanoes / bar plots. | **`06_phase_ActivePassive_FDR/`** (see `trap_config.m` → `phase_AP_root`) |
| **6b** | **Phase-delta screening** — how **\|Δ_Rein − Δ_With\|** ranks regions (exploratory). | Under phase AP / follow-up roots (see console + folder READMEs) |
| **7** | **Directional scenarios** — separate folders for **scenario** contrasts (volcano + directional bar-style summaries). | **`07_directional_AP_scenarios/`** |
| **8** | **Within-group Rein vs Withdrawal** — for **Active-only** and **Passive-only** mice: regional **mean difference** Rein−With + **Wilcoxon**. | **`08_within_group_Rein_vs_Withdrawal_delta/`** |
| **9** | **Same analyses as 6–8** on a **forebrain-focused** region set (excludes brainstem + cerebellum per Step 9 mask). | **`09_forebrain_no_brainstem_cerebellum/`** |
| **10** | **Five-phase timeline** — within delivery: means per phase, **Δ vs baseline**, heatmaps / line plots; cross-group: **Active vs Passive within each phase**. | **`10_five_phase_timeline/`** (`phase5_timeline_root`) |

**Step 4 vs Step 6:** Step 6 is “**one phase at a time**” (A vs P in Rein **or** in With). Step 4 is “**both phases together**” for a **pattern** (e.g. opposite signs). **Top-region lists need not match** — that is expected.

---

## What the figures mean

Most PNGs live in a **`figures_described/`** folder. For many plots there is a **same-base-name `.txt`** file next to the PNG with **methods and interpretation** — read that first.

### Step 1 — BRANCH (`01_…/figures_described/`)

- **Atlas / tree views** — which regions show **strong global** Active−Passive separation or pass FDR-style thresholds (whole matrix context; **not** the same as single-phase Step 6).
- **PCA / UMAP / dendrogram** — **low-dimensional view of samples or regions**; useful for outliers and global structure.
- **Tables (CSV)** — per-region statistics, fold changes, optional bootstrap CIs (`trap_config` → `bootstrap_B`).

### Step 2 — Clustering sweep (`02_…/figures_described/`)

- **Silhouette vs K** — how well a given **K** separates clusters (higher is generally better; used to pick K).
- **Stability heatmap** — how often **pairs of regions** end up in the same cluster across resamples / K (robustness).
- **Sample PCA** — mice as points: do **Active/Passive** or **phase** separate in density space?

### Step 3 — v2 (`03_region_clustering_v2/figures_described/`)

- **Embedding / density heatmaps** — clusters of **regions** with similar TRAP patterns; **representative region** tracks per cluster.
- **Z-score / density panels** — often **within-phase z** per region across mice (matches `phase_AP_z_within_phase` when enabled).
- **`TRAP_downstream_input.mat`** — **no figure**; input to Step 4.

### Step 4 — Flip (`04_…/figures_described/`)

- **Flip summaries** — regions where **Reinstatement (A−P)** and **Withdrawal (A−P)** jointly satisfy a **condition** (magnitude/sign rules); includes **permutation**-style null summaries in CSVs.
- **Bar / ranking figures** — highlight regions with **consistent cross-phase** patterns (again: **not** the Step 6 ranking).

### Steps 6–7 — Phase AP & scenarios

- **Volcano** — x-axis: **mean Active − mean Passive** (or z-scale equivalent); y-axis: **−log10(p)**. Dots = regions; **highlighted** = pass your **raw p** or **FDR** rule (`trap_config`).
- **FDR / significance trees** — Allen hierarchy colored by **significance** or effect direction (whole-brain context).
- **Horizontal bar charts** — **top regions** by effect or **smallest p** among significant (or ranked) set.
- **Four-way / multi-panel plots** — **side-by-side** Active vs Passive structure across **phases** or conditions (exact layout depends on script).
- **Step 7 scenario folders** — same **spirit** as Step 6 but for **predefined directional stories** (see folder READMEs and scenario names).

### Step 8 — Within-group Rein vs With

- **Delta bar plots** — per region: **mean(Rein) − mean(With)** within **Active** or within **Passive**; **Wilcoxon p** across mice in that delivery group.
- **Volcano / ranked lists** — which regions **increase** or **decrease** from Withdrawal to Reinstatement **within** one delivery group.

### Step 9 — Forebrain subset

- **Same figure types as Steps 6–8**, restricted to **forebrain** (excludes cerebellum + brainstem); use when **global brainstem signal** should not dominate.

### Step 10 — Five-phase timeline (`10_five_phase_timeline/`)

- **Δ vs baseline heatmap** — rows = regions, columns = phases; values = **mean(phase) − mean(baseline)** within **Active** or **Passive** (raw density or z, per config).
- **Row-wise z across phases** — **shape** of change across the timeline (magnitude-invariant).
- **Line plots (top-N regions)** — **mean density (or z)** vs phase for the most **variable** regions.
- **Volcano per phase** — **Active vs Passive** **within that phase only** (same statistical idea as Step 6).

---

## Methods & statistics (EN + KR) — 추가 문서

**[`PIPELINE_DATA_STATISTICS_EN_KR.md`](PIPELINE_DATA_STATISTICS_EN_KR.md)** — Mermaid flow, **processing**, **group definitions**, **tests** (Steps 1–9).  
**[`RUN_MATLAB_AND_GITHUB_EN_KR.md`](RUN_MATLAB_AND_GITHUB_EN_KR.md)** — MATLAB run + GitHub push (EN + KR).

**Step 10 (five-phase):** [`STEP10_NEW_DATA_FIVE_PHASE_WORKFLOW.md`](STEP10_NEW_DATA_FIVE_PHASE_WORKFLOW.md)

**Run order (script index):** [`WORKFLOW.md`](WORKFLOW.md)  
**Where code vs outputs live:** [`FOLDERS_GUIDE.md`](FOLDERS_GUIDE.md)  
**Step 3 vs 6–9 sample alignment:** [`STEP_CONSISTENCY_3_vs_6_8.md`](STEP_CONSISTENCY_3_vs_6_8.md)  
**Technical index:** [`PIPELINE_ROADMAP.md`](PIPELINE_ROADMAP.md)  
**Warnings:** [`WARNINGS_EXPLAINED_EN_KR.md`](WARNINGS_EXPLAINED_EN_KR.md)  
**Output paths (legacy notes):** [`OUTPUT_GUIDE_EN_KR.md`](OUTPUT_GUIDE_EN_KR.md)

---

## Before you run (minimal checklist)

1. **`init_TRAP_pipeline`** once per MATLAB session.
2. **`trap_config.m`** — `csvPath`, outputs, `runMode`, optional `v2_sample_source`, `phase_AP_z_within_phase`, FDR/raw *p*.
3. **`TRAP_cohort_CSVs.txt`** — one density CSV path per line (cohort 1, 2, …); same `id` universe across files.
4. **`TRAP_sample_manifest.csv`** — `cohort_id`, `column_name`, `delivery`, `phase`, `include`.
5. **`MOUSE_COHORT.txt`** — optional human notes.

---

## Folder layout (code in repo)

| Folder | Role |
|--------|------|
| **`shared/`** | Loaders, FDR, export, phase normalization, AP helpers |
| **`Step_01_…` – `Step_05_…`** | BRANCH, correlations, clustering, flip, utilities |
| **`Step_06_phase_AP_contrasts/`** | Steps 6–10 drivers (`trap_run_*`) |

**Generated outputs** — under **`TRAP_OUTPUT/`** as in the step table above; many figures sit in **`figures_described/`** with **paired `.txt`** captions.

---

## Requirements

- MATLAB (Statistics Toolbox)
- Optional: UMAP (`run_umap`)
- Optional: Bioinformatics Toolbox — `BRANCH_analysis_TRAP_density.m` only (`mafdr`)

---

## License

Add a `LICENSE` file if you want open reuse.
