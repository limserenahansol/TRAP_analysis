# Major-division labels on region plots

When **`phase_AP_plot_major_class = true`** (default), tick labels look like **`ACA (cortex, prefrontal / anterior cingulate)`**, **`CP (striatum)`**, **`B (brainstem)`**.

| Examples | Region group |
|----------|----------------|
| **thalamus** / **hypothalamus** | Allen TH- / HY- |
| **midbrain** | MB- |
| **brainstem** | Pons, medulla, etc. (not TH/HY) |
| **cerebellum** | CB |
| **hippocampus** / **olfactory** | HPF, OLF |
| **amygdala** | BLA, CEA, sAMY, … |
| **striatum** / **basal ganglia (pallidum)** | STR, CP, ACB / PAL |
| **cortex, motor** / **somatosensory** / **visual** / **auditory** / **retrosplenial** | MO, SS, VIS, AUD, RSP |
| **cortex, prefrontal / anterior cingulate** | ORB, FRP, PL, ILA, ACA |
| **cortex, posterior parietal** / **temporal / association** | PTLp, TEa, ENT, … |
| **cortex, gustatory** / **insular** | GU, VISC |

Turn off: **`C.phase_AP_plot_major_class = false`**.

## Steps 1–4 (same toggle)

- **Step 1** (`run_BRANCH_TRAP_density2.m` tree labels on significant nodes).
- **Step 3** v2 / Top50: representative-region plots + embedding text; Top50 rank plots.
- **Step 4** (`trap_run_flip_advanced`, `TRAP_run_downstream`): flip Top‑N x-axis labels.

CSV columns still use **ParentD4** where applicable; only **plot** tick/text labels use the major class.

## Step 9

**`trap_run_step9_forebrain_exclude_bs_cb`** (also at end of `RUN_PIPELINE_ALL`) writes under  
`TRAP_OUTPUT/09_forebrain_no_brainstem_cerebellum/` the same figures/CSVs as Steps 6–8 but **only cerebrum + thalamus + hypothalamus** (cerebellum and rest of brainstem removed).
