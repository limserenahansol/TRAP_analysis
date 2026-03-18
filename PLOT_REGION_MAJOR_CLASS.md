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

## Step 9

**`trap_run_step9_forebrain_exclude_bs_cb`** (also at end of `RUN_PIPELINE_ALL`) writes under  
`TRAP_OUTPUT/09_forebrain_no_brainstem_cerebellum/` the same figures/CSVs as Steps 6–8 but **only cerebrum + thalamus + hypothalamus** (cerebellum and rest of brainstem removed).
