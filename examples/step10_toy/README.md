# Step 10 toy / smoketest (no real TRAP data)

Generates **synthetic** density CSV + manifest (20 mice × 5 phases × Active/Passive), runs **`trap_run_phase5_timeline_analysis`** with outputs under **`TRAP_OUTPUT/step10_toy_smoketest/`** so your real **`10_five_phase_timeline/`** folder is untouched.

## Run

```matlab
cd('C:\path\to\TRAP_analysis')
addpath(fullfile(pwd, 'examples', 'step10_toy'))
init_TRAP_pipeline
run_step10_toy_smoketest
```

Optional: also run the forebrain duplicate suite (may drop many toy regions):

```matlab
run_step10_toy_smoketest(struct('phase5_run_forebrain_duplicate', true))
```

Generated files (safe to delete) live in **`examples/step10_toy/generated/`** (gitignored).

## Fix applied for real data too

If **Baseline** has no mice, within-group Δ vs baseline is undefined and the old code could error in the delta tree plot. The pipeline now handles empty / all-NaN highlights without the `&&` scalar error.

When you have real five-phase data, add **Baseline** rows to **`TRAP_sample_manifest.csv`** so Q1/Q2 within-group baselines are meaningful.
