# TRAP pipeline — Python (Steps 1–7)

Mirrors **`trap_pipeline_matlab/`** logic. **Step 6 is the core question:**

> For **each** brain region, is TRAP density **significantly different between Active and Passive mice**, tested **separately** in Reinstatement and in Withdrawal?

- **Default test:** Wilcoxon rank-sum (`ranksum`) — robust for small *n* and skewed density.
- **Alternative:** `--test welch` (unequal-variance *t*-test).

## Layout

```
trap_pipeline_py/
  run_pipeline.py          # entry point
  trap_pipeline/
    config.py
    io_loader.py
    fdr.py
    stats_ap.py
    plots.py
    steps/
      step01_branch.py
      step02_clustering.py
      step03_region_clusters.py
      step04_flip.py
      step05_export.py
      step06_regionwise.py   # CORE
      step07_followup.py
```

## Run

```bash
cd trap_pipeline_py
pip install -r requirements.txt
python run_pipeline.py --csv ../path/to/density.csv --manifest ../path/to/TRAP_sample_manifest.csv --out TRAP_PYTHON_OUTPUT
python run_pipeline.py --steps 6 --csv ... --manifest ...   # only Step 6
```

Outputs: `06_regionwise_Active_vs_Passive/` (tables + figures), `07_followup/` (phase delta + exploratory CSVs).

The legacy script `run_trap_pipeline.py` is superseded by `run_pipeline.py`.
