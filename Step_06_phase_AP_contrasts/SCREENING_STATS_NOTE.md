# Screening many brain regions: Active vs Passive (short reference note)

## The problem
You test **hundreds of regions** × **one contrast** (Active vs Passive). Uncorrected **p < 0.05** gives **~5% false positives** among null regions → many spurious hits if you only use raw p.

## Common approaches in practice

1. **FDR (Benjamini–Hochberg)**  
   Control expected fraction of false discoveries among called “significant” regions. Standard in genomics/imaging-style mass univariate screens. With **very small n per group**, power is low → **few** regions pass strict **q ≤ 0.05**; some papers use **q = 0.10–0.20** for **exploratory** screens.

2. **Effect size + FDR (recommended hybrid)**  
   Rank regions by **|mean(Active) − mean(Passive)|**, **Cohen’s d**, or **Cliff’s delta**, and among them apply **FDR** (or the reverse: FDR first, then report largest effects). This matches reports that care about **biological magnitude**, not only p-values.

3. **Welch t-test vs Wilcoxon**  
   Welch: parametric, good if distributions are roughly normal. **Wilcoxon / ranksum**: robust, often preferred for **small n** and skewed density data. Both are fine; the critical choice is **multiple-comparison control** (or explicit exploratory labeling).

4. **Two-stage / validation**  
   Screen on one cohort or relaxed FDR → **replicate** in independent mice or with focused ROI analysis (common in behavioral neuroscience).

## What this pipeline does
- **`phase_AP_use_fdr = false`**: **Exploratory** screen (raw p) — good for **hypothesis generation**, like many bar-plot figures; **not** strong family-wise control.  
- **`phase_AP_use_fdr = true`**: **Conservative** screen across all atlas regions.  
- **`trap_run_phase_delta_screening`**: **Additional** criterion — large **|Δ_Rein − Δ_With|** flags regions where the **Active–Passive gap changes strongly** between Reinstatement and Withdrawal (interaction-style importance), **complementary** to simple Active vs Passive tests.

Always report **n mice per group per phase** next to any screen.
