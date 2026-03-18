# Phase labels in `TRAP_sample_manifest.csv`

| You type in manifest | Pipeline uses |
|----------------------|---------------|
| **Reexposure**, Re-exposure, reexposure, … | **Reinstatement** |
| Reinstatement | Reinstatement |
| Withdrawal | Withdrawal |

So **Reexposure and Reinstatement columns are pooled** as one phase for Steps 6–8 (and any code using `GroupPhase`).

### Your design (example)

- **Reinstatement / Reexposure:** Active = **3** mice, Passive = **2** mice → **5** density columns in that phase.
- **Withdrawal:** Active = **1**, Passive = **1** → **2** columns.

Total **7** manifest rows (`include=1`), each pointing to a column in the cohort CSV.

After editing the manifest, re-run the pipeline. Figures will show  
`This phase (Reinstatement): 3 Active, 2 Passive` (Reexposure rows count toward Reinstatement).
