# Steps 1–8: rule alignment

| Rule | Step 1 (BRANCH) | Step 2 | Step 3 / 4 (v2) | Steps 6–8 |
|------|-----------------|--------|-----------------|-----------|
| **Data load** | `trap_load_pooled_density_LR` | — | manifest or `all_csv` | `trap_load_pooled_density_LR` (manifest + L+R pooled) ✓ |
| **Reexposure vs Reinstatement** | normalized in loader | — | — | same (`trap_normalize_manifest_phase`) ✓ |
| **Brain regions** | full L/R-pooled atlas | — | `v2_depth_rule` subset | **same subset** if `phase_AP_region_mask_step3=true` ✓ |
| **Exclude-phase samples** | **kept** in matrix | — | **dropped** in manifest mode | **dropped** if `phase_AP_drop_exclude_samples=true` ✓ |
| **v2_sample_source = all_csv** | — | — | extra CSV columns possible | Steps 6–8 still **manifest-only** — **not** sample-identical to v2 unless v2 uses manifest |

**Step 1 Active vs Passive** mixes all phases in one ranksum; **Step 6** is Active vs Passive **within Reinstatement** (or Withdrawal) separately — same test type on mice, different hypothesis.

**Config:** `trap_config.m` — `phase_AP_region_mask_step3`, `phase_AP_drop_exclude_samples`, `v2_depth_rule`, `v2_sample_source`.
