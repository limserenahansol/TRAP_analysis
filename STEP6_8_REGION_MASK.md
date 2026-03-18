# Steps 6–8 use the same brain regions as Step 3

With **`phase_AP_region_mask_step3 = true`** in `trap_config.m` (default):

- Region list = **exactly the Step 3 v2 depth/parent rule** (`C.v2_depth_rule`):
  - **`hierarchy567`** — depth-7 (non-“layer”), or depth-6 unless replaced by d7 under same d5 parent, or depth-5 leaves.
  - **`depth56_fixed`** — depth-6 + depth-5 nodes without a depth-6 child.

Weird / deep atlas rows outside that set are **dropped** in Steps 6, 7, and 8 so statistics and plots match the regions you already use downstream.

To use the full ~800+ L/R-pooled rows again: set **`C.phase_AP_region_mask_step3 = false`**.

## Exclude-phase samples

With **`phase_AP_drop_exclude_samples = true`** (default), rows with **phase = Exclude** are removed before Steps 6–8 — matching Step 3 when samples come from the manifest. Step 1 BRANCH keeps those columns.

## Step 3 `all_csv` vs manifest

If **`v2_sample_source = 'all_csv'`**, Step 3/4 can include samples not in the manifest; Steps 6–8 remain manifest-based. See **`PIPELINE_CONSISTENCY_STEPS1_through_8.md`**.
