# Step 00 vs Steps 1–11 — which CSV columns are “mice”?

Exports often repeat **one mouse** across several columns:

| Kind of column | Example header fragment | Use as a mouse sample? |
|----------------|-------------------------|-------------------------|
| Cell **count** | `… count` | No |
| **Density** | `… density (cells/mm^3)` | **Yes** (this is what TRAP analyzes) |
| **Volume** | `… volume (mm^3)` | No |
| **Aggregate** | `AVERAGE …` | No |

## Step 00 (`trap_run_mouse_qc_density`)

If the manifest lists a sample with **`include=0`**, that column is **dropped** in Step 00 (same intent as Steps 1+).

With default **`C.mouse_qc_use_all_csv_columns = true`**, Step 00 **does not** treat every numeric column as a mouse.

- It keeps only columns whose name **contains** **`C.mouse_qc_density_column_header_substring`** (default: `density (cells/mm^3)`), **case-insensitive**.
- It **drops** columns whose name contains **`average`** (so pooled AVERAGE columns are excluded).

So count / volume / AVERAGE columns are **not** duplicated as extra “mice.”

Optional: set **`mouse_qc_use_all_csv_columns = false`** to load **only** manifest rows with **`include=1`** (same column list as Steps 1+).

## Steps 1–11 (`trap_load_pooled_density_LR`)

These steps **never** auto-scan all numeric columns.

- They load **only** columns named in **`TRAP_sample_manifest.csv`** (`include=1`), with **`column_name`** matching the file header **exactly**.
- The “wrong” behavior (treating count + density + volume each as a separate mouse) **cannot** happen unless **you** add multiple manifest rows pointing at those different column types.

**Recommendation:** Each manifest row’s **`column_name`** should be the **density** column for that mouse (e.g. `HaLi_… density (cells/mm^3)`), not the count or volume column.

## Cohort files can have different numbers of atlas rows

Exports may differ slightly in row count (e.g. CSV 1681 rows vs Excel 1677). Step 00 accepts density columns if they match **that file’s** row count; values are aligned to cohort 1’s atlas via `id` when building the pooled matrix (same as `trap_load_pooled_density_LR`).

## Summary

| | Step 00 (default) | Steps 1–11 |
|--|-------------------|------------|
| How columns are chosen | Auto: density substring + no “average” | Manifest: explicit `column_name` per sample |
| Can count/volume be mistaken for mice? | No (filtered out) | Only if listed separately in manifest |

See also: **`trap_config.m`** (`mouse_qc_density_column_header_substring`), **`README.md`**, **`RUN_MATLAB_AND_GITHUB_EN_KR.md`**.
