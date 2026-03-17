# Adding mice / cohorts / 쥐·코호트 추가

## English

### One cohort (single CSV)

1. **`TRAP_cohort_CSVs.txt`** — one line: your CSV filename (or delete file to use `csvPath` in `trap_config.m` only).
2. **`TRAP_sample_manifest.csv`** — `cohort_id` = **1** for every row; set `column_name`, `delivery`, `phase`, `include`.

### New cohort = new CSV (2nd file, 3rd file, …)

1. Export the new cohort as the **same Allen atlas table** (same `id` list as cohort 1). Row order can differ; rows are matched by **`id`**.
2. **`TRAP_cohort_CSVs.txt`** — add a **new line** with the new CSV path (line 2 = cohort 2, line 3 = cohort 3, …).
3. **`TRAP_sample_manifest.csv`** — add rows with **`cohort_id` = 2** (or 3, …) and the **exact** density column name from **that** CSV.
4. Run (same as always):

```matlab
cd('...\TRAP_pipeline')
init_TRAP_pipeline
RUN_PIPELINE_ALL
```

All steps (BRANCH, clustering, v2, flip) use **pooled** mice from every listed cohort.

---

## 한국어

### 코호트 하나 (CSV 하나)

- **`TRAP_cohort_CSVs.txt`**에 CSV 한 줄.
- **`TRAP_sample_manifest.csv`**에서 **`cohort_id` = 1**, 열 이름·그룹 정보 입력.

### CSV를 추가할 때마다 새 코호트

1. 새 CSV는 **코호트1과 동일한 Allen `id` 목록**이어야 함 (행 순서는 달라도 됨, **id**로 맞춤).
2. **`TRAP_cohort_CSVs.txt`**에 **한 줄 추가** (2번째 줄 = 코호트 2).
3. **`TRAP_sample_manifest.csv`**에 **`cohort_id` = 2** 와 해당 CSV의 **density 열 이름**으로 행 추가.
4. MATLAB에서 **`init_TRAP_pipeline`** → **`RUN_PIPELINE_ALL`** 그대로 실행.

분석은 나열된 모든 CSV에서 **쥐를 합쳐서** 진행됩니다.
