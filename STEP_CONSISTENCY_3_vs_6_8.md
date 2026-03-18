# Step 3 vs Steps 6–9 — 통일 규칙

**Step 4 vs Step 6:** 질문이 다릅니다 (한 phase A vs P vs 두 phase 조합). **상위 영역이 안 겹쳐도 정상**입니다. 아래는 **데이터 파이프라인**만 맞출 때의 옵션입니다.

## 1) 스케일 (Z-score, phase 내)

- **Step 3** `03_rep_regions_ZSCORED_within_phase_*.png`: 각 **phase** 안에서, **영역(행)마다** 그 phase의 모든 마우스에 대해 `zscore(..., 0, 2)` → 해당 phase에서 평균 0, 표준편차 1.
- **Steps 6–9**는 동일하게 **`trap_zscore_within_phase_columns`** 적용 후 분석·플롯 (설정: `C.phase_AP_z_within_phase = true`, 기본값).

`phase_AP_z_within_phase = false` 이면 예전처럼 raw density (cells/mm³).

## 2) 샘플(열) 통일

- Steps **6–8–9**는 항상 **매니페스트**(`TRAP_sample_manifest.csv`, include=1) 열만 사용 (`trap_load_pooled_density_LR`).
- **Step 3 v2**는 `trap_config.m`의 **`v2_sample_source`**에 따라 달라짐:
  - **`'manifest'`** → Step 6–8과 **동일한 마우스 집합**.
  - **`'all_csv'`** → CSV의 추가 열 포함 → Step 6–8과 **샘플이 다를 수 있음**.

**완전히 맞추려면:** `C.v2_sample_source = 'manifest';`

## 3) 뇌영역 라벨

- Step 6–8 바 플롯: `ACA (cortex, …)`, `PVT (thalamus)` 등 **`trap_region_plot_tick_labels`** (Allen 부모 + 일반명).
- Step 3 대표영역 플롯은 depth-4 부모 약어 `(MBmot)` 형태 — 원하면 Step 3도 동일 스타일로 바꿀 수 있음(별도 작업).

## 4) Step 4 vs Step 6 — 질문은 의도적으로 다름

| | Step 6 | Step 4 Flip |
|--|--------|-------------|
| 질문 | 한 phase에서 A vs P | Rein·With **조합** 조건 (A/B/C) |
| Top 영역 | phase별 유의/분리 | 플립 규칙 만족 — **6과 겹칠 필요 없음** |

### (선택) 같은 **데이터**로 3·6·8을 맞출 때만

| 항목 | 설정 |
|------|------|
| 샘플 | `v2_sample_source = 'manifest'` |
| 영역 | 동일 `v2_depth_rule` + Step 6 region mask |
| 스케일 | `phase_AP_z_within_phase = true` (기본) |

**같은 Top 목록**이 필요하면 별도로 CSV에서 id 필터·교집합 후처리.
