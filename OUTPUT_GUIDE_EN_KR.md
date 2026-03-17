# Where outputs go / 새 출력이 어디에 생기는지

## English

### What “v2 + downstream follow `trap_config`” means

- **Before:** Some scripts saved next to your CSV under `Downloads\…` or other hard-coded paths.
- **Now:** Paths come from **`trap_config.m`**. The pipeline root is the folder that contains `trap_config.m` (your `TRAP_pipeline` folder).
- **v2 outputs** go to:  
  **`TRAP_pipeline\TRAP_region_clusters_by_phase_density_v2\`**  
  including **`TRAP_downstream_input.mat`**.
- **Downstream / flip** read that `.mat` via **`C.downstream_mat`** — always the same folder unless you change **`C.v2_outDir`** in `trap_config.m`.

So everything stays **next to the code + CSV**, not under Downloads, unless you edit the config.

---

### Output folders from Step 1 → end (after `RUN_PIPELINE_ALL` or manual steps)

| Location | Contents |
|----------|----------|
| **`TRAP_OUTPUT\BRANCH_advanced\`** | **`BRANCH_stats_density.csv`** — per-region p, q, Cliff’s δ, fold change, mean Active−Passive; if `bootstrap_B>0`, also bootstrap CI columns. **`PairedTests_Withdrawal_density.csv`** — paired test Withdrawal Passive vs Active across regions. **PNG figures** (tree, PCA, UMAP, dendrogram) when MATLAB can write graphics (see note below). |
| **`TRAP_OUTPUT\clustering_sweep\`** | **New analysis:** silhouette vs K, stability heatmaps, sample PCA — **PNG** + console summary. |
| **`TRAP_region_clusters_by_phase_density_v2\`** | **v2:** embedding/density **PNGs**, **`RepRegions_*_Cluster*_fullnames.csv`**, **`TRAP_downstream_input.mat`**, **`Depth56_region_acronym_fullname.csv`**. |
| **`TRAP_OUTPUT\flip_downstream\`** | **New:** **`ConditionA/B/C_flip_full.csv`**, **`Flip_permutation_summary.csv`**, flip **PNGs** (top regions). |

**Note (figures):** If you run from **`matlab -batch`** (no display), some systems **do not save PNGs**. Run **`RUN_PIPELINE_ALL`** or each step in the **MATLAB desktop app** to generate all plots.

---

### New methods & plots (summary)

| Method | What it does |
|--------|----------------|
| **Manifest (`TRAP_sample_manifest.csv`)** | Delivery / phase / include per column; no need to encode rules in filenames. |
| **Bootstrap CI (Step 1)** | Optional 95% CI for mean(Active)−mean(Passive) per region (`bootstrap_B` in config). |
| **BY-FDR** | Optional stricter FDR (`C.fdrMethod = 'BY'`). |
| **Clustering sweep** | Tries K=2…Kmax, **mean silhouette** curve; **stability matrix** (how often two regions co-cluster); **PCA of mice** (samples as points). |
| **Flip + min \|Δ\|** | Only counts a region as “flip” if Active−Passive difference in Reinstatement **and** Withdrawal exceeds **`flip_min_abs_delta`**. |
| **Permutation null** | Shuffles Active/Passive labels within each phase; **p_perm** = how often random labelings give as many flip regions as observed. |

---

## 한국어 (Korean)

### “v2 + downstream이 trap_config를 따른다”는 뜻

- **예전:** 결과가 `Downloads` 옆 CSV 폴더 등 **경로가 코드에 박혀 있던** 곳에 저장됨.
- **지금:** **`trap_config.m`** 한 곳에서 경로 지정. 기본적으로 **`trap_config.m`이 있는 폴더(TRAP_pipeline)** 가 기준.
- **v2 결과:**  
  **`TRAP_pipeline\TRAP_region_clusters_by_phase_density_v2\`**  
  여기에 **`TRAP_downstream_input.mat`** 포함.
- **다운스트림(flip 등)** 은 **`C.downstream_mat`** 로 그 `.mat`을 읽음 → **v2와 같은 폴더** (원하면 `C.v2_outDir`만 바꾸면 됨).

즉, **코드/CSV와 같은 프로젝트 폴더 안**에 결과가 모이고, Downloads에 흩어지지 않습니다.

---

### 단계별 출력 폴더

| 경로 | 내용 |
|------|------|
| **`TRAP_OUTPUT\BRANCH_advanced\`** | **`BRANCH_stats_density.csv`**: 뇌영역별 Active vs Passive 검정, q값, Cliff’s δ, fold change, 평균 차이; bootstrap 켜면 **신뢰구간 열** 추가. **`PairedTests_Withdrawal_density.csv`**: Withdrawal에서 Paired 검정. **그림(PNG)** 은 MATLAB이 그래픽 저장 가능할 때 생성. |
| **`TRAP_OUTPUT\clustering_sweep\`** | **새 분석:** K에 따른 **실루엣 곡선**, **클러스터 안정성 히트맵**, **쥐(샘플) PCA 산점도**. |
| **`TRAP_region_clusters_by_phase_density_v2\`** | **v2:** 위상 임베딩/밀도 **PNG**, 대표뇌역 **CSV**, **`TRAP_downstream_input.mat`**, depth 5–6 **이름 CSV**. |
| **`TRAP_OUTPUT\flip_downstream\`** | **새:** 조건 A/B/C **flip 뇌역 목록 CSV**, **순열 검정 요약 CSV**, 상위 뇌역 **PNG**. |

**그림 주의:** 명령줄 **`matlab -batch`** 로 돌리면 환경에 따라 **PNG가 안 생길 수 있음**. **MATLAB GUI**에서 `RUN_PIPELINE_ALL` 또는 각 스크립트를 실행하면 그림까지 생성되는 경우가 많습니다.

---

### 새로 들어간 방법·그림 요약

| 방법 | 설명 |
|------|------|
| **매니페스트** | 샘플 열 이름마다 Active/Passive, phase, 포함 여부 지정 — 쥐 늘릴 때 **CSV만 수정**. |
| **Bootstrap CI** | (옵션) 영역별 Active−Passive **평균 차이의 95% CI**. |
| **BY-FDR** | (옵션) 다중비교 보수적 보정. |
| **Clustering sweep** | K 바꿔가며 **실루엣**, **부트스트랩식 공동클러스터 비율**, **샘플(쥐) PCA**. |
| **최소 \|Δ\|** | Reinstatement·Withdrawal 각각에서 Active−Passive 차이가 **임계값 이상**일 때만 flip으로 카운트. |
| **순열 검정** | phase 안에서 라벨을 섞어 **우연히 flip이 이만큼 나올 확률(p_perm)** 추정. |

---

### 한 번에 돌리기

```matlab
cd('...\TRAP_pipeline')
init_TRAP_pipeline
RUN_PIPELINE_ALL
```

빠른 테스트: `trap_config.m`에서 **`C.runMode = 'quick'`** (bootstrap 끔, 순열 500회, k-means 반복 감소).

---

*This guide matches a completed run on your machine (CSVs + `.mat` present). Regenerate PNGs in MATLAB GUI if needed.*
