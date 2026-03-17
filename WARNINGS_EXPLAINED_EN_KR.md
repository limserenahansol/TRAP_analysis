# MATLAB warnings you may see / 경고 설명

## English

| Message | Meaning | What to do |
|---------|---------|------------|
| **UMAP not found — using PCA** | The optional UMAP toolbox (`run_umap`) is not on the path. | Install UMAP for MATLAB if you want nonlinear embedding; PCA is still valid. With **only 2 mice in Withdrawal**, PCA has at most 1–2 meaningful directions — **adding mice helps a lot**. |
| **Columns linearly dependent … PCA** | Few samples (e.g. 2) → design matrix is rank-deficient. MATLAB still draws PC1/PC2 but PC2 may be weak or redundant. | Expected with small n per phase. Interpret Withdrawal embedding cautiously; Reinstatement with more mice is more stable. |
| **k-means Failed to converge in 100 iterations** | Some random k-means starts did not finish in 100 steps (high-dimensional region space). | **Mitigation:** code now uses **MaxIter 500**. Warnings can still appear occasionally; results use the best replicate. Increasing `v2_kmeans_replicates` in `trap_config` can help. |
| **Exported image displays axes toolbar** | `exportgraphics` captured the small axes UI. | Figures now turn **toolbar off** before save. Re-run to get clean PNGs. |

---

## 한국어

| 메시지 | 의미 | 대응 |
|--------|------|------|
| **UMAP not found** | UMAP 도구가 없어 **PCA**로 대체. | UMAP 설치 가능; **Withdrawal에 쥐가 2마리뿐이면** PCA도 정보가 제한적 → **쥐 수 늘리기**가 가장 좋음. |
| **PCA linearly dependent** | 샘플 수가 매우 적어 **랭크 부족**. | Withdrawal 위상도는 해석 조심. |
| **k-means 수렴 실패** | 일부 반복이 100회 안에 안 끝남. | **MaxIter 500**으로 완화됨. 가끔 경고는 남을 수 있음. |
| **axes toolbar** | 그림에 툴바가 찍힘. | 저장 전 툴바 끔 처리함. 다시 실행하면 깨끗한 PNG. |
