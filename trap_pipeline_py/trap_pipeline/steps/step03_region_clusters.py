"""Step 3: k-means on region profiles (samples as features), depth 5–6."""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

from trap_pipeline.config import TrapConfig


def run(cfg: TrapConfig, dens: np.ndarray, node: pd.DataFrame) -> None:
    out = cfg.v2_dir
    tab = out / "tables"
    tab.mkdir(parents=True, exist_ok=True)
    m = (node["depth"] >= 5) & (node["depth"] <= 6)
    X = dens[m.values, :]
    K = 6
    labels = KMeans(n_clusters=K, random_state=0, n_init=30).fit_predict(X)
    sub = node.loc[m].copy()
    sub["cluster"] = labels
    sub.to_csv(tab / "region_clusters_kmeans.csv", index=False)
