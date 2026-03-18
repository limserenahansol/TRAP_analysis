"""Step 2: sample PCA + silhouette sweep (depth 5-6 regions)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from trap_pipeline.config import TrapConfig


def run(cfg: TrapConfig, dens: np.ndarray, node: pd.DataFrame, delivery, phase) -> None:
    out = cfg.cluster_dir / "figures_described"
    tab = cfg.cluster_dir / "tables"
    out.mkdir(parents=True, exist_ok=True)
    tab.mkdir(parents=True, exist_ok=True)

    m = (node["depth"] >= 5) & (node["depth"] <= 6)
    X = dens[m.values, :].T
    X = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2)
    Z = pca.fit_transform(X)
    fig, ax = plt.subplots(figsize=(7, 6))
    for ph in np.unique(phase):
        sel = phase == ph
        ax.scatter(Z[sel, 0], Z[sel, 1], label=str(ph), alpha=0.75, s=40)
    ax.set_title("Step2 PCA samples (depth 5-6 regions)")
    ax.legend()
    plt.tight_layout()
    fig.savefig(out / "01_PCA_samples.png", dpi=160)
    plt.close()

    rows = []
    n_s = X.shape[0]
    for k in range(2, min(9, n_s)):
        km = KMeans(n_clusters=k, random_state=0, n_init=20).fit(X)
        if len(np.unique(km.labels_)) < 2:
            continue
        sil = silhouette_score(X, km.labels_)
        rows.append({"K": k, "silhouette": sil})
    pd.DataFrame(rows).to_csv(tab / "silhouette_sweep.csv", index=False)
