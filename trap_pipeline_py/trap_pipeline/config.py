"""Paths and analysis options (align with trap_pipeline_matlab/trap_config.m)."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class TrapConfig:
    csv_path: Path
    manifest_path: Path
    out_root: Path
    # Step 6: per-region Active vs Passive — ranksum (Wilcoxon) default; welch optional
    ap_test: str = "ranksum"  # "ranksum" | "welch"
    use_fdr: bool = False
    p_raw: float = 0.05
    fdr_alpha: float = 0.05
    flip_n_perm: int = 2000
    flip_min_abs_delta: float = 0.5

    @property
    def branch_dir(self) -> Path:
        return self.out_root / "01_BRANCH_tables_and_figures"

    @property
    def cluster_dir(self) -> Path:
        return self.out_root / "02_clustering_sweep"

    @property
    def v2_dir(self) -> Path:
        return self.out_root / "03_region_clustering_v2"

    @property
    def flip_dir(self) -> Path:
        return self.out_root / "04_flip_downstream"

    @property
    def step06_dir(self) -> Path:
        return self.out_root / "06_regionwise_Active_vs_Passive"

    @property
    def step07_dir(self) -> Path:
        return self.out_root / "07_followup"
