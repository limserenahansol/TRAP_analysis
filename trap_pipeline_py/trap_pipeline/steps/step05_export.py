"""Step 5: export region list depth ≥ 5."""
from __future__ import annotations

from trap_pipeline.config import TrapConfig


def run(cfg: TrapConfig, node) -> None:
    p = cfg.out_root / "05_export_regions"
    p.mkdir(parents=True, exist_ok=True)
    sub = node[node["depth"] >= 5][["id", "acronym", "name", "depth"]]
    sub.to_csv(p / "regions_depth5plus.csv", index=False)
