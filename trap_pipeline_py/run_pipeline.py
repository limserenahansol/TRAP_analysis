#!/usr/bin/env python3
"""
TRAP pipeline Steps 1-7 (Python), aligned with trap_pipeline_matlab.

  cd behavior_task/trap_pipeline_py
  pip install -r requirements.txt
  python run_pipeline.py --csv PATH --manifest PATH [--out DIR] [--test ranksum|welch]

Step 6 (CORE): each brain region, Active vs Passive in Reinstatement and Withdrawal separately.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from trap_pipeline.config import TrapConfig
from trap_pipeline.io_loader import load_lr_averaged_density
from trap_pipeline.steps import step01_branch, step02_clustering, step03_region_clusters
from trap_pipeline.steps import step04_flip, step05_export, step06_regionwise, step07_followup


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--out", default="TRAP_PYTHON_OUTPUT")
    ap.add_argument("--test", choices=("ranksum", "welch"), default="ranksum")
    ap.add_argument("--fdr", action="store_true", help="Use FDR q<=0.05 instead of raw p<=0.05")
    ap.add_argument("--steps", default="1-7", help="e.g. 1-7 or 6 only")
    args = ap.parse_args()

    cfg = TrapConfig(
        csv_path=Path(args.csv),
        manifest_path=Path(args.manifest),
        out_root=Path(args.out).resolve(),
        ap_test=args.test,
        use_fdr=args.fdr,
    )
    dens, node, _, delivery, phase = load_lr_averaged_density(
        str(cfg.csv_path), str(cfg.manifest_path)
    )

    def parse_range(s: str) -> set[int]:
        if "-" in s:
            a, b = s.split("-", 1)
            return set(range(int(a), int(b) + 1))
        return {int(x.strip()) for x in s.split(",")}

    steps = parse_range(args.steps)
    runners = {
        1: lambda: step01_branch.run(cfg, dens, node, delivery, phase),
        2: lambda: step02_clustering.run(cfg, dens, node, delivery, phase),
        3: lambda: step03_region_clusters.run(cfg, dens, node),
        4: lambda: step04_flip.run(cfg, dens, node, delivery, phase),
        5: lambda: step05_export.run(cfg, node),
        6: lambda: step06_regionwise.run(cfg, dens, node, delivery, phase),
        7: lambda: step07_followup.run(cfg, dens, node, delivery, phase),
    }
    for k in sorted(steps):
        if k in runners:
            print(f"--- Step {k} ---")
            runners[k]()
    print(f"Done. Output: {cfg.out_root}")


if __name__ == "__main__":
    main()
