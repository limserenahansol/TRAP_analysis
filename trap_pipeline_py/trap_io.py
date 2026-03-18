"""Load TRAP density CSV + manifest; L/R average (same logic as MATLAB)."""
from __future__ import annotations

import re
import numpy as np
import pandas as pd


def _read_manifest(path: str) -> pd.DataFrame:
    m = pd.read_csv(path, encoding="utf-8-sig")
    cols = {c.lower().strip("\ufeff"): c for c in m.columns}
    def pick(*names):
        for n in names:
            for k, v in cols.items():
                if k.replace(" ", "_") == n.replace(" ", "_") or k == n:
                    return m[v]
        raise KeyError(names[0])
    out = pd.DataFrame({
        "column_name": pick("column_name", "ColumnName").astype(str).str.strip(),
        "delivery": pick("delivery", "Delivery").astype(str).str.strip(),
        "phase": pick("phase", "Phase").astype(str).str.strip(),
    })
    inc_col = None
    for c in m.columns:
        if str(c).lower().strip() == "include":
            inc_col = m[c]
            break
    if inc_col is not None:
        inc = pd.to_numeric(inc_col, errors="coerce").fillna(0)
        out["include"] = (inc == 1) | (inc_col.astype(str).str.strip() == "1")
    else:
        out["include"] = True
    out = out[out["include"]].reset_index(drop=True)
    return out


def load_lr_averaged_density(csv_path: str, manifest_path: str):
    """
    Returns:
        dens: (n_regions, n_samples) float
        node: DataFrame id, acronym, depth, parent_structure_id, name
        sample_names: list[str]
        delivery: np.ndarray str
        phase: np.ndarray str
    """
    df = pd.read_csv(csv_path, encoding="utf-8-sig")
    meta = ["id", "name", "acronym", "parent_structure_id", "depth"]
    for c in meta:
        if c not in df.columns:
            raise ValueError(f"CSV missing column {c}")

    man = _read_manifest(manifest_path)
    col_names = []
    for _, row in man.iterrows():
        c = row["column_name"]
        matches = [x for x in df.columns if x == c or x.strip('"') == c.strip('"')]
        if not matches:
            raise FileNotFoundError(f"Manifest column not in CSV: {c!r}")
        col_names.append(matches[0])

    acr = df["acronym"].astype(str)
    is_l = acr.str.endswith("-L")
    is_r = acr.str.endswith("-R")
    is_g = ~(is_l | is_r)
    keep = is_l | is_g
    node_full = df.loc[keep, meta].reset_index(drop=True)
    acr_full = acr[keep].values

    D = df.loc[keep, col_names].to_numpy(dtype=float)
    n_r, n_s = D.shape
    dens = np.full((n_r, n_s), np.nan)
    acr_list = acr_full.tolist()
    acr_to_idx = {a: i for i, a in enumerate(acr_list)}

    for ii in range(n_r):
        ac = acr_list[ii]
        if ac.endswith("-L"):
            base = ac[:-2]
            r_ac = base + "-R"
            j = acr_to_idx.get(r_ac)
            if j is not None:
                dens[ii, :] = (D[ii, :] + D[j, :]) / 2.0
            else:
                dens[ii, :] = D[ii, :]
        else:
            dens[ii, :] = D[ii, :]

    node = node_full.copy()
    node["acronym"] = node["acronym"].astype(str).str.replace("-L$", "", regex=True)

    delivery = man["delivery"].str.strip().str.capitalize()
    delivery = delivery.replace({"Active": "Active", "Passive": "Passive"})
    phase = man["phase"].str.strip()

    return dens, node, col_names, delivery.values, phase.values
