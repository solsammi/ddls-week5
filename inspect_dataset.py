"""Inspect the supplied PBMC dataset and summarise clusters 6 and 7."""
from pathlib import Path
import zipfile
import scanpy as sc
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
archive = ROOT / "ddls-week5-s1-junk-or-signal-dataset.zip"
data_path = ROOT / "data" / "pbmc3k.h5ad"
data_path.parent.mkdir(exist_ok=True)
if not data_path.exists():
    with zipfile.ZipFile(archive) as z:
        z.extract("data/pbmc3k.h5ad", ROOT)

ad = sc.read_h5ad(data_path)
obs = ad.obs.copy()
labels = obs["leiden"].astype(str)

print("DATASET")
print(f"path={data_path}")
print(f"shape={ad.n_obs} cells x {ad.n_vars} genes")
print(f"X_type={type(ad.X).__name__}")
print(f"layers={list(ad.layers.keys())}")
print(f"obs_columns={list(obs.columns)}")
print(f"obsm_keys={list(ad.obsm.keys())}")
print(f"cluster_labels={sorted(labels.unique(), key=lambda x: int(x))}")
print(f"missing_obs_total={int(obs.isna().sum().sum())}")
print(f"duplicate_cell_ids={int(obs.index.duplicated().sum())}")
print(f"duplicate_expression_rows={int(pd.DataFrame(ad.X.toarray() if hasattr(ad.X, 'toarray') else np.asarray(ad.X)).duplicated().sum())}")

print("\nCLUSTER_SUMMARY")
summary = []
for cluster in sorted(labels.unique(), key=lambda x: int(x)):
    group = obs.loc[labels == cluster]
    row = {"cluster": cluster, "n_cells": len(group)}
    for field in ("n_genes", "total_counts", "pct_mito"):
        values = pd.to_numeric(group[field], errors="coerce")
        row[f"{field}_min"] = values.min()
        row[f"{field}_median"] = values.median()
        row[f"{field}_mean"] = values.mean()
        row[f"{field}_max"] = values.max()
    summary.append(row)
print(pd.DataFrame(summary).to_string(index=False))

print("\nTARGET_CLUSTERS")
for cluster in ("6", "7"):
    group = obs.loc[labels == cluster]
    print(f"cluster={cluster} n_cells={len(group)}")
    print(group[["leiden", "n_genes", "total_counts", "pct_mito"]].to_string())
