"""Preliminary clustering stability analysis across 3 seeds, neighbours, resolutions."""
from pathlib import Path
from itertools import product
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
base = sc.read_h5ad(ROOT / "data/pbmc3k.h5ad")
original = base.obs["leiden"].astype(str).to_numpy()

# Preprocess once; vary graph construction and Leiden randomness. No repeated UMAP.
ad = base.copy()
if "counts" in ad.layers:
    ad.X = ad.layers["counts"].copy()
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat", subset=True)
sc.pp.scale(ad, max_value=10)
sc.tl.pca(ad, n_comps=40, svd_solver="arpack", random_state=0)

seeds = (0, 1, 2)
neighbours = (10, 15, 30)
resolutions = (0.5, 0.8, 1.0)
assignments = {}
run_rows = []
for seed, n_neighbors, resolution in product(seeds, neighbours, resolutions):
    graph = f"neighbors_{seed}_{n_neighbors}"
    sc.pp.neighbors(ad, n_neighbors=n_neighbors, n_pcs=30, random_state=seed, key_added=graph)
    key = f"cluster_{seed}_{n_neighbors}_{str(resolution).replace('.', '_')}"
    sc.tl.leiden(ad, resolution=resolution, neighbors_key=graph, key_added=key, random_state=seed, flavor="igraph", n_iterations=2)
    labels = ad.obs[key].astype(str).to_numpy()
    assignments[key] = labels
    run_rows.append({"run": key, "seed": seed, "n_neighbors": n_neighbors, "resolution": resolution, "n_clusters": len(set(labels)), "ARI_vs_original": adjusted_rand_score(original, labels), "NMI_vs_original": normalized_mutual_info_score(original, labels)})

run_df = pd.DataFrame(run_rows)
run_df.to_csv(OUT / "stability_runs.csv", index=False)
assignment_df = pd.DataFrame(assignments, index=ad.obs_names)
assignment_df.insert(0, "original_cluster", original)
assignment_df.to_csv(OUT / "stability_cell_assignments.csv")

# Label-invariant cell consistency: pairwise co-clustering frequency across all runs.
label_matrix = assignment_df.drop(columns="original_cluster").to_numpy()
co_cluster = np.zeros((ad.n_obs, ad.n_obs), dtype=np.float32)
for labels in label_matrix.T:
    co_cluster += labels[:, None] == labels[None, :]
co_cluster /= label_matrix.shape[1]
np.save(OUT / "stability_co_cluster_frequency.npy", co_cluster)

# For each cell, frequency of co-clustering with its original-cluster peers.
cell_rows = []
for i, cell in enumerate(ad.obs_names):
    peers = np.flatnonzero(original == original[i])
    peers = peers[peers != i]
    same_original = float(co_cluster[i, peers].mean()) if len(peers) else np.nan
    other = np.flatnonzero(original != original[i])
    outside = float(co_cluster[i, other].mean()) if len(other) else np.nan
    cell_rows.append({"cell_id": str(cell), "original_cluster": original[i], "same_original_cocluster_frequency": same_original, "other_cluster_cocluster_frequency": outside})
cell_df = pd.DataFrame(cell_rows)
cell_df.to_csv(OUT / "stability_cell_scores.csv", index=False)

cluster_rows = []
for cluster in sorted(set(original), key=int):
    x = cell_df[cell_df.original_cluster == cluster]
    cluster_rows.append({"cluster": cluster, "n_cells": len(x), "median_same_original_cocluster_frequency": x.same_original_cocluster_frequency.median(), "min_same_original_cocluster_frequency": x.same_original_cocluster_frequency.min(), "max_same_original_cocluster_frequency": x.same_original_cocluster_frequency.max(), "median_outside_cocluster_frequency": x.other_cluster_cocluster_frequency.median()})
cluster_df = pd.DataFrame(cluster_rows)
cluster_df.to_csv(OUT / "stability_cluster_scores.csv", index=False)

report = ["# Preliminary clustering stability analysis", "", "## Design", "- Three random seeds: 0, 1, 2.", "- Three neighbour counts: 10, 15, 30.", "- Three Leiden resolutions: 0.5, 0.8, 1.0.", "- Total: 27 clustering runs.", "- Preprocessing and PCA were performed once; neighbour graphs and Leiden assignments were varied. UMAP was not recomputed because this analysis evaluates cluster assignment stability.", "", "## Run summary", run_df.to_string(index=False), "", "## Cluster stability", cluster_df.to_string(index=False), "", "## Interpretation", "- Co-clustering frequency is label-invariant: it measures how often cells are assigned to the same cluster across runs, regardless of numeric cluster names.", "- `same_original_cocluster_frequency` asks whether a cell repeatedly stays with cells from its supplied cluster; `other_cluster_cocluster_frequency` identifies possible mixing with other supplied clusters.", "- This is preliminary evidence, not proof of biological identity or a definitive clustering ground truth.", "", "## Provenance", "- DATASET-DERIVED: all assignments, agreement metrics, and stability scores.", "- TRANSCRIPT-DERIVED: the supplied clusters and focus on clusters 6 and 7.", "- MODEL PRETRAINED KNOWLEDGE: the use of PCA, nearest-neighbour graphs, Leiden, ARI, NMI, and co-clustering as validation tools follows standard single-cell analysis practice and is not stated in the transcript.", "", "## Output files", "- `results/stability_runs.csv`", "- `results/stability_cell_assignments.csv`", "- `results/stability_cell_scores.csv`", "- `results/stability_cluster_scores.csv`", "- `results/stability_co_cluster_frequency.npy`"]
(OUT / "stability_report.md").write_text("\n".join(report) + "\n")
print("\n".join(report))
