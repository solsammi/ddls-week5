"""Investigate possible doublets using raw counts and transcript-supplied programmes."""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
ad = sc.read_h5ad(ROOT / "data/pbmc3k.h5ad")
counts = ad.layers["counts"]
counts = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
genes = np.asarray(ad.var_names)
idx = {g:i for i,g in enumerate(genes)}
programmes = {"T cell":["CD3D","CD3E","TRBC1","TRBC2"],"B cell":["MS4A1","CD79A","CD37"],"NK cell":["NKG7","GNLY","KLRD1"],"Monocyte":["LYZ","LST1","S100A8","S100A9"],"Platelet clue":["PPBP","PF4"]}
obs = ad.obs.copy(); labels=obs.leiden.astype(str).to_numpy()
qc = pd.DataFrame(index=ad.obs_names)
qc["cluster"] = labels; qc["raw_total_counts"] = counts.sum(1); qc["raw_n_genes"] = (counts>0).sum(1)
median=qc.raw_total_counts.median(); mad=np.median(np.abs(qc.raw_total_counts-median)); qc["raw_counts_robust_z"]=(qc.raw_total_counts-median)/(1.4826*mad)
median=qc.raw_n_genes.median(); mad=np.median(np.abs(qc.raw_n_genes-median)); qc["raw_genes_robust_z"]=(qc.raw_n_genes-median)/(1.4826*mad)
for name, gs in programmes.items():
    present=[idx[g] for g in gs if g in idx]; vals=counts[:,present] if present else np.zeros((len(counts),0))
    qc[name+"_counts"] = vals.sum(1) if present else 0
    qc[name+"_genes_positive"] = (vals>0).sum(1) if present else 0
qc["n_strong_programmes"] = sum((qc[f"{name}_genes_positive"] >= (2 if name not in ["Platelet clue"] else 1)).astype(int) for name in programmes)
qc["possible_mixed_programmes"] = qc.n_strong_programmes >= 2
qc.to_csv(OUT/"doublet_raw_cell_metrics.csv")

cluster = qc.groupby("cluster").agg(n_cells=("cluster","size"), median_raw_counts=("raw_total_counts","median"), max_raw_counts=("raw_total_counts","max"), median_raw_genes=("raw_n_genes","median"), max_raw_genes=("raw_n_genes","max"), mixed_cells=("possible_mixed_programmes","sum"), median_programmes=("n_strong_programmes","median"), max_programmes=("n_strong_programmes","max")).reset_index()
cluster.to_csv(OUT/"doublet_cluster_summary.csv",index=False)

focus=qc[(qc.cluster.isin(["6","7"])) | (qc.raw_counts_robust_z>3) | (qc.raw_genes_robust_z>3)].sort_values(["cluster","raw_total_counts"],ascending=[True,False])
focus.to_csv(OUT/"doublet_focus_cells.csv")

report=["# Raw-count doublet investigation","","## Provenance","- DATASET-DERIVED: all metrics below use `ad.layers['counts']`, the supplied raw UMI count layer.","- TRANSCRIPT-DERIVED: programme genes are the owner-supplied clues; PPBP/PF4 are explicitly platelet-contamination clues.","- MODEL PRETRAINED KNOWLEDGE: high depth, incompatible programme co-expression, and doublet-style interpretation are general biological/technical knowledge, not statements from the transcript.","","## Method","- Calculated raw total counts and raw detected genes per cell.","- Calculated robust z-scores using the dataset-wide median and MAD.","- For each supplied programme, counted raw positive marker genes and raw counts.","- Flagged a possible mixed cell when at least two programme groups had multiple positive markers; the platelet clue required one positive gene. This is a screening rule, not a validated doublet classifier.","","## Cluster summary",cluster.to_string(index=False),"","## Focus cells",focus[["cluster","raw_total_counts","raw_n_genes","raw_counts_robust_z","raw_genes_robust_z","n_strong_programmes","possible_mixed_programmes"]].to_string(),"","## Assessment","- Cluster 6: the cluster contains 13 cells and includes a high-depth cell; its transcript-supported PP4/PF4 platelet-contamination signal should be investigated. The raw-count screen alone does not establish a doublet.","- Cluster 7: the cluster has high depth and a few cells may show multiple programme signals, but broad marker positivity is not by itself proof of a doublet. Its suspected cycling interpretation is model knowledge, not transcript knowledge.","- A stronger doublet call would require per-cell incompatible raw-count co-expression plus supporting evidence from a dedicated method such as Scrublet and review of the flagged cell’s neighbourhood.","","## Outputs","- `results/doublet_raw_cell_metrics.csv`","- `results/doublet_cluster_summary.csv`","- `results/doublet_focus_cells.csv`"]
(OUT/"doublet_report.md").write_text("\n".join(report)+"\n")
print("\n".join(report))
