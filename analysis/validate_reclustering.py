"""Validate the supplied Leiden clustering with an independent Scanpy workflow."""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parent.parent  # project root (scripts live in analysis/)
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
ad = sc.read_h5ad(ROOT / "data/pbmc3k.h5ad")
original = ad.obs["leiden"].astype(str).to_numpy()

# Independent validation: use raw counts, normalize/log-transform, select HVGs,
# PCA, neighbors, UMAP, and Leiden at several resolutions.
if "counts" in ad.layers:
    ad.X = ad.layers["counts"].copy()
sc.pp.normalize_total(ad, target_sum=1e4)
sc.pp.log1p(ad)
sc.pp.highly_variable_genes(ad, n_top_genes=2000, flavor="seurat", subset=True)
sc.pp.scale(ad, max_value=10)
sc.tl.pca(ad, n_comps=40, svd_solver="arpack", random_state=0)
sc.pp.neighbors(ad, n_neighbors=15, n_pcs=30, random_state=0)
sc.tl.umap(ad, random_state=0)

rows = []
for resolution in (0.3, 0.5, 0.8, 1.0):
    key = f"recluster_{str(resolution).replace('.', '_')}"
    sc.tl.leiden(ad, resolution=resolution, key_added=key, random_state=0, flavor="igraph", n_iterations=2)
    labels = ad.obs[key].astype(str).to_numpy()
    contingency = pd.crosstab(pd.Series(original, name="original"), pd.Series(labels, name="recluster"))
    agreement = np.mean(original == labels) if len(set(original)) == len(set(labels)) else np.nan
    # Best one-to-one-independent overlap: for each recluster, largest original fraction.
    purity = np.mean([contingency[col].max() / contingency[col].sum() for col in contingency.columns])
    rows.append({"resolution": resolution, "n_reclustered": len(set(labels)), "same_label_fraction_unadjusted": agreement, "mean_recluster_cluster_purity": purity})
    contingency.to_csv(OUT / f"reclustering_contingency_{str(resolution).replace('.', '_')}.csv")

summary = pd.DataFrame(rows)
summary.to_csv(OUT / "reclustering_summary.csv", index=False)
ad.obs[["leiden", "recluster_0_3", "recluster_0_5", "recluster_0_8", "recluster_1_0"]].to_csv(OUT / "reclustering_cell_assignments.csv")

# A reproducible report with provenance and limitations.
report = ["# Independent reclustering validation", "", "## Method", "- DATASET-DERIVED: started from the supplied raw UMI layer, normalised to 10,000 counts per cell, log-transformed, selected 2,000 highly variable genes, scaled, computed 40-component PCA, built a 15-nearest-neighbour graph using 30 PCs, computed a fresh UMAP, and ran Leiden at four resolutions.", "- The supplied clustering was not overwritten; original labels remain in `leiden`.", "- This is a validation sensitivity analysis, not proof that one clustering is biologically correct.", "", "## Results", summary.to_string(index=False), "", "## Interpretation", "- Compare the contingency tables and cell assignments rather than relying on raw label equality: Leiden labels are arbitrary names, so identical numeric labels are not required for agreement.", "- A cluster is more robust when its cells remain together across resolutions and when reclustering produces a similar broad partition.", "- Small clusters, especially clusters 6 and 7, are expected to be sensitive to neighbours, HVG choices, resolution, and random seed.", "", "## Provenance", "- TRANSCRIPT-DERIVED: the owner says clusters 0–7 were already computed and the immediate decisions concern clusters 6 and 7.", "- MODEL PRETRAINED KNOWLEDGE: the methodological description above follows standard single-cell analysis practice; it is not stated in the transcript. It is disclosed separately and is not evidence of biological identity.", "- No cell-type labels are assigned from this validation alone.", "", "## Output files", "- `results/reclustering_summary.csv`", "- `results/reclustering_cell_assignments.csv`", "- `results/reclustering_contingency_*.csv`"]
(OUT / "reclustering_report.md").write_text("\n".join(report) + "\n")
print("\n".join(report))
