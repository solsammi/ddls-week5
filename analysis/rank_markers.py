"""Rank cluster markers against all other cells and report clusters 6 and 7."""
from pathlib import Path
import scanpy as sc
import numpy as np

ROOT = Path(__file__).resolve().parent.parent  # project root (scripts live in analysis/)
ad = sc.read_h5ad(ROOT / "data/pbmc3k.h5ad")
ad.obs["cluster_for_ranking"] = ad.obs["leiden"].astype(str)
sc.tl.rank_genes_groups(ad, groupby="cluster_for_ranking", method="wilcoxon", n_genes=30, key_added="cluster_markers")

out = ["\n\nCORRECTED MARKER RANKING (CLUSTER VS ALL OTHER CELLS)", "PROVENANCE: marker rankings and expression summaries are DATASET-DERIVED. Marker meanings explicitly supplied by the TRANSCRIPT are limited to its listed clues. Any additional biological meaning is MODEL PRETRAINED KNOWLEDGE, was not retrieved from the supplied files, and must not be presented as established evidence."]
for cluster in ("6", "7"):
    names = ad.uns["cluster_markers"]["names"][cluster]
    scores = ad.uns["cluster_markers"]["scores"][cluster]
    pvals = ad.uns["cluster_markers"]["pvals_adj"][cluster]
    logfc = ad.uns["cluster_markers"]["logfoldchanges"][cluster]
    out.append(f"\nCLUSTER {cluster} TOP 30 MARKERS")
    out.append("rank\tgene\tscore\tadjusted_p\tlog2_fold_change")
    for i, (gene, score, pval, fc) in enumerate(zip(names, scores, pvals, logfc), 1):
        out.append(f"{i}\t{gene}\t{score:.3f}\t{pval:.3g}\t{fc:.3f}")

    cells = ad.obs["cluster_for_ranking"] == cluster
    out.append("\nTRANSCRIPT-DERIVED MARKER EVIDENCE")
    marker_sets = {
        "T_cell": ["CD3D", "CD3E", "TRBC1", "TRBC2"],
        "B_cell": ["MS4A1", "CD79A", "CD37"],
        "NK_cell": ["NKG7", "GNLY", "KLRD1"],
        "Monocyte": ["LYZ", "LST1", "S100A8", "S100A9"],
    }
    for label, genes in marker_sets.items():
        present = [g for g in genes if g in ad.var_names]
        x = ad[cells, present].X
        x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
        signal = x.sum(axis=1)
        out.append(f"{label}: genes={present}; positive_cells={int((signal > 0).sum())}/{len(signal)}; mean_sum={signal.mean():.4f}; median_sum={np.median(signal):.4f}")

text = "\n".join(out) + "\n"
print(text)
(ROOT / "results" / "inspection.txt").open("a").write(text)
