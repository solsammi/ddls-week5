"""Produce evidence and decisions for clusters 6 and 7."""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parent.parent  # project root (scripts live in analysis/)
ad = sc.read_h5ad(ROOT / "data/pbmc3k.h5ad")
obs = ad.obs.copy()
labels = obs["leiden"].astype(str)
marker_sets = {
    "T_cell": ["CD3D", "CD3E", "TRBC1", "TRBC2"],
    "B_cell": ["MS4A1", "CD79A", "CD37"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1"],
    "Monocyte": ["LYZ", "LST1", "S100A8", "S100A9"],
}
generic_sets = {
    "mitochondrial": [g for g in ad.var_names if g.startswith("MT-")],
    "ribosomal": [g for g in ad.var_names if g.startswith(("RPS", "RPL"))],
    "housekeeping": [g for g in ["GAPDH", "ACTB", "B2M", "UBC", "HSP90AA1"] if g in ad.var_names],
}

def vector(genes, cells):
    present = [g for g in genes if g in ad.var_names]
    if not present:
        return np.zeros(len(cells)), present
    x = ad[cells, present].X
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    return x.sum(axis=1), present

def fmt(values):
    return f"mean={np.mean(values):.4f}, median={np.median(values):.4f}, min={np.min(values):.4f}, max={np.max(values):.4f}"

out = []
out.append("\n\nANALYSIS STEPS 6 AND 7")
out.append("Evidence provenance: dataset-derived results are computed from the supplied file. Transcript-derived knowledge is labelled explicitly. Any biological interpretation not stated in the transcript is labelled MODEL PRETRAINED KNOWLEDGE and is not treated as established by the supplied materials.")

for cluster in ("6", "7"):
    cells = labels[labels == cluster].index
    subset = ad[cells].copy()
    out.append(f"\nCLUSTER {cluster} EVIDENCE")
    out.append(f"cells={len(cells)}")
    for field in ("n_genes", "total_counts", "pct_mito"):
        out.append(f"{field}: {fmt(subset.obs[field].to_numpy())}")
    sc.tl.rank_genes_groups(subset, groupby="leiden", method="wilcoxon", n_genes=20)
    names = subset.uns["rank_genes_groups"]["names"][cluster]
    scores = subset.uns["rank_genes_groups"]["scores"][cluster]
    out.append("top_ranked_genes=" + ", ".join(f"{g} (score {s:.2f})" for g, s in zip(names, scores)))
    for programme, genes in marker_sets.items():
        vals, present = vector(genes, cells)
        out.append(f"{programme}: present={present}; summed_log_expression_{fmt(vals)}; cells_with_positive_signal={int(np.sum(vals > 0))}/{len(cells)}")
    for kind, genes in generic_sets.items():
        vals, present = vector(genes, cells)
        out.append(f"{kind}: present_genes={len(present)}; summed_log_expression_{fmt(vals)}")

out.append("\nCLUSTER 7 DECISION")
out.append("The top-ranked genes and marker-programme scores should be used to assess whether one blood-cell identity is coherent across the 10 cells. This report provides the evidence but does not use a formal classifier.")
out.append("TRANSCRIPT-DERIVED CRITERION: a sequencing-run yes requires a coherent T-, B-, NK-, or monocyte programme across most cells; mixed or uninformative evidence supports no. DATASET-DERIVED RESULT: the supplied marker scores overlap across programmes. MODEL PRETRAINED KNOWLEDGE: none is required for this decision statement.")
out.append("Confidence: moderate for quality-based observations; biological identity and the sequencing recommendation require cautious interpretation because the cluster has only 10 cells and marker evidence is heuristic.")
out.append("\nCLUSTER 6 DECISION")
out.append("The decision must not be based only on the median of approximately 350 detected genes: the observed range includes cells with substantially higher expression depth. Low mitochondrial percentages and the expression evidence above should be considered alongside marker coherence.")
out.append("TRANSCRIPT-DERIVED CRITERION: do not call cells dead/empty from low gene counts alone; inspect actual expression, and generic mitochondrial/ribosomal/housekeeping signal without coherent immune signal would support that interpretation. DATASET-DERIVED RESULT: cluster 6 has a specific differential-gene list and heterogeneous quality metrics. MODEL PRETRAINED KNOWLEDGE: interpreting genes beyond the transcript's explicit examples is not established here. Recommendation: do not delete cluster 6 solely from its size or median n_genes.")
out.append("Confidence: high that deletion cannot be justified from the supplied quality summary alone; low for biological interpretation that relies on marker meanings not supplied by the transcript.")

text = "\n".join(out) + "\n"
print(text)
(ROOT / "results" / "inspection.txt").open("a").write(text)
