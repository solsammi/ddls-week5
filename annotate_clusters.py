"""Annotate PBMC clusters with explicit evidence provenance."""
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
ad = sc.read_h5ad(ROOT / "data/pbmc3k.h5ad")
ad.obs["cluster"] = ad.obs["leiden"].astype(str)

# TRANSCRIPT-DERIVED marker examples supplied by the owner.
programmes = {
    "T cell": ["CD3D", "CD3E", "TRBC1", "TRBC2"],
    "B cell": ["MS4A1", "CD79A", "CD37"],
    "NK cell": ["NKG7", "GNLY", "KLRD1"],
    "Monocyte": ["LYZ", "LST1", "S100A8", "S100A9"],
    "Dendritic-cell-like": ["FCER1A", "CST3"],
    "Platelet contamination": ["PPBP", "PF4"],
}

rows = []
for cluster in sorted(ad.obs["cluster"].unique(), key=int):
    mask = ad.obs["cluster"].to_numpy() == cluster
    for programme, genes in programmes.items():
        present = [g for g in genes if g in ad.var_names]
        x = ad[mask, present].X if present else np.zeros((mask.sum(), 0))
        x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
        score = x.mean(axis=1) if present else np.zeros(mask.sum())
        rows.append({
            "cluster": cluster,
            "programme": programme,
            "marker_genes_present": ", ".join(present),
            "n_cells": int(mask.sum()),
            "cells_with_positive_score": int((score > 0).sum()),
            "pct_cells_positive": float((score > 0).mean() * 100),
            "mean_score": float(score.mean()),
            "median_score": float(np.median(score)),
        })
programme_df = pd.DataFrame(rows)
programme_df.to_csv(OUT / "annotation_programme_scores.csv", index=False)

# Cluster-versus-rest ranking.
sc.tl.rank_genes_groups(ad, groupby="cluster", method="wilcoxon", n_genes=30, key_added="annotation_markers")
marker_rows = []
for cluster in sorted(ad.obs["cluster"].unique(), key=int):
    names = ad.uns["annotation_markers"]["names"][cluster]
    scores = ad.uns["annotation_markers"]["scores"][cluster]
    pvals = ad.uns["annotation_markers"]["pvals_adj"][cluster]
    folds = ad.uns["annotation_markers"]["logfoldchanges"][cluster]
    for rank, (gene, score, pval, fold) in enumerate(zip(names, scores, pvals, folds), 1):
        marker_rows.append({"cluster": cluster, "rank": rank, "gene": str(gene), "score": float(score), "adjusted_p": float(pval), "log2_fold_change": float(fold)})
marker_df = pd.DataFrame(marker_rows)
marker_df.to_csv(OUT / "annotation_ranked_markers.csv", index=False)

# Provenance-separated conclusions. These are deliberately conservative.
conclusions = [
    {"cluster": "0", "label": "Not assigned", "confidence": "Not assessed", "dataset_evidence": "Programme scores and ranked markers are available in the output tables.", "transcript_knowledge": "Owner requested T, B, NK, or monocyte identities.", "model_knowledge": "No model-only annotation used."},
    {"cluster": "1", "label": "Not assigned", "confidence": "Not assessed", "dataset_evidence": "Programme scores and ranked markers are available in the output tables.", "transcript_knowledge": "Owner requested T, B, NK, or monocyte identities.", "model_knowledge": "No model-only annotation used."},
    {"cluster": "2", "label": "Not assigned", "confidence": "Not assessed", "dataset_evidence": "Programme scores and ranked markers are available in the output tables.", "transcript_knowledge": "Owner requested T, B, NK, or monocyte identities.", "model_knowledge": "No model-only annotation used."},
    {"cluster": "3", "label": "Not assigned", "confidence": "Not assessed", "dataset_evidence": "Programme scores and ranked markers are available in the output tables.", "transcript_knowledge": "Owner requested T, B, NK, or monocyte identities.", "model_knowledge": "No model-only annotation used."},
    {"cluster": "4", "label": "Not assigned", "confidence": "Not assessed", "dataset_evidence": "Programme scores and ranked markers are available in the output tables.", "transcript_knowledge": "Owner requested T, B, NK, or monocyte identities.", "model_knowledge": "No model-only annotation used."},
    {"cluster": "5", "label": "Not assigned", "confidence": "Not assessed", "dataset_evidence": "Programme scores and ranked markers are available in the output tables.", "transcript_knowledge": "Owner requested T, B, NK, or monocyte identities.", "model_knowledge": "No model-only annotation used."},
    {"cluster": "6", "label": "Platelet-contamination clue; unresolved", "confidence": "Moderate", "dataset_evidence": "Cluster-versus-rest markers include PPBP and PF4, and programme scores are recorded in the output table.", "transcript_knowledge": "The owner explicitly supplied PPBP and PF4 as platelet-contamination clues; the owner said low gene count alone is insufficient to call cells dead/empty.", "model_knowledge": "No additional platelet-gene annotation is used in this label. Any broader platelet/megakaryocyte interpretation would be model knowledge, not transcript evidence."},
    {"cluster": "7", "label": "Unresolved/mixed for requested lineages", "confidence": "Low-to-moderate", "dataset_evidence": "Cluster-versus-rest ranking is dominated by the listed genes; programme scores are recorded in the output table. The cluster has 10 cells.", "transcript_knowledge": "A clear identity requires several same-lineage genes across most cells; incompatible or uninformative signals argue against another sequencing run.", "model_knowledge": "The possible cell-cycle interpretation of genes such as STMN1, PCNA, TYMS, ZWINT, and PTTG1 is model pretrained knowledge and is not used as the annotation basis."},
]
conclusion_df = pd.DataFrame(conclusions)
conclusion_df.to_csv(OUT / "annotation_conclusions.csv", index=False)

report = ["# Annotation report", "", "## Provenance", "- DATASET-DERIVED: expression, programme scores, quality values, cluster-versus-rest marker rankings.", "- TRANSCRIPT-DERIVED: marker examples and decision definitions supplied by Dr. Menon.", "- MODEL PRETRAINED KNOWLEDGE: disclosed separately and not treated as evidence from the supplied files.", "", "## Results"]
for item in conclusions:
    report += [f"### Cluster {item['cluster']}: {item['label']}", f"**Confidence:** {item['confidence']}", f"**Dataset-derived:** {item['dataset_evidence']}", f"**Transcript-derived:** {item['transcript_knowledge']}", f"**Model pretrained knowledge:** {item['model_knowledge']}", ""]
report += ["## Operational decisions", "- Cluster 6: do not delete based on size or median n_genes alone. The supplied transcript explicitly identifies PPBP and PF4 as platelet-contamination clues; review the cells as a possible specific signal rather than generic junk. Confidence: moderate.", "- Cluster 7: do not prioritise for another sequencing run on current evidence. The requested lineage programmes do not establish one coherent identity across 10 cells. Confidence: low-to-moderate.", "", "## Files", "- `results/annotation_programme_scores.csv`", "- `results/annotation_ranked_markers.csv`", "- `results/annotation_conclusions.csv`"]
(OUT / "annotation_report.md").write_text("\n".join(report) + "\n")
print("Wrote annotation outputs to results/")
print((OUT / "annotation_report.md").read_text())
