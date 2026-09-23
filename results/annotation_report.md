# Annotation report

## Provenance
- DATASET-DERIVED: expression, programme scores, quality values, cluster-versus-rest marker rankings.
- TRANSCRIPT-DERIVED: marker examples and decision definitions supplied by Dr. Menon.
- MODEL PRETRAINED KNOWLEDGE: disclosed separately and not treated as evidence from the supplied files.

## Results
### Cluster 0: Not assigned
**Confidence:** Not assessed
**Dataset-derived:** Programme scores and ranked markers are available in the output tables.
**Transcript-derived:** Owner requested T, B, NK, or monocyte identities.
**Model pretrained knowledge:** No model-only annotation used.

### Cluster 1: Not assigned
**Confidence:** Not assessed
**Dataset-derived:** Programme scores and ranked markers are available in the output tables.
**Transcript-derived:** Owner requested T, B, NK, or monocyte identities.
**Model pretrained knowledge:** No model-only annotation used.

### Cluster 2: Not assigned
**Confidence:** Not assessed
**Dataset-derived:** Programme scores and ranked markers are available in the output tables.
**Transcript-derived:** Owner requested T, B, NK, or monocyte identities.
**Model pretrained knowledge:** No model-only annotation used.

### Cluster 3: Not assigned
**Confidence:** Not assessed
**Dataset-derived:** Programme scores and ranked markers are available in the output tables.
**Transcript-derived:** Owner requested T, B, NK, or monocyte identities.
**Model pretrained knowledge:** No model-only annotation used.

### Cluster 4: Not assigned
**Confidence:** Not assessed
**Dataset-derived:** Programme scores and ranked markers are available in the output tables.
**Transcript-derived:** Owner requested T, B, NK, or monocyte identities.
**Model pretrained knowledge:** No model-only annotation used.

### Cluster 5: Not assigned
**Confidence:** Not assessed
**Dataset-derived:** Programme scores and ranked markers are available in the output tables.
**Transcript-derived:** Owner requested T, B, NK, or monocyte identities.
**Model pretrained knowledge:** No model-only annotation used.

### Cluster 6: Platelet-contamination clue; unresolved
**Confidence:** Moderate
**Dataset-derived:** Cluster-versus-rest markers include PPBP and PF4, and programme scores are recorded in the output table.
**Transcript-derived:** The owner explicitly supplied PPBP and PF4 as platelet-contamination clues; the owner said low gene count alone is insufficient to call cells dead/empty.
**Model pretrained knowledge:** No additional platelet-gene annotation is used in this label. Any broader platelet/megakaryocyte interpretation would be model knowledge, not transcript evidence.

### Cluster 7: Unresolved/mixed for requested lineages
**Confidence:** Low-to-moderate
**Dataset-derived:** Cluster-versus-rest ranking is dominated by the listed genes; programme scores are recorded in the output table. The cluster has 10 cells.
**Transcript-derived:** A clear identity requires several same-lineage genes across most cells; incompatible or uninformative signals argue against another sequencing run.
**Model pretrained knowledge:** The possible cell-cycle interpretation of genes such as STMN1, PCNA, TYMS, ZWINT, and PTTG1 is model pretrained knowledge and is not used as the annotation basis.

## Operational decisions
- Cluster 6: do not delete based on size or median n_genes alone. The supplied transcript explicitly identifies PPBP and PF4 as platelet-contamination clues; review the cells as a possible specific signal rather than generic junk. Confidence: moderate.
- Cluster 7: do not prioritise for another sequencing run on current evidence. The requested lineage programmes do not establish one coherent identity across 10 cells. Confidence: low-to-moderate.

## Files
- `results/annotation_programme_scores.csv`
- `results/annotation_ranked_markers.csv`
- `results/annotation_conclusions.csv`
