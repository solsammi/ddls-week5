# Independent reclustering validation

## Method
- DATASET-DERIVED: started from the supplied raw UMI layer, normalised to 10,000 counts per cell, log-transformed, selected 2,000 highly variable genes, scaled, computed 40-component PCA, built a 15-nearest-neighbour graph using 30 PCs, computed a fresh UMAP, and ran Leiden at four resolutions.
- The supplied clustering was not overwritten; original labels remain in `leiden`.
- This is a validation sensitivity analysis, not proof that one clustering is biologically correct.

## Results
 resolution  n_reclustered  same_label_fraction_unadjusted  mean_recluster_cluster_purity
        0.3              6                             NaN                       0.945980
        0.5              7                             NaN                       0.981208
        0.8             10                             NaN                       0.966907
        1.0             10                             NaN                       0.959105

## Interpretation
- Compare the contingency tables and cell assignments rather than relying on raw label equality: Leiden labels are arbitrary names, so identical numeric labels are not required for agreement.
- A cluster is more robust when its cells remain together across resolutions and when reclustering produces a similar broad partition.
- Small clusters, especially clusters 6 and 7, are expected to be sensitive to neighbours, HVG choices, resolution, and random seed.

## Provenance
- TRANSCRIPT-DERIVED: the owner says clusters 0–7 were already computed and the immediate decisions concern clusters 6 and 7.
- MODEL PRETRAINED KNOWLEDGE: the methodological description above follows standard single-cell analysis practice; it is not stated in the transcript. It is disclosed separately and is not evidence of biological identity.
- No cell-type labels are assigned from this validation alone.

## Output files
- `results/reclustering_summary.csv`
- `results/reclustering_cell_assignments.csv`
- `results/reclustering_contingency_*.csv`
