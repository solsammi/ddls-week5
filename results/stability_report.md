# Preliminary clustering stability analysis

## Design
- Three random seeds: 0, 1, 2.
- Three neighbour counts: 10, 15, 30.
- Three Leiden resolutions: 0.5, 0.8, 1.0.
- Total: 27 clustering runs.
- Preprocessing and PCA were performed once; neighbour graphs and Leiden assignments were varied. UMAP was not recomputed because this analysis evaluates cluster assignment stability.

## Run summary
             run  seed  n_neighbors  resolution  n_clusters  ARI_vs_original  NMI_vs_original
cluster_0_10_0_5     0           10         0.5           8         0.940768         0.925157
cluster_0_10_0_8     0           10         0.8          10         0.660454         0.819551
cluster_0_10_1_0     0           10         1.0          10         0.630023         0.808246
cluster_0_15_0_5     0           15         0.5           7         0.950587         0.937387
cluster_0_15_0_8     0           15         0.8          10         0.638145         0.814802
cluster_0_15_1_0     0           15         1.0          10         0.632353         0.813415
cluster_0_30_0_5     0           30         0.5           6         0.910473         0.904153
cluster_0_30_0_8     0           30         0.8           9         0.685212         0.825764
cluster_0_30_1_0     0           30         1.0           9         0.640311         0.818130
cluster_1_10_0_5     1           10         0.5           8         0.949168         0.932798
cluster_1_10_0_8     1           10         0.8          10         0.678470         0.813635
cluster_1_10_1_0     1           10         1.0          10         0.629069         0.812417
cluster_1_15_0_5     1           15         0.5           6         0.902578         0.899966
cluster_1_15_0_8     1           15         0.8          10         0.670671         0.824613
cluster_1_15_1_0     1           15         1.0          10         0.635105         0.814073
cluster_1_30_0_5     1           30         0.5           6         0.912520         0.903561
cluster_1_30_0_8     1           30         0.8           9         0.638816         0.815659
cluster_1_30_1_0     1           30         1.0           9         0.639911         0.817662
cluster_2_10_0_5     2           10         0.5           8         0.943716         0.928805
cluster_2_10_0_8     2           10         0.8          10         0.646803         0.817215
cluster_2_10_1_0     2           10         1.0          10         0.636972         0.814360
cluster_2_15_0_5     2           15         0.5           8         0.954305         0.942484
cluster_2_15_0_8     2           15         0.8          10         0.634640         0.815276
cluster_2_15_1_0     2           15         1.0          10         0.631900         0.810352
cluster_2_30_0_5     2           30         0.5           6         0.912632         0.906140
cluster_2_30_0_8     2           30         0.8           9         0.652576         0.831050
cluster_2_30_1_0     2           30         1.0           9         0.640587         0.819734

## Cluster stability
cluster  n_cells  median_same_original_cocluster_frequency  min_same_original_cocluster_frequency  max_same_original_cocluster_frequency  median_outside_cocluster_frequency
      0     1197                                  0.669763                               0.000836                               0.684535                            0.004633
      1      489                                  0.969110                               0.175850                               0.969110                            0.012480
      2      445                                  0.738155                               0.012513                               0.738155                            0.016129
      3      347                                  1.000000                               1.000000                               1.000000                            0.001228
      4      163                                  0.980338                               0.162551                               0.980338                            0.034263
      5       36                                  1.000000                               1.000000                               1.000000                            0.000584
      6       13                                  0.978395                               0.740741                               0.978395                            0.000207
      7       10                                  0.794239                               0.213992                               0.794239                            0.034325

## Interpretation
- Co-clustering frequency is label-invariant: it measures how often cells are assigned to the same cluster across runs, regardless of numeric cluster names.
- `same_original_cocluster_frequency` asks whether a cell repeatedly stays with cells from its supplied cluster; `other_cluster_cocluster_frequency` identifies possible mixing with other supplied clusters.
- This is preliminary evidence, not proof of biological identity or a definitive clustering ground truth.

## Provenance
- DATASET-DERIVED: all assignments, agreement metrics, and stability scores.
- TRANSCRIPT-DERIVED: the supplied clusters and focus on clusters 6 and 7.
- MODEL PRETRAINED KNOWLEDGE: the use of PCA, nearest-neighbour graphs, Leiden, ARI, NMI, and co-clustering as validation tools follows standard single-cell analysis practice and is not stated in the transcript.

## Output files
- `results/stability_runs.csv`
- `results/stability_cell_assignments.csv`
- `results/stability_cell_scores.csv`
- `results/stability_cluster_scores.csv`
- `results/stability_co_cluster_frequency.npy`
