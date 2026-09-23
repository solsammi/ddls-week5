# Specification

## Problem and decision

Dr. Ravi Menon has a processed PBMC single-cell RNA-seq dataset containing 2,700 cells assigned to clusters 0–7. He needs a clear decision about whether to prioritise cluster 7 for the next sequencing run, and whether deleting cluster 6 is safe.

The immediate decision is:

- **Cluster 7:** Decide yes or no on asking the facility to prioritise this cluster for the next sequencing run. A yes requires a clear, coherent likely blood-cell identity; a no is appropriate if the genes are mixed or uninformative.
- **Cluster 6:** Decide whether its 13 cells are genuinely empty/dead-looking or instead express a consistent immune programme, before deleting them.

## Clusters

There are 8 clusters, numbered 0 through 7, covering all 2,700 cells. The transcript reports these cluster-level `genes/cell` values, but does not establish whether they are means, medians, or another statistic:

- Cluster 0: 809
- Cluster 1: 850
- Cluster 2: 828
- Cluster 3: 673
- Cluster 4: 1,263
- Cluster 5: 1,570
- Cluster 6: 350; 13 cells; summary mitochondrial fraction 1.6%
- Cluster 7: 2,363; 10 cells

Only clusters 6 and 7 are in scope for the requested decision. The owner does not want a full audit of clusters 0–5.

## Files

- `ddls-week5-interview.md`: the interview transcript and the source of the owner's goals, definitions, requested evidence, and known uncertainties.
- `ddls-week5-s1-junk-or-signal-dataset.zip`: the supplied data archive.
- Inside the archive, `data/pbmc3k.h5ad`: processed single-cell RNA-seq data for 2,700 human PBMCs and 13,714 genes.
- Inside the archive, `data/ABOUT_THIS_FILE.txt`: file description and access instructions. It states that `ad.X` contains log-normalised expression, `ad.layers["counts"]` raw UMI counts, `ad.var_names` gene symbols, `ad.obs["leiden"]` cluster labels, `ad.obs["n_genes"]`, `ad.obs["total_counts"]`, and `ad.obs["pct_mito"]` per-cell quality fields, and `ad.obsm["X_umap"]` 2-D map coordinates.
- `.env`: local environment/secrets file; it is excluded by `.gitignore` and must not be committed.

## Exact owner ask and evidence required

### Cluster 7

The owner asks for:

1. The genes highlighted on the map for cluster 7.
2. The top expressed/defining genes for the 10 cells in cluster 7, including names and expression values.
3. Expression across those 10 cells, so consistency can be assessed rather than relying on cluster size or gene count.
4. A likely identity among T cell, B cell, NK cell, or monocyte, or a finding that the cluster is mixed/uninformative.
5. A yes/no recommendation on whether to prioritise cluster 7 for the next sequencing run.

A clear identity means several genes indicating the same programme across most cells, not one isolated gene. Relevant transcript-provided marker clues are:

- T cells: `CD3D`, `CD3E`, `TRBC1/2`
- B cells: `MS4A1`, `CD79A`, `CD37`
- NK cells: `NKG7`, `GNLY`, `KLRD1`
- Monocytes: `LYZ`, `LST1`; inflammatory monocytes may show `S100A8/A9`
- Dendritic-cell-like clues: `FCER1A`, `CST3`
- Platelet contamination clues: `PPBP`, `PF4`

Strong incompatible programmes, such as convincing T-cell and monocyte signals together in the same cells, count against a coherent identity. Mostly housekeeping, ribosomal, stress, or mitochondrial genes count as generic/uninformative signals.

### Cluster 6

The owner asks whether the 13 cells are:

- dead-cell-like,
- empty-droplet-like, or
- expressing something specific and consistently immune-related.

The evidence must inspect the actual expression, not infer junk solely from approximately 350 genes per cell. Dead/empty-looking evidence would be very little overall expression, mostly mitochondrial, ribosomal, or housekeeping genes, and no coherent immune-cell pattern. A conspicuously high mitochondrial fraction relative to other cells could support a dying-cell interpretation, but there is no fixed threshold supplied. The reported 1.6% mitochondrial fraction alone is not enough to call the cells junk.

Examples supplied by the owner include mitochondrial genes (`MT-ND1`, `MT-ND2`, `MT-CO1`, `MT-CO2`, `MT-ATP6`, `MT-CYB`), ribosomal genes (`RPLP0`, `RPS3`, `RPL13`, `RPS18`), and housekeeping genes (`GAPDH`, `ACTB`, `B2M`, `UBC`, `HSP90AA1`).

## Definition of done

Done means the owner can give the facility a defensible yes-or-no answer about spending the next sequencing run on cluster 7, supported by the genes highlighted on the map, top genes, their expression across the 10 cells, a likely identity or a mixed/uninformative conclusion, and an explicit confidence level matching that claim.

Done also means the owner has an evidence-based answer about whether deleting cluster 6 is safe, based on its actual expression across all 13 cells and its quality signals—not just the 350-gene summary. Any limitations, unknowns, and confidence must be stated. No answer about either cluster is complete without its confidence.
