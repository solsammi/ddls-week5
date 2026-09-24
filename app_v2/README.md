# PBMC Evidence Navigator — version 2

Version 2 of the navigator. Version 1 is kept unchanged in the project root (`app.py`, `README.md`, git tag `app-v1`). Version 2 lives only in this `app_v2/` folder and reads the shared `data/` and `results/` folders in the project root without changing them.

## Changes from version 1

No new analyses were added. Version 2 only surfaces results that already exist in `data/pbmc3k.h5ad` and `results/`.

- **Find a cell:** type a cell ID, or click any dot on the map, to highlight that cell and see its cluster, UMAP position, quality values, per-cell clustering stability (`results/stability_cell_scores.csv`), and raw-count screening result (`results/doublet_raw_cell_metrics.csv`). An unknown ID returns the closest matching IDs.
- **Colour by per-cell clustering stability:** a new "Colour map by" option shows how often each cell stayed with its original cluster across the 27 stability runs.
- **Doublet card read from results:** the doublet-screening card now reads `results/doublet_cluster_summary.csv` instead of fixed text.
- **Updated decision text for clusters 6 and 7:** matches the seminar slides, with confidence stated for each.
- **Faster marker lists:** the Wilcoxon cluster-versus-rest ranking is computed once for all clusters at startup instead of on each first click. Scores are identical to version 1; genes with exactly tied scores may appear in a different order.


A small FastAPI and Scanpy application for exploring the supplied processed PBMC single-cell RNA-seq dataset. The navigator is designed to help Dr. Ravi Menon make two decisions: whether cluster 7 merits another sequencing run, and whether deleting cluster 6 is safe.

## What the navigator does

- Loads `data/pbmc3k.h5ad` once when the application starts.
- Displays the existing UMAP with all eight computational clusters in separate colours.
- Lets you click cluster buttons to highlight one cluster while greying the others.
- Shows cluster quality summaries (`n_genes`, `total_counts`, and `pct_mito`).
- Colours the UMAP by any gene or by a quality metric.
- Shows two distinct gene lists for a selected cluster:
  - **Highest average expression:** genes with the highest mean log-normalised expression within that cluster.
  - **Ranked marker genes:** genes ranked by Wilcoxon comparison of the selected cluster against all other cells.
- Shows owner-supplied marker programmes for T cells, B cells, NK cells, monocytes, dendritic-cell-like cells, and platelet contamination clues.
- Shows a size-stable expression dot plot and table. Dot colour represents mean expression; dot size represents the percentage of cells expressing the gene.
- Shows clustering-stability results from 27 validation runs and a raw-count doublet-screening summary.
- Shows provisional annotations with separate columns for dataset-derived observations, transcript-derived knowledge, and model pretrained biological knowledge.
- Provides a downloadable CSV report for each selected cluster.
- Looks up a single cell by ID or by clicking the map (new in v2).
- Colours the map by per-cell clustering stability (new in v2).

The app is exploratory and evidence-oriented. Cluster numbers are computational labels, not cell-type names. UMAP coordinates are a visual projection: being above or below an axis has no biological meaning.

## Data

The supplied archive is:

```text
ddls-week5-s1-junk-or-signal-dataset.zip
```

It contains `data/pbmc3k.h5ad` and `data/ABOUT_THIS_FILE.txt`. The application can extract the `.h5ad` file from the archive automatically if the extracted file is not present.

The `.h5ad` file contains 2,700 cells and 13,714 genes, with:

- log-normalised expression in `ad.X`
- raw UMI counts in `ad.layers["counts"]`
- cluster labels in `ad.obs["leiden"]`
- quality fields in `ad.obs["n_genes"]`, `ad.obs["total_counts"]`, and `ad.obs["pct_mito"]`
- UMAP coordinates in `ad.obsm["X_umap"]`

## Environment setup

`uv` is required. The environment is already configured for this project; recreate it if necessary:

```bash
uv venv
uv pip install scanpy fastapi 'uvicorn[standard]' python-multipart igraph
```

Run all Python through `uv run`.

## Run the navigator

From the project root (the folder containing `data/` and `results/`), start version 2 on port 8010:

```bash
uv run uvicorn app2:app --app-dir app_v2 --host 127.0.0.1 --port 8010
```

Open [http://127.0.0.1:8010](http://127.0.0.1:8010). Startup takes a little longer than version 1 because marker genes are ranked for all clusters at startup.

To compare with version 1, run it at the same time on port 8000:

```bash
uv run uvicorn app:app --host 127.0.0.1 --port 8000
```

If a port is already in use, find and stop the existing process:

```bash
lsof -nP -iTCP:8010 -sTCP:LISTEN
kill <PID>
```

## API endpoints

- `GET /api/umap` — coordinates, cluster labels, and cell identifiers.
- `GET /api/gene/{gene}` — expression values for a gene across cells.
- `GET /api/quality/{field}` — values for `n_genes`, `total_counts`, or `pct_mito`.
- `GET /api/clusters/{cluster}/genes` — highest average expression and ranked markers.
- `GET /api/clusters/{cluster}/quality` — quality values and summaries.
- `GET /api/clusters/{cluster}/programmes` — owner-supplied programme scores.
- `GET /api/expression-summary/{cluster}` — dot-plot and table data.
- `GET /api/stability` — clustering stability summary.
- `GET /api/stability/cells` — per-cell clustering stability for colouring the map (new in v2).
- `GET /api/cell/{cell_id}` — one cell's cluster, UMAP position, quality values, stability, and raw-count screen (new in v2).
- `GET /api/doublets` — doublet-screening summary for clusters 6 and 7, read from `results/doublet_cluster_summary.csv`.
- `GET /api/annotations` — provenance-separated provisional annotations.
- `GET /api/report/{cluster}.csv` — downloadable selected-cluster report.

## Analysis scripts and results

Analysis code is kept in reviewable Python scripts and outputs are written under `results/`:

- `inspect_dataset.py` — dataset and cluster quality inspection.
- `rank_markers.py` — cluster-versus-rest marker ranking.
- `analysis_report.py` — initial cluster evidence report.
- `annotate_clusters.py` — programme scores and provisional annotation tables.
- `validate_reclustering.py` — independent reclustering validation.
- `stability_analysis.py` — 27-run clustering stability analysis.
- `doublet_investigation.py` — raw-count doublet screening.

## Answer to the owner's question

Cluster 7 should **not currently be prioritised for another sequencing run**: although its 10 cells have high expression depth, the supplied lineage programmes overlap and do not establish one coherent T-cell, B-cell, NK-cell, or monocyte identity across the cluster; confidence is low-to-moderate because the cluster is very small and computationally less stable. Cluster 6 should **not be deleted as dead or empty based on its size or median detected-gene count alone**: its 13 cells are heterogeneous, include a high-depth outlier, and show a dataset-derived `PPBP`/`PF4` signal that the transcript explicitly identifies as a platelet-contamination clue; confidence is moderate for retaining it for review, but its final biological identity remains unresolved. This conclusion separates dataset observations from the transcript’s criteria; broader gene-function interpretations are model pretrained knowledge and are not treated as evidence from the supplied files.

## Provenance and limitations

- **Dataset-derived:** expression, raw counts, quality values, cluster assignments, UMAP, marker rankings, stability metrics, and doublet-screening metrics.
- **Transcript-derived:** the owner’s decision criteria, marker examples, and definitions of coherent, mixed, dead-looking, and empty-droplet-like signals.
- **Model pretrained knowledge:** biological interpretations not explicitly stated by the transcript or file description. These are disclosed separately and should not be confused with direct evidence.
- Small clusters 6 and 7 contain only 13 and 10 cells, respectively.
- Marker rankings and doublet flags are evidence for review, not proof of cell identity or doublet status.
- The clustering stability analysis used three seeds, three neighbour settings, and three resolutions; it is preliminary rather than definitive.
