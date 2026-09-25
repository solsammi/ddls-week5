# PBMC Evidence Navigator: junk or signal?

DDLS 2026, week 5. Presented by Sammi Baudin. Data owner: Dr. Ravi Menon, immunologist.

A FastAPI and Scanpy app for exploring a processed PBMC single-cell RNA-seq dataset (2,700 human peripheral blood mononuclear cells, 13,714 genes). It was built to help the data owner make two decisions about two tiny clusters, and it also works as the seminar presentation.

## The owner's question

The core facility returned the cells already grouped into clusters 0–7, together with a 2-D map (UMAP). Before making a figure, the owner has to decide:

- **Cluster 6** (13 cells, ~350 genes per cell): is it safe to delete as dead cells or empty droplets?
- **Cluster 7** (10 cells, ~2,363 genes per cell): should the next sequencing run go to this cluster? The answer is yes only if the cluster has a clear, likely blood-cell identity (T, B, NK or monocyte).

## The answer

| Cluster | Answer | Confidence |
|---|---|---|
| 6 | **Platelet contamination, not dead or empty cells.** All 13 cells express `PPBP` and `PF4` in raw counts, which is the owner's platelet-contamination clue, and the mitochondrial share is low (0.7–3.2%). Deleting it as dead or empty junk would be wrong; removing it as platelet contamination is the owner's call. | Moderate–high that it is platelet contamination rather than dead or empty material |
| 7 | **Don't prioritise the next run: probably cycling cells.** Its top markers include `STMN1`, `PCNA`, `TYMS` and `PTTG1`, which suggests cycling (dividing) cells, not one of the T, B, NK or monocyte identities the owner asked for. Those programmes overlap in the same cells. The cluster is an 8-cell core plus 2 loosely attached cells. | Low–moderate for the recommendation and for the cycling reading; moderate for the 8 + 2 structure |

Both answers rest on very few cells. See [Limits](#limits).

## Project layout

```text
app.py                 The evidence navigator (FastAPI + Scanpy + Plotly)
analysis/              Analysis scripts (write their outputs to results/)
results/               Analysis outputs read by the app (CSV, Markdown, NPY)
data/pbmc3k.h5ad       The processed dataset
ddls-week5-s1-junk-or-signal-dataset.zip   The dataset archive as supplied
week5_slides.html      Seminar slide deck (self-contained HTML)
ddls-week5-interview.md  Transcript of the interview with the data owner
AGENTS.md, spec.md     Working rules and problem specification
```

## Setup

The project uses [`uv`](https://docs.astral.sh/uv/). Create the environment and install the packages:

```bash
uv venv
uv pip install scanpy fastapi 'uvicorn[standard]' python-multipart igraph
```

Run all Python through `uv run`.

## Run the app

From the project folder:

```bash
uv run uvicorn app:app --host 127.0.0.1 --port 8000
```

Open <http://127.0.0.1:8000>. At startup the app ranks marker genes for all clusters, so it takes a few seconds to be ready. If `data/pbmc3k.h5ad` is missing, the app extracts it from the dataset archive.

The page loads Tailwind and Plotly from online CDNs, so it needs an internet connection. If port 8000 is already in use, find and stop the old process with `lsof -nP -iTCP:8000 -sTCP:LISTEN`, then `kill <PID>`.

### What the app shows

The header shows the talk title and names. A **Larger text** button enlarges everything for a projector. Two cards at the top give the answer for clusters 6 and 7. The five tabs follow the order of the talk:

1. **Data & question**: what data was supplied, the owner's two decisions, their definitions of "coherent", "mixed" and "junk", and their marker clues.
2. **Explore**:
   - the UMAP, coloured by cluster, by any gene, by a quality metric or by per-cell clustering stability;
   - clicking a cluster keeps the current colouring and greys out the other clusters;
   - **Find a cell**: type an ID or click a dot to see that cell's quality values, stability and raw-count screening result;
   - lasso-select cells for a summary;
   - a cluster card with the highest-expressed and ranked marker genes (click a gene to colour the map by it), the owner's marker programmes, and a CSV download;
   - a dot plot of the owner's marker genes across all clusters, where you can add more genes.
   - provisional annotations for all clusters (6 and 7 first), with dataset, transcript and model-knowledge evidence in separate columns.
3. **Clusters 6 & 7**: whether either cluster should be merged, split or binned, with the evidence, and a side-by-side comparison of any two clusters.
4. **Methods & validation**: why each analysis was done and what it showed, the independent reclustering, clustering stability, and the raw-count doublet screen.
5. **Answers, limits & AI use**: the answers, the limits, where each piece of knowledge came from, and how the work was split between Pi, Claude and me.

Cluster numbers are computational labels, not cell types. UMAP is a 2-D projection, so a cell's position on it is not a measurement.

### API endpoints

| Endpoint | Returns |
|---|---|
| `GET /api/umap` | Coordinates, cluster labels, cell IDs and per-cell quality values |
| `GET /api/genes` | All gene names (used for autocomplete) |
| `GET /api/gene/{gene}` | One gene's expression across all cells |
| `GET /api/quality/{field}` | `n_genes`, `total_counts` or `pct_mito` for all cells |
| `GET /api/cell/{cell_id}` | One cell's cluster, UMAP position, quality values, stability and raw-count screen |
| `GET /api/clusters/{cluster}/genes` | Highest average expression and Wilcoxon-ranked markers |
| `GET /api/clusters/{cluster}/quality` | Quality values and summaries for one cluster |
| `GET /api/clusters/{cluster}/programmes` | Owner-supplied marker programme scores for one cluster |
| `GET /api/expression-summary/{cluster}` | Marker expression table for one cluster |
| `GET /api/dotplot?genes=A,B` | Mean expression and % positive cells per cluster (defaults to the owner's marker clues) |
| `GET /api/stability` | Clustering stability summary per cluster |
| `GET /api/stability/cells` | Per-cell clustering stability |
| `GET /api/reclustering` | Independent reclustering summary |
| `GET /api/doublets` | Raw-count doublet screen for clusters 6 and 7 |
| `GET /api/annotations` | Provisional annotations, with sources kept separate |
| `GET /api/report/{cluster}.csv` | CSV report for one cluster |

## Analysis scripts

The scripts in `analysis/` produced everything in `results/`. The app only reads those results; it does not re-run the analysis. To reproduce a result, run its script from the project folder, for example:

```bash
uv run python analysis/inspect_dataset.py
```

| Script | What it does |
|---|---|
| `inspect_dataset.py` | File structure, missing values, duplicates, and quality summaries per cluster |
| `rank_markers.py` | Wilcoxon cluster-versus-rest marker ranking |
| `analysis_report.py` | Evidence report for clusters 6 and 7 |
| `annotate_clusters.py` | Marker programme scores and provisional annotations |
| `validate_reclustering.py` | Independent reclustering from the raw counts at 4 resolutions |
| `stability_analysis.py` | 27-run clustering stability (3 seeds × 3 neighbour settings × 3 resolutions) |
| `doublet_investigation.py` | Raw-count doublet screening |

Some scripts append to `results/inspection.txt`, so re-running them adds duplicate sections to that file.

## Data

`data/pbmc3k.h5ad` contains:

- log-normalised expression in `X`;
- raw UMI counts in `layers["counts"]`;
- gene symbols in `var_names`;
- cluster labels in `obs["leiden"]`;
- per-cell quality fields in `obs["n_genes"]`, `obs["total_counts"]` and `obs["pct_mito"]`;
- UMAP coordinates in `obsm["X_umap"]`.

The facility computed the clusters and the map, and they were not changed.

## Where each piece of knowledge came from

- **Dataset-derived:** expression, quality values, clusters, UMAP, marker rankings, and the stability and screening metrics.
- **Transcript-derived:** the owner's goals, marker examples and decision criteria, from `ddls-week5-interview.md`.
- **Model pretrained knowledge:** biological interpretations not stated in the transcript. These are kept in a separate column in the annotation table and are not treated as evidence.

## Limits

- **Very few cells, and best readings rather than proof:** both answers rest on 10 and 13 cells. Platelet contamination rests on the owner's own `PPBP`/`PF4` clue. "Cycling" comes from the model's general biology knowledge, not from the transcript or the data, and cluster 7 could also be doublets. Cluster 7 is the least stable cluster (2 of its cells stay with it in only 21% of runs), and the doublet flag is a screening rule, not a validated classifier.
- **Read-level quality was never checked:** only the processed file was supplied, with no raw reads or run report.
- **The facility's filtering is unknown:** its cell filtering before delivery is undocumented, so which cells it removed, and with what thresholds, is unknown.
- **No mitochondrial cut-off:** none was supplied, and marker rankings do not prove identity.
