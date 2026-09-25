"""PBMC single-cell viewer with provenance-aware cluster exploration."""
from __future__ import annotations

import csv
import io
import os
import zipfile
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.responses import HTMLResponse, StreamingResponse

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / "ddls-week5-s1-junk-or-signal-dataset.zip"
DATA_PATH = ROOT / "data" / "pbmc3k.h5ad"
RESULTS = ROOT / "results"
MAX_MARKERS = 50
MARKERS = {
    "T-cell": ["CD3D", "CD3E", "TRBC1", "TRBC2"],
    "B-cell": ["MS4A1", "CD79A", "CD37"],
    "NK-cell": ["NKG7", "GNLY", "KLRD1"],
    "Monocyte": ["LYZ", "LST1", "S100A8", "S100A9"],
    "Dendritic-cell-like": ["FCER1A", "CST3"],
    "Platelet clue": ["PPBP", "PF4"],
}


def find_data_file() -> Path:
    configured = os.environ.get("H5AD_FILE")
    if configured:
        path = Path(configured)
        return path if path.is_absolute() else ROOT / path
    if DATA_PATH.exists():
        return DATA_PATH
    if ARCHIVE.exists():
        DATA_PATH.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(ARCHIVE) as archive:
            archive.extract("data/pbmc3k.h5ad", ROOT)
        return DATA_PATH
    matches = list(ROOT.glob("*.h5ad"))
    if matches:
        return matches[0]
    raise FileNotFoundError("No .h5ad file found")


def read_indexed(path: Path, id_column: str) -> pd.DataFrame | None:
    return pd.read_csv(path).set_index(id_column) if path.exists() else None


def rank_all_clusters(ad: Any) -> dict:
    """Wilcoxon cluster-vs-rest ranking for every cluster, computed once (same method as v1)."""
    ad.obs["cluster_for_ranking"] = ad.obs["leiden"].astype(str)
    sc.tl.rank_genes_groups(ad, groupby="cluster_for_ranking", method="wilcoxon", n_genes=MAX_MARKERS, key_added="markers")
    names, scores = ad.uns["markers"]["names"], ad.uns["markers"]["scores"]
    return {c: [{"gene": str(g), "score": float(v)} for g, v in zip(names[c], scores[c])] for c in names.dtype.names}


@asynccontextmanager
async def lifespan(app: FastAPI):
    app.state.ad = sc.read_h5ad(find_data_file())
    app.state.markers = rank_all_clusters(app.state.ad)
    stability_path = RESULTS / "stability_cluster_scores.csv"
    app.state.stability = pd.read_csv(stability_path).to_dict(orient="records") if stability_path.exists() else []
    app.state.cell_stability = read_indexed(RESULTS / "stability_cell_scores.csv", "cell_id")
    app.state.cell_screen = read_indexed(RESULTS / "doublet_raw_cell_metrics.csv", "index")
    summary_path = RESULTS / "doublet_cluster_summary.csv"
    app.state.doublet_summary = pd.read_csv(summary_path) if summary_path.exists() else None
    yield
    del app.state.ad


app = FastAPI(title="PBMC Cluster Viewer", lifespan=lifespan)


def dataset(request: Request) -> Any:
    return request.app.state.ad


def json_value(value: Any) -> Any:
    return value.item() if isinstance(value, np.generic) else value


def valid_cluster(ad: Any, cluster: str) -> np.ndarray:
    labels = ad.obs["leiden"].astype(str)
    if cluster not in set(labels):
        raise HTTPException(404, f"Unknown cluster: {cluster}")
    return labels.to_numpy() == cluster


def expression_vector(ad: Any, gene: str) -> np.ndarray:
    if gene not in ad.var_names:
        raise HTTPException(404, f"Unknown gene: {gene}")
    index = int(np.flatnonzero(ad.var_names == gene)[0])
    values = ad.X[:, index]
    return values.toarray().ravel() if hasattr(values, "toarray") else np.asarray(values).ravel()


def cluster_result(ad: Any, cluster: str, n_genes: int) -> dict:
    mask = valid_cluster(ad, cluster)
    subset = ad[mask]
    values = subset.X.toarray() if hasattr(subset.X, "toarray") else np.asarray(subset.X)
    means = np.asarray(values.mean(axis=0)).ravel()
    indices = np.argsort(means)[::-1][:n_genes]
    average = [{"gene": str(ad.var_names[i]), "value": float(means[i])} for i in indices]
    return {"cluster": cluster, "average_expression": average, "markers": app.state.markers[cluster][:n_genes]}


def row_dict(table: pd.DataFrame | None, cell_id: str) -> dict | None:
    if table is None or cell_id not in table.index:
        return None
    return {k: json_value(v) for k, v in table.loc[cell_id].items()}


def nearest_cell_ids(ids: Any, cell_id: str, n: int = 3) -> list[str]:
    def shared_prefix(other: str) -> int:
        length = 0
        for a, b in zip(other, cell_id):
            if a != b:
                break
            length += 1
        return length
    return sorted(map(str, ids), key=shared_prefix, reverse=True)[:n]


@app.get("/api/umap")
def umap(request: Request):
    ad = dataset(request)
    coords = np.asarray(ad.obsm["X_umap"])
    quality = {f: [json_value(v) for v in ad.obs[f].to_numpy()] for f in ("n_genes", "total_counts", "pct_mito") if f in ad.obs}
    return {"x": coords[:, 0].tolist(), "y": coords[:, 1].tolist(), "cluster": ad.obs["leiden"].astype(str).tolist(), "cell_id": [str(x) for x in ad.obs_names], **quality}


@app.get("/api/genes")
def gene_names(request: Request):
    return {"genes": [str(g) for g in dataset(request).var_names]}


@app.get("/api/dotplot")
def dotplot(request: Request, genes: str | None = Query(None)):
    """Mean log-normalised expression and % positive cells per cluster for a gene list (defaults to the owner's marker clues)."""
    ad = dataset(request)
    programme_by_gene = {gene: programme for programme, gene_list in MARKERS.items() for gene in gene_list}
    requested = [g.strip() for g in genes.split(",") if g.strip()] if genes else list(programme_by_gene)
    requested = list(dict.fromkeys(requested))
    present = [g for g in requested if g in ad.var_names]
    unknown = [g for g in requested if g not in ad.var_names and g not in programme_by_gene]
    if not present:
        raise HTTPException(404, "None of the requested genes are in the dataset")
    matrix = ad[:, present].X
    matrix = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
    labels = ad.obs["leiden"].astype(str).to_numpy()
    clusters = sorted(set(labels), key=int)
    mean = [matrix[labels == c].mean(axis=0).tolist() for c in clusters]
    pct = [((matrix[labels == c] > 0).mean(axis=0) * 100).tolist() for c in clusters]
    return {"genes": [{"gene": g, "programme": programme_by_gene.get(g, "Added")} for g in present], "clusters": clusters, "mean": mean, "pct": pct, "max_mean": float(max(max(row) for row in mean)), "unknown": unknown}


@app.get("/api/gene/{gene}")
def gene_expression(gene: str, request: Request):
    ad = dataset(request)
    return {"gene": gene, "values": expression_vector(ad, gene).tolist(), "cell_id": [str(x) for x in ad.obs_names]}


@app.get("/api/quality/{field}")
def quality_values(field: str, request: Request):
    ad = dataset(request)
    if field not in {"n_genes", "total_counts", "pct_mito"}:
        raise HTTPException(404, f"Unknown quality field: {field}")
    return {"field": field, "values": [json_value(x) for x in ad.obs[field].to_numpy()], "cell_id": [str(x) for x in ad.obs_names]}


@app.get("/api/clusters/{cluster}/genes")
def cluster_genes(cluster: str, request: Request, n_genes: int = Query(10, ge=1, le=50)):
    return cluster_result(dataset(request), cluster, n_genes)


@app.get("/api/stability")
def stability(request: Request):
    return {"method": "27 runs: 3 seeds × 3 neighbour settings × 3 Leiden resolutions", "clusters": request.app.state.stability}


@app.get("/api/stability/cells")
def stability_cells(request: Request):
    ad = dataset(request)
    table = request.app.state.cell_stability
    if table is None:
        raise HTTPException(404, "Per-cell stability results not found in results/")
    values = table["same_original_cocluster_frequency"].reindex(ad.obs_names)
    return {"field": "same_original_cocluster_frequency", "values": [None if pd.isna(v) else float(v) for v in values], "cell_id": [str(x) for x in ad.obs_names]}


@app.get("/api/cell/{cell_id}")
def cell_lookup(cell_id: str, request: Request):
    ad = dataset(request)
    if cell_id not in ad.obs_names:
        raise HTTPException(404, f"Unknown cell {cell_id}. Closest IDs: {', '.join(nearest_cell_ids(ad.obs_names, cell_id))}")
    x, y = (float(v) for v in np.asarray(ad.obsm["X_umap"])[ad.obs_names.get_loc(cell_id)])
    obs = ad.obs.loc[cell_id]
    return {
        "cell_id": cell_id,
        "cluster": str(obs["leiden"]),
        "umap": {"x": x, "y": y},
        "quality": {f: json_value(obs[f]) for f in ("n_genes", "total_counts", "pct_mito") if f in ad.obs},
        "stability": row_dict(request.app.state.cell_stability, cell_id),
        "raw_count_screen": row_dict(request.app.state.cell_screen, cell_id),
    }


@app.get("/api/reclustering")
def reclustering_summary():
    path = RESULTS / "reclustering_summary.csv"
    if not path.exists():
        raise HTTPException(404, "Reclustering summary not found in results/")
    table = pd.read_csv(path)
    return {"method": "Independent reclustering from the raw UMI counts: normalised to 10,000 counts per cell, log-transformed, 2,000 highly variable genes, scaled, 40-component PCA, 15-nearest-neighbour graph on 30 PCs, then Leiden at four resolutions. The facility's clusters were not overwritten.",
            "runs": [{"resolution": float(r.resolution), "n_clusters": int(r.n_reclustered), "mean_purity": float(r.mean_recluster_cluster_purity)} for r in table.itertuples()]}


@app.get("/api/doublets")
def doublet_summary(request: Request):
    table = request.app.state.doublet_summary
    if table is None:
        raise HTTPException(404, "Doublet summary not found in results/")
    rows = table[table["cluster"].astype(str).isin(["6", "7"])]
    others = table[~table["cluster"].astype(str).isin(["6", "7"])]["median_raw_counts"]
    return {"others_median_raw_counts": [float(others.min()), float(others.max())], "clusters": [{"cluster": str(r.cluster), "n_cells": int(r.n_cells), "median_raw_counts": float(r.median_raw_counts), "max_raw_counts": float(r.max_raw_counts), "mixed_cells": int(r.mixed_cells)} for r in rows.itertuples()],
            "method": "From results/doublet_cluster_summary.csv: raw UMI depth and co-expression of the owner-supplied marker programmes."}


@app.get("/api/annotations")
def annotations():
    """Provisional labels for all clusters, clusters 6 and 7 first. Dataset, transcript and model-knowledge evidence are kept in separate fields."""
    return {"rows": [
        {"cluster": "6", "label": "Platelet contamination", "basis": "Transcript clue, supported by model knowledge", "confidence": "Moderate–high", "dataset": "Top ranked markers include SDPR, GNG11, PF4, PPBP, NRGN, SPARC and GP9. All 13 cells are PPBP- and PF4-positive in raw counts.", "transcript": "The owner named PPBP and PF4 as platelet-contamination clues.", "model": "PF4, PPBP, GP9, ITGA2B, TUBB1, GNG11 and SPARC are classic platelet/megakaryocyte genes. Low gene counts fit small platelets or platelets stuck to other cells."},
        {"cluster": "7", "label": "Probably cycling cells", "basis": "Model knowledge", "confidence": "Low–moderate", "dataset": "Top ranked markers include ACTG1, CFL1, GAPDH, KIAA0101, STMN1, PCNA, TUBB and HMGB2. The T, B, NK and monocyte programmes are each positive in 8–10 of the 10 cells. High depth.", "transcript": "No coherent T, B, NK or monocyte programme; overlapping programmes count as “mixed”.", "model": "KIAA0101, STMN1, PCNA, TYMS, ZWINT, PTTG1 and HMGB2 are cell-cycle and proliferation genes, so these are likely dividing cells, possibly proliferating lymphocytes. High depth with mixed programmes could also mean doublets."},
        {"cluster": "0", "label": "Likely T cells (mainly CD4-like)", "basis": "Transcript clue + model knowledge", "confidence": "Moderate for T cells; low for the subtype", "dataset": "Top ranked markers are mostly ribosomal genes plus LDHB and CD3D. Largest cluster (1,197 cells) and the least stable one (median stability 0.67).", "transcript": "CD3D, one of the owner’s T-cell clues, is among the top 7 markers.", "model": "High ribosomal genes, LDHB and CD3D are typical of resting (naive or memory) CD4 T cells. Because it is large and less stable, it may also contain CD8 T cells."},
        {"cluster": "1", "label": "Classical (CD14+) monocytes", "basis": "Transcript clue + model knowledge", "confidence": "High for monocytes; moderate for the classical subtype", "dataset": "Top ranked markers: LYZ, S100A9, S100A8, TYROBP, CST3, FCN1.", "transcript": "LYZ and S100A8/A9 are the owner’s monocyte clues (S100A8/A9 for more inflammatory monocytes).", "model": "LYZ, S100A8/9 and FCN1 are typical of classical (CD14+) monocytes."},
        {"cluster": "2", "label": "Cytotoxic lymphocytes (NK and CD8 T-like)", "basis": "Transcript clue + model knowledge", "confidence": "Moderate", "dataset": "Top ranked markers: NKG7, GZMA, CST7, CTSW, CCL5, PRF1, GZMM, GZMB, FGFBP2.", "transcript": "NKG7 is one of the owner’s NK-cell clues.", "model": "Granzymes, perforin and NKG7 mark cytotoxic lymphocytes. CCL5 and GZMM lean towards CD8 T cells, GZMB and FGFBP2 towards NK cells, so this is likely a mix of both."},
        {"cluster": "3", "label": "B cells (likely naive)", "basis": "Transcript clue + model knowledge", "confidence": "High for B cells; moderate for naive", "dataset": "Top ranked markers: CD74, CD79A, HLA-DRA, CD79B, MS4A1, CD37, TCL1A, LINC00926. Perfectly stable (1.00).", "transcript": "CD79A, MS4A1 and CD37 are the owner’s B-cell clues.", "model": "CD79A/B and MS4A1 (CD20) define B cells; TCL1A and LINC00926 suggest naive B cells."},
        {"cluster": "4", "label": "Non-classical (CD16+) monocytes", "basis": "Transcript clue + model knowledge", "confidence": "Moderate–high", "dataset": "Top ranked markers: LST1, FCER1G, FCGR3A, COTL1, AIF1, IFITM2, IFITM3.", "transcript": "LST1 is one of the owner’s monocyte clues.", "model": "FCGR3A (CD16) together with LST1 and AIF1 marks non-classical (CD16+) monocytes."},
        {"cluster": "5", "label": "Dendritic cells (cDC2-like)", "basis": "Transcript clue + model knowledge", "confidence": "Moderate–high", "dataset": "Top ranked markers: HLA class II genes, CD74, CST3, FCER1A, CLEC10A. 36 cells, perfectly stable (1.00).", "transcript": "FCER1A and CST3 are the owner’s dendritic-cell-like clues.", "model": "FCER1A and CLEC10A with high HLA class II expression mark conventional dendritic cells (cDC2)."}
    ]}


@app.get("/api/clusters/{cluster}/quality")
def cluster_quality(cluster: str, request: Request):
    ad = dataset(request)
    mask = valid_cluster(ad, cluster)
    fields = [f for f in ("n_genes", "total_counts", "pct_mito") if f in ad.obs]
    return {"cluster": cluster, "n_cells": int(mask.sum()), "quality": {f: {"values": [json_value(x) for x in ad.obs.loc[mask, f].to_numpy()], "mean": float(np.mean(ad.obs.loc[mask, f])), "median": float(np.median(ad.obs.loc[mask, f]))} for f in fields}}


@app.get("/api/clusters/{cluster}/programmes")
def cluster_programmes(cluster: str, request: Request):
    ad = dataset(request)
    mask = valid_cluster(ad, cluster)
    result = []
    for name, genes in MARKERS.items():
        present = [g for g in genes if g in ad.var_names]
        x = ad[mask, present].X
        x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
        signal = x.sum(axis=1) if present else np.zeros(mask.sum())
        result.append({"programme": name, "genes": present, "positive_cells": int((signal > 0).sum()), "n_cells": int(mask.sum()), "mean_signal": float(signal.mean())})
    return {"cluster": cluster, "programmes": result}


@app.get("/api/expression-summary/{cluster}")
def expression_summary(cluster: str, request: Request):
    ad = dataset(request)
    mask = valid_cluster(ad, cluster)
    programme_by_gene = {gene: programme for programme, genes in MARKERS.items() for gene in genes}
    genes = [gene for gene in dict.fromkeys(sum(MARKERS.values(), [])) if gene in ad.var_names]
    if not genes:
        raise HTTPException(404, "No requested marker genes are present")
    matrix = ad[mask, genes].X
    matrix = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
    return {"cluster": cluster, "n_cells": int(mask.sum()), "genes": [{"gene": gene, "programme": programme_by_gene[gene], "mean": float(matrix[:, i].mean()), "median": float(np.median(matrix[:, i])), "pct_positive": float((matrix[:, i] > 0).mean() * 100)} for i, gene in enumerate(genes)]}


@app.get("/api/report/{cluster}.csv")
def report_csv(cluster: str, request: Request):
    result = cluster_result(dataset(request), cluster, 10)
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(["cluster", "section", "gene", "value"])
    for row in result["average_expression"]:
        writer.writerow([cluster, "highest_average_expression", row["gene"], row["value"]])
    for row in result["markers"]:
        writer.writerow([cluster, "ranked_marker", row["gene"], row["score"]])
    return StreamingResponse(iter([output.getvalue()]), media_type="text/csv", headers={"Content-Disposition": f'attachment; filename="cluster-{cluster}-report.csv"'})


HTML = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PBMC evidence viewer v2</title><script src="https://cdn.tailwindcss.com"></script><script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
.tab-btn{background:#fffefa}
.tab-btn[aria-selected="true"]{background:#3f3a34;color:#fff;border-color:#3f3a34}
.gene-link{cursor:pointer;text-decoration:underline dotted;text-underline-offset:3px}
.gene-link:hover{color:#2563eb}
.err{color:#b91c1c;font-size:.75rem;margin-top:.25rem}
#bigText[aria-pressed="true"]{background:#3f3a34;color:#fff}
.no-split li{break-inside:avoid}
</style></head>
<body class="bg-[#fcfaf5] text-[#3f3a34]"><main class="mx-auto max-w-7xl p-4 sm:p-6">

<header class="flex flex-wrap items-start justify-between gap-3"><div><p class="text-xs font-semibold uppercase tracking-widest text-slate-500">DDLS 2026 · Week 5 · Junk or signal?</p><h1 class="mt-1 text-3xl font-bold">A PBMC dataset, two tiny clusters, two decisions</h1><p class="mt-1 text-sm text-slate-600">Presented by <b>Sammi Baudin</b> · Data owner: <b>Dr. Ravi Menon</b>, immunologist · PBMC evidence viewer v2</p></div><button id="bigText" class="rounded border bg-white px-3 py-1.5 text-sm" aria-pressed="false">Larger text</button></header>

<section id="decisionStrip" class="mt-4 grid gap-3 sm:grid-cols-2" aria-label="Decisions for clusters 6 and 7"></section>

<nav class="mt-5 flex flex-wrap gap-2" role="tablist">
<button class="tab-btn rounded-full border px-4 py-1.5 text-sm font-medium" role="tab" data-tab="data" aria-selected="true">Data &amp; question</button>
<button class="tab-btn rounded-full border px-4 py-1.5 text-sm font-medium" role="tab" data-tab="explore" aria-selected="false">Explore</button>
<button class="tab-btn rounded-full border px-4 py-1.5 text-sm font-medium" role="tab" data-tab="compare" aria-selected="false">Clusters 6 &amp; 7</button>
<button class="tab-btn rounded-full border px-4 py-1.5 text-sm font-medium" role="tab" data-tab="validation" aria-selected="false">Methods &amp; validation</button>
<button class="tab-btn rounded-full border px-4 py-1.5 text-sm font-medium" role="tab" data-tab="answers" aria-selected="false">Answers, limits &amp; AI use</button>
</nav>

<!-- ================= EXPLORE ================= -->
<!-- ================= DATA & QUESTION ================= -->
<section data-panel="data" class="mt-4">
<div class="grid gap-4 lg:grid-cols-2">
<div class="rounded-xl bg-white p-5 shadow-sm">
<h2 class="text-xl font-bold">The data</h2>
<div class="mt-3 grid grid-cols-3 gap-2 text-center"><div class="rounded bg-slate-50 p-2"><div class="text-2xl font-bold tabular-nums">2,700</div><div class="text-xs text-slate-500">human PBMCs (cells)</div></div><div class="rounded bg-slate-50 p-2"><div class="text-2xl font-bold tabular-nums">13,714</div><div class="text-xs text-slate-500">genes</div></div><div class="rounded bg-slate-50 p-2"><div class="text-2xl font-bold tabular-nums">8</div><div class="text-xs text-slate-500">clusters, 0–7</div></div></div>
<ul class="mt-4 list-disc space-y-2 pl-5 text-sm">
<li><b>PBMC</b> = peripheral blood mononuclear cells: mainly lymphocytes and monocytes isolated from blood.</li>
<li>The blood was collected in a university hospital lab. The core facility isolated the cells, ran single-cell RNA sequencing and returned a <b>processed</b> file with the cells already grouped into clusters 0–7, a 2-D map (UMAP) and summary numbers.</li>
<li>The file, <code class="rounded bg-slate-100 px-1">pbmc3k.h5ad</code>, holds log-normalised expression, raw UMI counts, the facility’s cluster label (<code class="rounded bg-slate-100 px-1">leiden</code>), per-cell quality fields (<code class="rounded bg-slate-100 px-1">n_genes</code>, <code class="rounded bg-slate-100 px-1">total_counts</code>, <code class="rounded bg-slate-100 px-1">pct_mito</code>) and the UMAP coordinates. Every cell belongs to one cluster.</li>
<li>No methods sheet came with the file. The owner did not know the file’s layout, which processing had been done, or whether the reported “genes/cell” figures were means or medians.</li>
<li>The owner trusts the facility’s numbers most and the biological meaning of the numbered clusters least, because nobody has identified them.</li>
</ul>
<p class="mt-4 text-xs text-slate-500">Source: the data owner interview and the file description supplied with the dataset. Clusters and the map were computed by the facility and are not changed by this app.</p>
</div>

<div class="rounded-xl bg-white p-5 shadow-sm">
<h2 class="text-xl font-bold">The owner’s question</h2>
<p class="mt-1 text-sm text-slate-600">Dr. Ravi Menon, immunologist, is cleaning up the dataset before making a figure and has two decisions to make.</p>
<div class="mt-3 rounded-lg border-l-4 bg-slate-50 p-3" style="border-color:#d04482"><div class="text-xs font-semibold uppercase tracking-wide text-slate-500">Cluster 6 · 13 cells · ~350 genes/cell · 1.6% mito</div><div class="mt-1 font-semibold">“Is it safe to delete?”</div><p class="mt-1 text-sm text-slate-600">He was about to bin it as dead cells or empty droplets based only on the low gene count, without having looked at what the cells express.</p></div>
<div class="mt-3 rounded-lg border-l-4 bg-slate-50 p-3" style="border-color:#829b32"><div class="text-xs font-semibold uppercase tracking-wide text-slate-500">Cluster 7 · 10 cells · ~2,363 genes/cell</div><div class="mt-1 font-semibold">“Should the next sequencing run go to this cluster?”</div><p class="mt-1 text-sm text-slate-600">That spends facility time and budget. <b>Yes</b> only if the cluster has a clear, likely blood-cell identity (T, B, NK or monocyte); otherwise leave it out of the figure.</p></div>
<h3 class="mt-4 text-sm font-semibold">The owner’s own definitions</h3>
<dl class="mt-2 grid gap-2 text-sm sm:grid-cols-3"><div class="rounded bg-slate-50 p-2"><dt class="text-xs font-semibold uppercase text-slate-500">Coherent</dt><dd>Several genes from the same lineage across most cells, e.g. CD3D/E with TRBC for T cells.</dd></div><div class="rounded bg-slate-50 p-2"><dt class="text-xs font-semibold uppercase text-slate-500">Mixed</dt><dd>Strong programmes that don’t fit together in the same cells, e.g. T-cell and monocyte.</dd></div><div class="rounded bg-slate-50 p-2"><dt class="text-xs font-semibold uppercase text-slate-500">Generic / junk</dt><dd>Mostly mitochondrial, ribosomal, housekeeping or stress genes.</dd></div></dl>
<h3 class="mt-4 text-sm font-semibold">The owner’s marker clues</h3>
<p class="mt-1 text-sm">T cells <span class="gene-link" data-gene="CD3D">CD3D</span>, <span class="gene-link" data-gene="CD3E">CD3E</span>, TRBC1/2 · B cells <span class="gene-link" data-gene="MS4A1">MS4A1</span>, <span class="gene-link" data-gene="CD79A">CD79A</span>, <span class="gene-link" data-gene="CD37">CD37</span> · NK cells <span class="gene-link" data-gene="NKG7">NKG7</span>, <span class="gene-link" data-gene="GNLY">GNLY</span>, <span class="gene-link" data-gene="KLRD1">KLRD1</span> · Monocytes <span class="gene-link" data-gene="LYZ">LYZ</span>, <span class="gene-link" data-gene="LST1">LST1</span>, <span class="gene-link" data-gene="S100A8">S100A8</span>/<span class="gene-link" data-gene="S100A9">A9</span> · Dendritic-cell-like <span class="gene-link" data-gene="FCER1A">FCER1A</span>, <span class="gene-link" data-gene="CST3">CST3</span> · Platelet contamination <span class="gene-link" data-gene="PPBP">PPBP</span>, <span class="gene-link" data-gene="PF4">PF4</span></p>
<p class="mt-2 text-xs text-slate-500">These are clues, not proof. Click a gene to see it on the map. TRBC1/2 are not in this dataset.</p>
<button class="mt-4 rounded bg-blue-600 px-4 py-2 text-sm text-white" data-goto="explore">Explore the data →</button>
</div>
</div>
</section>

<section data-panel="explore" class="mt-4" hidden>
<div class="rounded-xl border border-[#d7e4da] bg-[#f1f7f2] p-3 text-sm"><b>How to read this:</b> each dot is one cell; UMAP places cells with similar measured expression near each other. Cluster colours are computational groups, not identities. “Highest expression” means abundant within a cluster; “ranked markers” means enriched versus all other cells. Neither alone proves identity.</div>

<div class="mt-4 grid grid-cols-1 gap-4 lg:grid-cols-[24rem_1fr]">
<div class="flex min-w-0 flex-col gap-4">
<div class="rounded-xl bg-[#fffefa] p-4 shadow-sm">
<h2 class="font-semibold">Select a cluster</h2>
<div id="clusterButtons" class="mt-2 grid grid-cols-4 gap-2"></div>
<button id="resetBtn" class="mt-3 w-full rounded border px-3 py-2 text-sm">Show all clusters</button>
<hr class="my-4">
<label for="colourMode" class="text-sm font-medium">Colour map by</label>
<div class="mt-2 flex gap-2"><select id="colourMode" class="min-w-0 flex-1 rounded border p-2"><option value="cluster">Cluster</option><option value="gene">Gene</option><option value="n_genes">Genes detected</option><option value="total_counts">Total counts</option><option value="pct_mito">Mitochondrial %</option><option value="stability">Clustering stability (per cell)</option></select><button id="colourBtn" class="rounded bg-blue-600 px-3 py-2 text-white">Apply</button></div>
<div id="geneControls" hidden><input id="gene" class="mt-2 w-full rounded border p-2" placeholder="e.g. PPBP" list="geneList" autocomplete="off"></div>
<p id="colourErr" class="err" hidden></p>
<p class="mt-2 text-xs text-slate-500">Selecting a cluster keeps the current colouring and greys out the other clusters.</p>
<label class="mt-3 flex items-center gap-2 text-sm"><input id="looseToggle" type="checkbox"> Mark cluster 7’s 2 loosely attached cells</label>
<hr class="my-4">
<label for="cellId" class="text-sm font-medium">Find a cell</label>
<div class="mt-2 flex gap-2"><input id="cellId" class="min-w-0 flex-1 rounded border p-2 text-sm" placeholder="e.g. CACAGAACCCTTGC-1" autocomplete="off"><button id="findBtn" class="rounded bg-blue-600 px-3 py-2 text-white">Find</button></div>
<p id="cellErr" class="err" hidden></p>
<p class="mt-1 text-xs text-slate-500">Or click any dot on the map.</p>
<div id="cellCard" class="mt-2 rounded bg-slate-50 p-2 text-xs leading-relaxed" hidden></div>
</div>

<div id="selectionCard" class="rounded-xl bg-white p-4 shadow-sm" hidden></div>

<section id="card" class="rounded-xl bg-white p-4 shadow-sm" hidden>
<div class="flex items-start justify-between gap-2"><div><h2 id="title" class="text-2xl font-bold"></h2><p id="quality" class="text-sm text-slate-600"></p></div><a id="download" class="shrink-0 rounded bg-slate-800 px-3 py-2 text-sm text-white" download>CSV</a></div>
<div id="decision" class="mt-3 rounded border-l-4 bg-slate-50 p-3 text-sm"></div>
<p id="confidence" class="mt-2 text-sm"></p>
<p class="mt-3 text-xs text-slate-500">Click a gene to colour the map by it. <b>+</b> adds it to the dot plot.</p>
<div class="mt-2 grid gap-3">
<div class="rounded border p-3"><h3 class="font-semibold">Highest average expression</h3><p class="text-xs text-slate-500">DATASET-DERIVED: within-cluster mean log-normalised expression.</p><p class="mt-2 rounded bg-blue-50 p-2 text-xs text-slate-700"><b>Question answered:</b> Which genes have the highest average expression within this cluster?</p><ul id="average" class="mt-2 divide-y text-sm"></ul></div>
<div class="rounded border p-3"><h3 class="font-semibold">Ranked marker genes</h3><p class="text-xs text-slate-500">DATASET-DERIVED: Wilcoxon ranking versus all other cells.</p><p class="mt-2 rounded bg-blue-50 p-2 text-xs text-slate-700"><b>Question answered:</b> Which genes distinguish this cluster from all other cells?</p><ul id="markers" class="mt-2 divide-y text-sm"></ul></div>
</div>
<div class="mt-4"><h3 class="font-semibold">Owner-supplied marker programmes</h3><p class="text-xs text-slate-500">TRANSCRIPT-DERIVED marker clues; presence is not proof of identity.</p><div id="programmes" class="mt-2 grid gap-2"></div></div>
<details class="mt-4 text-sm"><summary class="cursor-pointer font-semibold">Marker expression table</summary><div class="mt-2 overflow-x-auto"><table class="w-full text-left text-sm"><thead><tr class="border-b"><th class="p-2">Gene</th><th class="p-2">Programme</th><th class="p-2">Mean</th><th class="p-2">Median</th><th class="p-2">Cells positive</th></tr></thead><tbody id="expressionTable"></tbody></table></div></details>
<p id="cardErr" class="err" hidden></p>
</section>
</div>

<div class="min-w-0"><div class="rounded-xl bg-[#fffefa] p-2 shadow-sm lg:sticky lg:top-4">
<div id="plot" class="h-[60vh] min-h-[28rem] w-full lg:h-[calc(100vh-7rem)]"></div>
<p class="px-2 pb-1 text-xs text-slate-500">Hover a dot for its quality values · click a dot to look it up · use the lasso tool (top right) to summarise a group of cells.</p>
</div></div>
</div>

<section class="mt-4 rounded-xl bg-white p-4 shadow-sm">
<div class="flex flex-wrap items-end justify-between gap-3">
<div><h2 class="text-xl font-bold">Marker dot plot, all clusters</h2><p class="max-w-3xl text-xs text-slate-500">DATASET-DERIVED values for TRANSCRIPT-DERIVED marker clues. Dot colour: mean log-normalised expression on one shared scale for every cluster. Dot size: % of cells in the cluster expressing the gene (0–100%). The selected cluster’s row is highlighted. Click a dot to colour the map by that gene.</p></div>
<div class="w-full sm:w-auto"><div class="flex flex-wrap gap-2"><input id="dotGene" class="min-w-0 flex-1 rounded border p-2 text-sm sm:w-40 sm:flex-none" placeholder="Add a gene" list="geneList" autocomplete="off"><button id="dotAddBtn" class="rounded bg-blue-600 px-3 py-2 text-sm text-white">Add</button><button id="dotResetBtn" class="rounded border px-3 py-2 text-sm">Reset</button></div><p id="dotErr" class="err" hidden></p></div>
</div>
<div class="mt-2 overflow-x-auto"><div id="dotplot" class="h-[26rem] min-w-[46rem]"></div></div>
</section>

<div class="mt-4 rounded-xl bg-white p-4 shadow-sm">
<div class="flex flex-wrap items-center justify-between gap-2"><h2 class="text-xl font-bold">Provisional cluster annotations</h2></div>
<p class="mt-1 text-xs text-slate-500">Dataset observations, the owner’s transcript clues, and the model’s general biology knowledge are kept in separate columns; “Label based on” says which of them each label rests on. Clusters 6 and 7, which the owner asked about, come first; the others are listed for context. All labels are provisional.</p>
<div class="mt-3 overflow-x-auto"><table class="w-full min-w-[56rem] text-left text-sm"><thead><tr class="border-b"><th class="p-2">Cluster</th><th class="p-2">Suggested label</th><th class="p-2">Label based on</th><th class="p-2">Confidence</th><th class="p-2">Dataset-derived</th><th class="p-2">Transcript-derived</th><th class="p-2">Model pretrained knowledge</th></tr></thead><tbody id="annotationRows"></tbody></table></div>
</div>
</section>

<!-- ================= CLUSTERS 6 & 7 ================= -->
<section data-panel="compare" class="mt-4" hidden>
<div class="rounded-xl bg-white p-4 shadow-sm">
<h2 class="text-xl font-bold">Should either cluster be merged, split or binned?</h2>
<p class="mt-1 text-xs text-slate-500">From the stability runs, the independent reclustering and the raw-count per-cell metrics already in <code>results/</code>.</p>
<div class="mt-3 grid gap-4 md:grid-cols-2">
<div class="rounded-lg border-t-4 p-3" style="border-color:#d04482"><div class="text-xs font-semibold uppercase tracking-wide text-slate-500">Cluster 6</div><div class="text-lg font-bold" style="color:#b0306a">Platelet contamination: a real, separate group, not junk.</div><div class="mt-2"><div class="grid gap-1 border-t py-2 text-sm sm:grid-cols-[9rem_1fr] sm:gap-3"><div class="font-semibold">Merge?<span class="block text-xs font-medium uppercase tracking-wide text-slate-500">No</span></div><div>All 13 cells stay together at all 4 reclustering resolutions, with no other cells joining them. 12/13 stay together in 98% of stability runs; the exception is one cell, at 74%.</div></div><div class="grid gap-1 border-t py-2 text-sm sm:grid-cols-[9rem_1fr] sm:gap-3"><div class="font-semibold">Split?<span class="block text-xs font-medium uppercase tracking-wide text-slate-500">No</span></div><div>No run ever divides it. One outlier cell has high depth (8,875 raw counts) and some monocyte and B-cell marker signal. It is a single cell, not a subgroup.</div></div><div class="grid gap-1 border-t py-2 text-sm sm:grid-cols-[9rem_1fr] sm:gap-3"><div class="font-semibold">Dead or empty?<span class="block text-xs font-medium uppercase tracking-wide text-slate-500">Unlikely</span></div><div>Mito is 0.7–3.2%, not high. All 13 cells express both <span class="gene-link" data-gene="PPBP">PPBP</span> and <span class="gene-link" data-gene="PF4">PF4</span> in raw counts (31–116 per cell), the owner’s platelet-contamination clue. The cells have low depth, but they are not blank.</div></div></div></div>
<div class="rounded-lg border-t-4 p-3" style="border-color:#829b32"><div class="text-xs font-semibold uppercase tracking-wide text-slate-500">Cluster 7</div><div class="text-lg font-bold" style="color:#5f7422">Probably cycling cells: an 8-cell core plus 2 loosely attached cells.</div><div class="mt-2"><div class="grid gap-1 border-t py-2 text-sm sm:grid-cols-[9rem_1fr] sm:gap-3"><div class="font-semibold">Split?<span class="block text-xs font-medium uppercase tracking-wide text-slate-500">Not into two clusters</span></div><div>The same 2 cells split off at 3 of 4 reclustering resolutions. They join a group made almost entirely of cluster 0 cells, not a new group of their own.</div></div><div class="grid gap-1 border-t py-2 text-sm sm:grid-cols-[9rem_1fr] sm:gap-3"><div class="font-semibold">Stability<span class="block text-xs font-medium uppercase tracking-wide text-slate-500">8 stable, 2 not</span></div><div>The 8 core cells stay together in 76–79% of runs. The 2 loosely attached cells stay only 21% of the time; they are the source of the “min 0.21” in the stability table.</div></div><div class="grid gap-1 border-t py-2 text-sm sm:grid-cols-[9rem_1fr] sm:gap-3"><div class="font-semibold">The 8 core cells<span class="block text-xs font-medium uppercase tracking-wide text-slate-500">Not uniform</span></div><div>Some core cells are high in NK markers and others high in B-cell markers, but no run separates them. With 8 cells, any sub-split would be guesswork.</div></div></div><label class="mt-2 flex items-center gap-2 text-xs"><input type="checkbox" data-show-loose> Show the 2 loosely attached cells on the Explore map</label></div>
</div>
</div>

<div class="mt-4 rounded-xl bg-white p-4 shadow-sm">
<div class="flex flex-wrap items-center gap-2 text-sm"><span class="font-semibold">Compare cluster</span><select id="cmpA" class="rounded border p-1.5"></select><span>with cluster</span><select id="cmpB" class="rounded border p-1.5"></select></div>
<p class="mt-1 text-xs text-slate-500">Stored quality values, owner-supplied marker programmes and top ranked markers, side by side. Click a gene to colour the map by it.</p>
<div id="compareGrid" class="mt-4 grid gap-4 md:grid-cols-2"></div>
</div>
</section>

<!-- ================= VALIDATION ================= -->
<section data-panel="validation" class="mt-4" hidden>
<div class="rounded-xl bg-white p-4 shadow-sm">
<h2 class="text-xl font-bold">How the analysis was done</h2>
<p class="mt-1 text-xs text-slate-500">Why each step was done, and what it showed. Clusters and UMAP came from the facility and were not changed.</p>
<div class="mt-2"><div class="grid gap-1 border-t py-2 text-sm md:grid-cols-[13rem_1fr] md:gap-4"><div class="font-semibold">Inspect the file<span class="block text-xs font-normal text-slate-500">The owner didn’t know the layout</span></div><div>Rows are cells, columns are genes. 0 missing values and 0 duplicate cells. Raw counts are kept in a separate layer.</div></div><div class="grid gap-1 border-t py-2 text-sm md:grid-cols-[13rem_1fr] md:gap-4"><div class="font-semibold">Per-cell quality<span class="block text-xs font-normal text-slate-500">The “350 genes” is a summary figure</span></div><div>Cluster 6: the median is 350 genes, but the range runs from 212 to 2,455, with one cell at 8,931 counts. Mito 0.7–3.2%, which is not high.</div></div><div class="grid gap-1 border-t py-2 text-sm md:grid-cols-[13rem_1fr] md:gap-4"><div class="font-semibold">Marker ranking<span class="block text-xs font-normal text-slate-500">The owner asked for “top genes”; Wilcoxon, cluster vs rest</span></div><div>Cluster 6: <span class="gene-link" data-gene="PPBP">PPBP</span> and <span class="gene-link" data-gene="PF4">PF4</span> are in the top 4, the owner’s platelet clue. Cluster 7: <span class="gene-link" data-gene="ACTG1">ACTG1</span>, <span class="gene-link" data-gene="GAPDH">GAPDH</span>, <span class="gene-link" data-gene="STMN1">STMN1</span>, <span class="gene-link" data-gene="PCNA">PCNA</span> are not lineage markers.</div></div><div class="grid gap-1 border-t py-2 text-sm md:grid-cols-[13rem_1fr] md:gap-4"><div class="font-semibold">Owner’s marker programmes<span class="block text-xs font-normal text-slate-500">Checks the coherent vs mixed criterion, cell by cell</span></div><div>Cluster 7 positive cells: T 8/10 · B 10/10 · NK 9/10 · Mono 10/10, so the programmes overlap. Cluster 6: Mono 7/13; T, B and NK each ≤3/13.</div></div><div class="grid gap-1 border-t py-2 text-sm md:grid-cols-[13rem_1fr] md:gap-4"><div class="font-semibold">Where each piece of knowledge came from<span class="block text-xs font-normal text-slate-500">Unsourced biology was challenged</span></div><div>Every report is split into <i>dataset</i>, <i>transcript</i> and <i>model pretrained</i> knowledge.</div></div><div class="grid gap-1 border-t py-2 text-sm md:grid-cols-[13rem_1fr] md:gap-4"><div class="font-semibold">This evidence navigator<span class="block text-xs font-normal text-slate-500">So the owner can see the evidence directly</span></div><div>A FastAPI + Scanpy app: the map, gene lookup, marker lists, dot plot, and the validation results below.</div></div></div>
</div>

<div class="mt-4 rounded-xl bg-white p-4 shadow-sm">
<h2 class="text-xl font-bold">Independent reclustering from raw counts</h2>
<p id="reclusterMethod" class="mt-1 text-xs text-slate-500"></p>
<div class="mt-3 overflow-x-auto"><table class="w-full text-left text-sm"><thead><tr class="border-b"><th class="p-2">Leiden resolution</th><th class="p-2">Clusters found</th><th class="p-2">Mean cluster purity</th></tr></thead><tbody id="reclusterRows"></tbody></table></div>
<p class="mt-2 text-sm">Mean cluster purity was <b id="purityRange"></b>: the new clusters mostly reproduce the facility’s groups. Clusters 6 and 7 are small, so they are the most sensitive to these settings (see the Clusters 6 &amp; 7 tab).</p>
<p class="mt-1 text-xs text-slate-500">A validation sensitivity check, not proof that one clustering is biologically correct.</p>
</div>

<div class="mt-4 rounded-xl bg-white p-4 shadow-sm">
<h2 class="text-xl font-bold">Clustering stability &amp; doublet screening</h2>
<p id="stabilityMethod" class="mt-1 text-xs text-slate-500"></p>
<div class="mt-3 overflow-x-auto"><table class="w-full text-left text-sm"><thead><tr class="border-b"><th class="p-2">Cluster</th><th class="p-2">Cells</th><th class="p-2">Median stability</th><th class="p-2">Range</th></tr></thead><tbody id="stabilityRows"></tbody></table></div>
<div class="mt-3 grid gap-2 sm:grid-cols-2"><div class="rounded bg-amber-50 p-3 text-sm"><b>Cluster 6:</b> Computationally stable; this does not establish identity or deletion safety.</div><div class="rounded bg-blue-50 p-3 text-sm"><b>Cluster 7:</b> Less stable; 2 of its 10 cells stay with it in only 21% of runs.</div></div>
<p class="mt-3 text-xs text-slate-500">DATASET-DERIVED stability scores. Conclusions are computational observations, not cell-type annotations.</p>
<div class="mt-4 border-t pt-4"><h3 class="font-semibold">Doublet screening summary</h3><p class="mt-1 text-sm"><b>Why:</b> unusually high counts plus markers from two lineages can mean a doublet, two cells captured together that look like a new cell type.</p><p id="doubletMethod" class="mt-1 text-xs text-slate-500"></p><div id="doubletRows" class="mt-2 grid gap-2 sm:grid-cols-2"></div><p id="doubletCompare" class="mt-2 text-sm"></p><p class="mt-1 text-sm"><b>Conclusion:</b> cluster 7 looks like it could be made of doublets rather than one real cell type. This is a screening flag, not proof.</p><p class="mt-2 text-xs text-slate-500">Screening flags are not doublet diagnoses. They use raw counts and marker co-expression; dedicated validation such as Scrublet would be needed for stronger evidence.</p></div>
</div>


</section>

<!-- ================= ANSWERS, LIMITS & AI USE ================= -->
<section data-panel="answers" class="mt-4" hidden>
<div class="grid gap-4 md:grid-cols-2">
<div class="rounded-xl p-4" style="background:#fbeaf2"><div class="text-xs font-semibold uppercase tracking-wide text-slate-500">Cluster 6</div><div class="mt-1 text-2xl font-bold" style="color:#b0306a">Platelet contamination, not dead or empty cells.</div><p class="mt-2 text-sm">All 13 cells express PPBP and PF4, the owner’s own platelet-contamination clue, and mito is low. Deleting it as dead or empty junk would be wrong; removing it as platelet contamination is the owner’s call.</p><p class="mt-2 text-xs"><b>Confidence:</b> moderate–high that this is platelet contamination (the owner’s own PPBP/PF4 clue) rather than dead or empty material.</p></div>
<div class="rounded-xl p-4" style="background:#eef3e1"><div class="text-xs font-semibold uppercase tracking-wide text-slate-500">Cluster 7</div><div class="mt-1 text-2xl font-bold" style="color:#5f7422">No. Don’t prioritise the next run.</div><p class="mt-2 text-sm">The cluster most likely holds cycling (dividing) cells: its top markers include STMN1, PCNA, TYMS and PTTG1. That is not one of the T, B, NK or monocyte identities the owner asked for, and those programmes overlap in the same cells.</p><p class="mt-2 text-xs"><b>Confidence:</b> low–moderate for the recommendation and for the cycling reading, which comes from model knowledge; moderate for the 8 + 2 structure.</p></div>
</div>

<div class="mt-4 rounded-xl bg-white p-4 shadow-sm">
<h2 class="text-xl font-bold">Limits</h2>
<ul class="mt-2 list-disc space-y-1.5 pl-5 text-sm">
<li>Both answers rest on <b>10 and 13 cells</b>, and both identities are best readings, not proof. Platelet contamination rests on the owner’s own PPBP/PF4 clue. “Cycling” comes from the model’s general biology knowledge, not from the transcript or the data, and cluster 7 could also be doublets: it is the least stable cluster (some cells stay with it in only 21% of runs), and the doublet flag is a screening rule, not a validated classifier.</li>
<li>Read-level quality (base quality, mapping rate, saturation) was never checked: only the processed file was supplied, with no raw reads or run report.</li>
<li>The facility’s cell filtering before delivery is undocumented, so the cells it removed, and the thresholds it used, are unknown.</li>
<li>UMAP is a 2-D projection, not a measurement; a cell’s position on the map can mislead (e.g. <span class="gene-link" data-cell="CACAGAACCCTTGC-1">CACAGAACCCTTGC-1</span>, a stable cluster 3 cell drawn between clusters 2 and 0). No universal mitochondrial cut-off was supplied. Marker rankings do not prove identity.</li>
</ul>
<h3 class="mt-4 text-sm font-semibold">Where each piece of knowledge came from</h3>
<ul class="mt-1 list-disc space-y-1 pl-5 text-sm text-slate-700"><li><b>Dataset-derived:</b> expression, quality, clusters, UMAP, rankings, stability and screening metrics.</li><li><b>Transcript-derived:</b> the owner’s goals, marker examples and decision criteria.</li><li><b>Model pretrained knowledge:</b> kept in a separate column in the annotation table and not used as established evidence.</li></ul>
</div>

<div class="mt-4 rounded-xl bg-white p-4 shadow-sm">
<h2 class="text-xl font-bold">How the work was split</h2>
<div class="mt-3 grid gap-4 lg:grid-cols-[1fr_1fr_1.4fr]">
<div><h3 class="font-semibold">Pi <span class="text-xs font-normal text-slate-500">coding agent</span></h3><ul class="mt-1 list-disc space-y-1 pl-5 text-sm"><li>Drafted <code>AGENTS.md</code> and <code>spec.md</code> from my interview transcript</li><li>Set up the <code>uv</code> environment; wrote all 7 analysis scripts</li><li>Built the FastAPI/Scanpy navigator (version 1) and made the git commits</li><li>Proposed the plan, annotation, reclustering and stability designs; I approved them step by step</li></ul></div>
<div><h3 class="font-semibold">Claude <span class="text-xs font-normal text-slate-500">Claude Code</span></h3><ul class="mt-1 list-disc space-y-1 pl-5 text-sm"><li>Built the seminar slide deck from the existing transcripts and result files</li><li>Built this version 2 of the app from version 1: cell lookup, stability colouring, tabs, the all-cluster dot plot and these summary pages</li><li>Looked up existing results for the cluster 6 and 7 closer look</li><li>No new analysis</li></ul></div>
<div class="lg:border-l lg:pl-4"><h3 class="font-semibold">Me</h3><ul class="no-split mt-1 list-disc pl-5 text-sm sm:columns-2 sm:gap-8 [&>li]:mb-1"><li>Ran the owner interview and pinned down “coherent”, “mixed” and “junk”</li><li>Set the rules: reviewable scripts, commit often</li><li>Asked for a plan first and ran it in stages</li><li>Caught the 0.00-score ranking bug</li><li>Challenged unsourced biology, which led to the provenance split</li><li>Chose 3 seeds to keep compute low</li><li>Asked for the raw-count doublet check</li><li>Designed the app layout</li></ul></div>
</div>
</div>
</section>

<datalist id="geneList"></datalist>
</main>
<script>
let data
const colors=['#3b82c4','#d6534f','#45a85a','#8956b8','#e58a2b','#168fa3','#d04482','#829b32']
const cmap=Object.fromEntries(colors.map((x,i)=>[String(i),x]))
const $=id=>document.getElementById(id)
const esc=s=>String(s).replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]))
const fmt=n=>Math.round(n).toLocaleString()
const S={selected:null,values:null,label:'Cluster',hl:null,stab:null,loose:[],showLoose:false,dotGenes:null,dot:null,annotations:[],compared:false}

async function get(url){const r=await fetch(url);if(!r.ok){let m=r.statusText;try{m=(await r.json()).detail||m}catch(e){}throw Error(m)}return r.json()}
function showErr(id,msg){const el=$(id);el.textContent=msg||'';el.hidden=!msg}
function median(a){const s=[...a].sort((x,y)=>x-y),m=s.length>>1;return s.length?(s.length%2?s[m]:(s[m-1]+s[m])/2):NaN}

function decisions(c){
 if(c==='6')return ["Platelet contamination, not dead or empty cells.","All 13 cells express PPBP and PF4 in raw counts, the owner's platelet-contamination clue, and the mitochondrial share is low (0.7–3.2%). Deleting it as dead or empty junk would be wrong; removing it as platelet contamination is the owner's call.","Moderate–high that this is platelet contamination (the owner’s own PPBP/PF4 clue) rather than dead or empty material."]
 if(c==='7')return ["Do not prioritise for the next sequencing run: probably cycling cells.","Its top markers include STMN1, PCNA, TYMS and PTTG1, which suggests cycling (dividing) cells. That is not one of the T, B, NK or monocyte identities the owner asked for, and those programmes overlap in the same cells. The cluster is an 8-cell core plus 2 loosely attached cells that group with cluster 0 at 3 of 4 reclustering resolutions.","Low–moderate for the recommendation and for the cycling reading, which comes from model knowledge; moderate for the 8 + 2 structure."]
 return ["Exploratory view.","No owner decision is assigned to this cluster.","Not assessed."]}

/* ---------- tabs ---------- */
function showTab(t){
 document.querySelectorAll('[data-tab]').forEach(b=>b.setAttribute('aria-selected',String(b.dataset.tab===t)))
 document.querySelectorAll('[data-panel]').forEach(p=>p.hidden=p.dataset.panel!==t)
 if(t==='explore'){try{Plotly.Plots.resize('plot')}catch(e){}try{Plotly.Plots.resize('dotplot')}catch(e){}}
 if(t==='compare'&&!S.compared){S.compared=true;renderCompare()}}

/* ---------- decision strip ---------- */
function renderStrip(){
 $('decisionStrip').innerHTML=['6','7'].map(c=>{const d=decisions(c),n=data.cluster.filter(x=>x===c).length
  return `<button class="rounded-xl border-l-4 bg-white p-4 text-left shadow-sm hover:shadow-md" style="border-color:${cmap[c]}" data-strip="${c}"><div class="text-xs font-semibold uppercase tracking-wide text-slate-500">Cluster ${c} · ${n} cells</div><div class="mt-1 text-lg font-bold">${d[0]}</div><div class="mt-1 text-xs text-slate-600"><b>Confidence:</b> ${d[2]}</div><div class="mt-1 text-xs text-blue-700">Show evidence →</div></button>`}).join('')
 document.querySelectorAll('[data-strip]').forEach(b=>b.onclick=()=>{showTab('explore');select(b.dataset.strip);$('card').scrollIntoView({behavior:'smooth',block:'start'})})}

/* ---------- UMAP ---------- */
function plot(){
 const n=data.x.length,idx=[...Array(n).keys()],traces=[]
 const pick=ids=>({x:ids.map(i=>data.x[i]),y:ids.map(i=>data.y[i]),text:ids.map(i=>data.hover[i]),customdata:ids})
 const base={mode:'markers',type:'scattergl',hovertemplate:'%{text}<extra></extra>'}
 if(!S.values){traces.push({...base,...pick(idx),marker:{size:7,color:idx.map(i=>S.selected===null||data.cluster[i]===S.selected?cmap[data.cluster[i]]:'#cbd5e1')}})}
 else{const vals=S.values.filter(v=>v!==null&&v!==undefined);let cmin=Infinity,cmax=-Infinity;vals.forEach(v=>{if(v<cmin)cmin=v;if(v>cmax)cmax=v})
  const inSel=idx.filter(i=>S.selected===null||data.cluster[i]===S.selected),out=idx.filter(i=>S.selected!==null&&data.cluster[i]!==S.selected)
  if(out.length)traces.push({...base,...pick(out),marker:{size:6,color:'#e2e8f0'}})
  traces.push({...base,...pick(inSel),marker:{size:7,color:inSel.map(i=>S.values[i]),colorscale:'Viridis',cmin,cmax,showscale:true,colorbar:{title:{text:S.label}}}})}
 if(S.showLoose&&S.loose.length)traces.push({mode:'markers',type:'scatter',...pick(S.loose),text:S.loose.map(i=>data.hover[i]+'<br><b>loosely attached</b>: stays with cluster 7 in '+Math.round(100*S.stab[i])+'% of runs'),hovertemplate:'%{text}<extra></extra>',marker:{size:18,color:'rgba(0,0,0,0)',line:{color:'#dc2626',width:3}}})
 if(S.hl!==null)traces.push({mode:'markers',type:'scatter',...pick([S.hl]),hovertemplate:'%{text}<extra></extra>',marker:{size:22,color:'rgba(0,0,0,0)',line:{color:'#111827',width:3}}})
 Plotly.react('plot',traces,{margin:{l:45,r:15,t:15,b:45},xaxis:{title:{text:'UMAP 1'}},yaxis:{title:{text:'UMAP 2'}},dragmode:'pan',showlegend:false,uirevision:'map'},{displaylogo:false,responsive:true})}

async function colourBy(mode,gene){showErr('colourErr');try{
 if(mode==='cluster'){S.values=null;S.label='Cluster'}
 else if(mode==='gene'){if(!gene)throw Error('Type a gene name first.');const d=await get('/api/gene/'+encodeURIComponent(gene));S.values=d.values;S.label=gene}
 else if(mode==='stability'){const d=await get('/api/stability/cells');S.values=d.values;S.label='Stability'}
 else{const d=await get('/api/quality/'+mode);S.values=d.values;S.label=mode}
 plot()}catch(e){showErr('colourErr',e.message)}}
function apply(){colourBy($('colourMode').value,$('gene').value.trim())}
function clickGene(g){showTab('explore');$('colourMode').value='gene';$('geneControls').hidden=false;$('gene').value=g;colourBy('gene',g);if(window.innerWidth<1024)$('plot').scrollIntoView({behavior:'smooth'})}

async function toggleLoose(on){S.showLoose=on;showErr('colourErr');if(on&&!S.stab){try{const d=await get('/api/stability/cells');S.stab=d.values;S.loose=d.values.map((v,i)=>i).filter(i=>data.cluster[i]==='7'&&S.stab[i]!==null&&S.stab[i]<0.5)}catch(e){showErr('colourErr',e.message);$('looseToggle').checked=false;S.showLoose=false}}plot()}

/* ---------- single cell ---------- */
async function findCell(){const id=$('cellId').value.trim();showErr('cellErr');if(!id)return
 try{const d=await get('/api/cell/'+encodeURIComponent(id));S.hl=data.index.get(d.cell_id);plot()
  const q=d.quality,s=d.stability,r=d.raw_count_screen,card=$('cellCard')
  card.innerHTML=`<b>${esc(d.cell_id)}</b><br>Cluster ${d.cluster} · UMAP (${d.umap.x.toFixed(2)}, ${d.umap.y.toFixed(2)})<br>${q.n_genes} genes · ${fmt(q.total_counts)} counts · ${q.pct_mito.toFixed(1)}% mito`+(s?`<br>Stays with its own cluster in ${Math.round(100*s.same_original_cocluster_frequency)}% of the 27 stability runs`:'')+(r?`<br>Raw-count screen: ${r.possible_mixed_programmes?'flagged as possibly mixed':'not flagged'} (${r.n_strong_programmes} marker programme${r.n_strong_programmes==1?'':'s'})`:'')+`<br><button id="clearCell" class="mt-1 underline">Clear</button>`
  card.hidden=false;$('clearCell').onclick=()=>{S.hl=null;plot();card.hidden=true}}
 catch(e){showErr('cellErr',e.message)}}

/* ---------- lasso selection ---------- */
function showSelection(ids){const card=$('selectionCard');if(!ids.length){card.hidden=true;return}
 const counts={};ids.forEach(i=>counts[data.cluster[i]]=(counts[data.cluster[i]]||0)+1)
 card.innerHTML=`<div class="flex items-center justify-between"><h3 class="font-semibold">Selected cells: ${ids.length}</h3><button id="clearSel" class="text-xs underline">Clear</button></div><p class="mt-1 text-sm">Median ${median(ids.map(i=>data.n_genes[i]))} genes · ${fmt(median(ids.map(i=>data.total_counts[i])))} counts · ${median(ids.map(i=>data.pct_mito[i])).toFixed(1)}% mito</p><div class="mt-2 flex flex-wrap gap-1 text-xs">${Object.entries(counts).sort((a,b)=>b[1]-a[1]).map(([c,n])=>`<span class="rounded px-2 py-0.5 text-white" style="background:${cmap[c]}">Cluster ${c}: ${n}</span>`).join('')}</div><p class="mt-2 text-xs text-slate-500">Summary of stored values for the selected cells.</p>`
 card.hidden=false;$('clearSel').onclick=()=>{card.hidden=true;plot()}}

/* ---------- cluster card ---------- */
const geneItem=(g,val)=>`<li class="flex items-center justify-between gap-2 py-1"><span><span class="gene-link" data-gene="${esc(g)}">${esc(g)}</span> <button class="add-gene rounded border px-1 text-xs leading-4" data-gene="${esc(g)}" title="Add ${esc(g)} to the dot plot" aria-label="Add ${esc(g)} to the dot plot">+</button></span><span class="tabular-nums">${val}</span></li>`
async function select(c){S.selected=c;plot();if(S.dot)drawDot(S.dot);showErr('cardErr')
 try{const [g,q,p,e]=await Promise.all([get(`/api/clusters/${c}/genes?n_genes=10`),get(`/api/clusters/${c}/quality`),get(`/api/clusters/${c}/programmes`),get(`/api/expression-summary/${c}`)])
  $('title').innerHTML=`<span class="mr-2 inline-block h-3 w-3 rounded-full align-middle" style="background:${cmap[c]}"></span>Cluster ${c}`
  $('quality').textContent=`${q.n_cells} cells · median ${q.quality.n_genes.median} genes · ${fmt(q.quality.total_counts.median)} counts · ${q.quality.pct_mito.median.toFixed(1)}% mito`
  $('average').innerHTML=g.average_expression.map(x=>geneItem(x.gene,x.value.toFixed(3))).join('')
  $('markers').innerHTML=g.markers.map(x=>geneItem(x.gene,x.score.toFixed(3))).join('')
  $('programmes').innerHTML=p.programmes.map(x=>`<div class="rounded bg-slate-50 p-2 text-sm"><div class="flex justify-between"><b>${x.programme}</b><span class="tabular-nums">${x.positive_cells}/${x.n_cells} cells</span></div><div class="mt-1 h-1.5 rounded bg-slate-200"><div class="h-1.5 rounded bg-slate-600" style="width:${100*x.positive_cells/x.n_cells}%"></div></div><div class="mt-1 text-xs text-slate-500">${x.genes.map(gn=>`<span class="gene-link" data-gene="${esc(gn)}">${esc(gn)}</span>`).join(', ')}</div></div>`).join('')
  const d=decisions(c);$('decision').textContent=d[0]+' '+d[1];$('decision').style.borderColor=cmap[c];$('confidence').innerHTML='<b>Confidence:</b> '+d[2]
  $('download').href=`/api/report/${c}.csv`
  $('expressionTable').innerHTML=e.genes.map(x=>`<tr class="border-b"><td class="p-2 font-medium"><span class="gene-link" data-gene="${esc(x.gene)}">${esc(x.gene)}</span></td><td class="p-2">${x.programme}</td><td class="p-2">${x.mean.toFixed(3)}</td><td class="p-2">${x.median.toFixed(3)}</td><td class="p-2">${x.pct_positive.toFixed(1)}%</td></tr>`).join('')
  $('card').hidden=false}
 catch(err){$('card').hidden=false;showErr('cardErr',err.message)}}

/* ---------- dot plot ---------- */
async function loadDot(){showErr('dotErr');try{const q=S.dotGenes?'?genes='+encodeURIComponent(S.dotGenes.join(',')):'';const d=await get('/api/dotplot'+q);S.dot=d;if(!S.dotGenes)S.defaultGenes=d.genes.map(x=>x.gene);if(d.unknown.length)showErr('dotErr','Not in the dataset: '+d.unknown.join(', '));drawDot(d)}catch(e){showErr('dotErr',e.message)}}
function drawDot(d){
 const G=d.genes,C=d.clusters,xs=[],ys=[],sz=[],col=[],tx=[],gs=[]
 C.forEach((c,ci)=>G.forEach((g,gi)=>{const m=d.mean[ci][gi],p=d.pct[ci][gi];xs.push(gi);ys.push(ci);sz.push(3+24*p/100);col.push(m);gs.push(g.gene);tx.push(`Cluster ${c} · ${g.gene} (${g.programme})<br>mean ${m.toFixed(2)} · ${p.toFixed(0)}% of cells`)}))
 const shapes=[],ann=[];let start=0
 G.forEach((g,gi)=>{if(gi===G.length-1||G[gi+1].programme!==g.programme){const k=shapes.length;shapes.push({type:'rect',xref:'x',yref:'paper',x0:start-0.45,x1:gi+0.45,y0:1.01,y1:1.1,fillcolor:k%2?'#e2e8f0':'#f1f5f9',line:{width:0}});ann.push({xref:'x',yref:'paper',x:(start+gi)/2,y:1.055,text:g.programme,showarrow:false,font:{size:11}});start=gi+1}})
 if(S.selected!==null){const ci=C.indexOf(S.selected);if(ci>=0)shapes.push({type:'rect',xref:'paper',yref:'y',x0:0,x1:1,y0:ci-0.5,y1:ci+0.5,fillcolor:'rgba(250,204,21,0.28)',line:{width:0},layer:'below'})}
 Plotly.react('dotplot',[{x:xs,y:ys,mode:'markers',type:'scatter',text:tx,customdata:gs,hovertemplate:'%{text}<extra></extra>',marker:{size:sz,color:col,colorscale:'Reds',cmin:0,cmax:d.max_mean,showscale:true,colorbar:{title:{text:'Mean expr.'},thickness:12},line:{color:'#94a3b8',width:0.5}}}],
  {margin:{l:80,r:20,t:40,b:70},xaxis:{tickvals:G.map((_,i)=>i),ticktext:G.map(g=>g.gene),tickangle:-45,range:[-0.6,G.length-0.4],showgrid:false,zeroline:false},yaxis:{tickvals:C.map((_,i)=>i),ticktext:C.map(c=>'Cluster '+c),range:[C.length-0.5,-0.5],showgrid:false,zeroline:false},shapes,annotations:ann,showlegend:false},{displaylogo:false,responsive:true})}
function addDotGene(g){g=(g||'').trim();if(!g)return;const cur=S.dotGenes||(S.dot?S.dot.genes.map(x=>x.gene):[]);if(!cur.includes(g))S.dotGenes=[...cur,g];loadDot();$('dotGene').value=''}

/* ---------- compare ---------- */
async function renderCompare(){const a=$('cmpA').value,b=$('cmpB').value;$('compareGrid').innerHTML='<p class="text-sm text-slate-500">Loading…</p>'
 try{const cols=await Promise.all([a,b].map(async c=>{const [q,p,g]=await Promise.all([get(`/api/clusters/${c}/quality`),get(`/api/clusters/${c}/programmes`),get(`/api/clusters/${c}/genes?n_genes=5`)]);const d=decisions(c)
  return `<div class="rounded-xl border-t-4 bg-slate-50 p-4" style="border-color:${cmap[c]}"><h3 class="text-xl font-bold">Cluster ${c}</h3><p class="mt-2 text-sm"><b>${d[0]}</b> ${d[1]}</p><p class="mt-1 text-xs text-slate-600"><b>Confidence:</b> ${d[2]}</p>
  <dl class="mt-3 grid grid-cols-2 gap-x-4 gap-y-1 text-sm"><dt class="text-slate-500">Cells</dt><dd class="tabular-nums">${q.n_cells}</dd><dt class="text-slate-500">Median genes</dt><dd class="tabular-nums">${q.quality.n_genes.median}</dd><dt class="text-slate-500">Median counts</dt><dd class="tabular-nums">${fmt(q.quality.total_counts.median)}</dd><dt class="text-slate-500">Median mito %</dt><dd class="tabular-nums">${q.quality.pct_mito.median.toFixed(1)}</dd></dl>
  <h4 class="mt-3 text-sm font-semibold">Marker programmes (cells positive)</h4><div class="mt-1 space-y-1">${p.programmes.map(x=>`<div class="text-xs"><div class="flex justify-between"><span>${x.programme}</span><span class="tabular-nums">${x.positive_cells}/${x.n_cells}</span></div><div class="h-1.5 rounded bg-slate-200"><div class="h-1.5 rounded" style="width:${100*x.positive_cells/x.n_cells}%;background:${cmap[c]}"></div></div></div>`).join('')}</div>
  <h4 class="mt-3 text-sm font-semibold">Top ranked markers</h4><p class="mt-1 text-sm">${g.markers.map(x=>`<span class="gene-link" data-gene="${esc(x.gene)}">${esc(x.gene)}</span>`).join(', ')}</p></div>`}))
  $('compareGrid').innerHTML=cols.join('')}catch(e){$('compareGrid').textContent=e.message}}

/* ---------- validation ---------- */
function renderAnnotations(){$('annotationRows').innerHTML=S.annotations.map(x=>`<tr class="border-b align-top${['6','7'].includes(String(x.cluster))?' bg-amber-50/60':''}"><td class="p-2 font-semibold"><span class="mr-1 inline-block h-2.5 w-2.5 rounded-full" style="background:${cmap[x.cluster]}"></span>${x.cluster}</td><td class="p-2 font-medium">${esc(x.label)}</td><td class="p-2 text-xs">${esc(x.basis)}</td><td class="p-2">${esc(x.confidence)}</td><td class="p-2">${esc(x.dataset)}</td><td class="p-2">${esc(x.transcript)}</td><td class="p-2">${esc(x.model)}</td></tr>`).join('')}
async function loadAnnotations(){const a=await get('/api/annotations');S.annotations=a.rows;renderAnnotations()}
async function loadReclustering(){const r=await get('/api/reclustering');$('reclusterMethod').textContent=r.method;$('reclusterRows').innerHTML=r.runs.map(x=>`<tr class="border-b"><td class="p-2">${x.resolution}</td><td class="p-2">${x.n_clusters}</td><td class="p-2 tabular-nums">${x.mean_purity.toFixed(2)}</td></tr>`).join('');const p=r.runs.map(x=>x.mean_purity);$('purityRange').textContent=Math.min(...p).toFixed(2)+'–'+Math.max(...p).toFixed(2)}
async function loadDoublets(){const d=await get('/api/doublets');$('doubletMethod').textContent=d.method;const c7=d.clusters.find(x=>x.cluster==='7');if(c7)$('doubletCompare').innerHTML=`Cluster 7 has a median of <b>${fmt(c7.median_raw_counts)}</b> raw counts, against ${fmt(d.others_median_raw_counts[0])}–${fmt(d.others_median_raw_counts[1])} in clusters 0–5.`;$('doubletRows').innerHTML=d.clusters.map(x=>`<div class="rounded bg-slate-50 p-3 text-sm"><b>Cluster ${x.cluster}</b> · ${x.n_cells} cells<br>Median raw counts ${fmt(x.median_raw_counts)} (max ${fmt(x.max_raw_counts)})<br>Cells with more than one marker programme: ${x.mixed_cells}/${x.n_cells}</div>`).join('')}
async function loadStability(){const s=await get('/api/stability');$('stabilityMethod').textContent=s.method;$('stabilityRows').innerHTML=s.clusters.filter(x=>['6','7'].includes(String(x.cluster))).map(x=>`<tr class="border-b"><td class="p-2">${x.cluster}</td><td class="p-2">${x.n_cells}</td><td class="p-2">${Number(x.median_same_original_cocluster_frequency).toFixed(3)}</td><td class="p-2">${Number(x.min_same_original_cocluster_frequency).toFixed(3)}–${Number(x.max_same_original_cocluster_frequency).toFixed(3)}</td></tr>`).join('')}

/* ---------- start ---------- */
async function init(){data=await get('/api/umap')
 data.hover=data.x.map((_,i)=>`cell ${data.cell_id[i]}<br>cluster ${data.cluster[i]}<br>${data.n_genes[i]} genes · ${fmt(data.total_counts[i])} counts · ${data.pct_mito[i].toFixed(1)}% mito`)
 data.index=new Map(data.cell_id.map((id,i)=>[id,i]))
 const clusters=[...new Set(data.cluster)].sort((a,b)=>+a-+b)
 clusters.forEach(c=>{const b=document.createElement('button');b.textContent=c;b.title='Cluster '+c;b.setAttribute('aria-label','Cluster '+c);b.className='rounded py-2 text-sm font-semibold';b.style.background=cmap[c];b.style.color='white';b.onclick=()=>select(c);$('clusterButtons').appendChild(b)})
 ;['cmpA','cmpB'].forEach((id,k)=>{$(id).innerHTML=clusters.map(c=>`<option value="${c}">${c}</option>`).join('');$(id).value=k?'7':'6';$(id).onchange=renderCompare})
 renderStrip();plot()
 $('plot').on('plotly_click',ev=>{const pt=ev.points[0];if(pt&&pt.customdata!==undefined&&typeof pt.customdata==='number'){$('cellId').value=data.cell_id[pt.customdata];findCell()}})
 $('plot').on('plotly_selected',ev=>showSelection(ev&&ev.points?[...new Set(ev.points.map(p=>p.customdata).filter(v=>typeof v==='number'))]:[]))
 $('plot').on('plotly_deselect',()=>showSelection([]))
 await loadDot()
 $('dotplot').on('plotly_click',ev=>{const pt=ev.points[0];if(pt&&pt.customdata)clickGene(pt.customdata)})
 get('/api/genes').then(g=>{$('geneList').innerHTML=g.genes.map(x=>`<option value="${esc(x)}">`).join('')}).catch(e=>console.warn('Gene list unavailable',e))}

document.addEventListener('click',ev=>{const add=ev.target.closest('.add-gene');if(add){addDotGene(add.dataset.gene);return}const link=ev.target.closest('.gene-link');if(link&&link.dataset.cell){showTab('explore');$('cellId').value=link.dataset.cell;findCell();return}if(link)clickGene(link.dataset.gene)})
document.querySelectorAll('[data-tab]').forEach(b=>b.onclick=()=>showTab(b.dataset.tab));document.querySelectorAll('[data-goto]').forEach(b=>b.onclick=()=>showTab(b.dataset.goto))
$('colourMode').onchange=()=>{$('geneControls').hidden=$('colourMode').value!=='gene';showErr('colourErr')}
$('colourBtn').onclick=apply;$('gene').onkeydown=e=>{if(e.key==='Enter')apply()}
$('findBtn').onclick=findCell;$('cellId').onkeydown=e=>{if(e.key==='Enter')findCell()}
$('looseToggle').onchange=e=>{document.querySelectorAll('[data-show-loose]').forEach(b=>b.checked=e.target.checked);toggleLoose(e.target.checked)}
document.querySelectorAll('[data-show-loose]').forEach(b=>b.onchange=e=>{$('looseToggle').checked=e.target.checked;toggleLoose(e.target.checked)})
$('resetBtn').onclick=()=>{S.selected=null;plot();if(S.dot)drawDot(S.dot);$('card').hidden=true}
$('dotAddBtn').onclick=()=>addDotGene($('dotGene').value);$('dotGene').onkeydown=e=>{if(e.key==='Enter')addDotGene($('dotGene').value)}
$('dotResetBtn').onclick=()=>{S.dotGenes=null;loadDot()}
function setBig(on){document.documentElement.style.fontSize=on?'19px':'';$('bigText').setAttribute('aria-pressed',String(on));$('bigText').textContent=on?'Normal text':'Larger text';try{localStorage.setItem('bigText',on?'1':'')}catch(e){}setTimeout(()=>{try{Plotly.Plots.resize('plot')}catch(e){}try{Plotly.Plots.resize('dotplot')}catch(e){}},50)}
$('bigText').onclick=()=>setBig($('bigText').getAttribute('aria-pressed')!=='true')
try{if(localStorage.getItem('bigText'))setBig(true)}catch(e){}
loadStability().catch(e=>console.warn('Stability summary unavailable',e))
loadDoublets().catch(e=>console.warn('Doublet summary unavailable',e))
loadReclustering().catch(e=>{$('reclusterMethod').textContent=e.message})
loadAnnotations().catch(e=>console.warn('Annotations unavailable',e))
init().catch(e=>{document.querySelector('main').insertAdjacentHTML('afterbegin',`<p class="mb-4 rounded bg-red-50 p-3 text-sm text-red-700">Could not load the dataset: ${esc(e.message)}</p>`)})
</script></body></html>'''


@app.get("/", response_class=HTMLResponse)
def index():
    return HTML
