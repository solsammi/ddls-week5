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
MARKERS = {
    "T-cell": ["CD3D", "CD3E", "TRBC1", "TRBC2"],
    "B-cell": ["MS4A1", "CD79A", "CD37"],
    "NK-cell": ["NKG7", "GNLY", "KLRD1"],
    "Monocyte": ["LYZ", "LST1", "S100A8", "S100A9"],
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


@asynccontextmanager
async def lifespan(app: FastAPI):
    app.state.ad = sc.read_h5ad(find_data_file())
    app.state.marker_cache = {}
    stability_path = ROOT / "results" / "stability_cluster_scores.csv"
    app.state.stability = pd.read_csv(stability_path).to_dict(orient="records") if stability_path.exists() else []
    yield
    del app.state.ad


app = FastAPI(title="PBMC Cluster Viewer", lifespan=lifespan)


def dataset(request: Request) -> Any:
    return request.app.state.ad


def json_value(value: Any) -> Any:
    return value.item() if isinstance(value, (np.integer, np.floating)) else value


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
    key = (cluster, n_genes)
    if key not in app.state.marker_cache:
        ranking = ad.copy()
        ranking.obs["cluster_for_ranking"] = ranking.obs["leiden"].astype(str)
        sc.tl.rank_genes_groups(ranking, groupby="cluster_for_ranking", method="wilcoxon", n_genes=n_genes, key_added="markers")
        names = ranking.uns["markers"]["names"][cluster]
        scores = ranking.uns["markers"]["scores"][cluster]
        app.state.marker_cache[key] = [{"gene": str(g), "score": float(s)} for g, s in zip(names, scores)]
    return {"cluster": cluster, "average_expression": average, "markers": app.state.marker_cache[key]}


@app.get("/api/umap")
def umap(request: Request):
    ad = dataset(request)
    coords = np.asarray(ad.obsm["X_umap"])
    return {"x": coords[:, 0].tolist(), "y": coords[:, 1].tolist(), "cluster": ad.obs["leiden"].astype(str).tolist(), "cell_id": [str(x) for x in ad.obs_names]}


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


@app.get("/api/annotations")
def annotations():
    return {"rows": [
        {"cluster": "0", "dataset": "Mostly ribosomal/housekeeping genes; CD3D is among the top markers.", "transcript": "T-cell clue is present, but no coherent programme was established.", "model": "Possible T-cell-rich or low-specificity lymphocyte group.", "label": "Unresolved lymphocyte-like / low-specificity", "confidence": "Low"},
        {"cluster": "1", "dataset": "LYZ, S100A8, S100A9, TYROBP, CST3, FCN1.", "transcript": "LYZ and S100A8/A9 support the owner’s monocyte clues.", "model": "Inflammatory/classical monocyte-like.", "label": "Inflammatory monocyte-like", "confidence": "Moderate-high"},
        {"cluster": "2", "dataset": "NKG7, GZMA, CST7, CTSW, CCL5, PRF1, GZMB, FGFBP2.", "transcript": "NKG7 supports the owner’s NK-cell clue.", "model": "Cytotoxic NK-like or cytotoxic lymphocyte-like.", "label": "Cytotoxic NK-like", "confidence": "Moderate"},
        {"cluster": "3", "dataset": "CD74, CD79A, HLA-DRA, CD79B, MS4A1, CD37, TCL1A.", "transcript": "CD79A, MS4A1, and CD37 support the owner’s B-cell clues.", "model": "B-cell-like; possible naïve/transitional subtype.", "label": "B-cell-like", "confidence": "High for broad identity"},
        {"cluster": "4", "dataset": "LST1, FCER1G, FCGR3A, AIF1, CTSS, SERPINA1.", "transcript": "LST1 supports the owner’s monocyte clue.", "model": "FCGR3A-associated/non-classical monocyte-like.", "label": "FCGR3A-associated monocyte-like", "confidence": "Moderate-high"},
        {"cluster": "5", "dataset": "HLA-DPA1, HLA-DPB1, HLA-DRA, CD74, FCER1A, CLEC10A, CST3, LYZ.", "transcript": "FCER1A and CST3 support the owner’s dendritic-cell-like clues.", "model": "Dendritic-cell-like antigen-presenting population.", "label": "Dendritic-cell-like", "confidence": "Moderate"},
        {"cluster": "6", "dataset": "PF4, PPBP, GP9, ITGA2B, TUBB1, GNG11, SPARC.", "transcript": "PPBP and PF4 are explicitly platelet-contamination clues.", "model": "Platelet/megakaryocyte-associated interpretation.", "label": "Platelet-contamination signal / unresolved", "confidence": "Moderate"},
        {"cluster": "7", "dataset": "KIAA0101, STMN1, PCNA, TYMS, ZWINT, PTTG1; broad lineage-marker overlap.", "transcript": "No coherent T-, B-, NK-, or monocyte programme established.", "model": "Possible cycling/proliferative lymphocyte-like population.", "label": "Cycling/proliferative signal / unresolved", "confidence": "Low-moderate"}
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


HTML = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>PBMC evidence viewer</title><script src="https://cdn.tailwindcss.com"></script><script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script></head><body class="bg-slate-50 text-slate-900"><main class="mx-auto max-w-7xl p-4 sm:p-6"><header><h1 class="text-3xl font-bold">PBMC evidence viewer</h1><p class="mt-1 text-slate-600">Explore computational clusters without treating them as cell-type labels.</p></header><section class="mt-4 rounded-xl border border-blue-100 bg-blue-50 p-4 text-sm"><b>How to read this:</b> each dot is one cell; UMAP places cells with similar measured expression near each other. Cluster colours are computational groups, not identities. “Highest expression” means abundant within a cluster; “ranked markers” means enriched versus all other cells. Neither alone proves identity.</section><div class="mt-5 grid grid-cols-1 gap-4 lg:grid-cols-[18rem_1fr]"><aside class="rounded-xl bg-white p-4 shadow-sm"><h2 class="font-semibold">Select a cluster</h2><div id="clusterButtons" class="mt-2 grid grid-cols-2 gap-2"></div><button id="resetBtn" class="mt-3 w-full rounded border px-3 py-2 text-sm">Show all clusters</button><hr class="my-4"><label class="text-sm font-medium">Colour map by</label><div class="mt-2 flex gap-2"><select id="colourMode" class="min-w-0 flex-1 rounded border p-2"><option value="cluster">Cluster</option><option value="gene">Gene</option><option value="n_genes">Genes detected</option><option value="total_counts">Total counts</option><option value="pct_mito">Mitochondrial %</option></select><button id="colourBtn" class="rounded bg-blue-600 px-3 py-2 text-white">Apply</button></div><div id="geneControls" class="mt-2 hidden"><input id="gene" class="w-full rounded border p-2" placeholder="e.g. CD3D" list="genes"><datalist id="genes"></datalist></div><p class="mt-4 text-xs text-slate-500">Markers and decisions are shown with their evidence source and confidence.</p></aside><div class="min-w-0"><div class="rounded-xl bg-white p-2 shadow-sm"><div id="plot" class="h-[60vh] min-h-[28rem] w-full"></div></div></div><section id="card" class="col-span-1 mt-0 hidden w-full rounded-xl bg-white p-4 shadow-sm lg:col-span-2"><div class="flex flex-wrap items-start justify-between gap-2"><div><h2 id="title" class="text-2xl font-bold"></h2><p id="quality" class="text-sm text-slate-600"></p></div><a id="download" class="rounded bg-slate-800 px-3 py-2 text-sm text-white" download>Download CSV</a></div><div id="decision" class="mt-3 rounded border-l-4 p-3 text-sm"></div><p id="confidence" class="mt-2 text-sm"></p><div class="mt-5 grid grid-cols-1 gap-4 xl:grid-cols-2"><div class="rounded border p-3"><h3 class="font-semibold">Highest average expression</h3><p class="text-xs text-slate-500">DATASET-DERIVED: within-cluster mean log-normalised expression.</p><p class="mt-2 rounded bg-blue-50 p-2 text-xs text-slate-700"><b>Question answered:</b> Which genes have the highest average expression within this cluster?</p><ul id="average" class="mt-2 divide-y text-sm"></ul></div><div class="rounded border p-3"><h3 class="font-semibold">Ranked marker genes</h3><p class="text-xs text-slate-500">DATASET-DERIVED: Wilcoxon ranking versus all other cells.</p><p class="mt-2 rounded bg-blue-50 p-2 text-xs text-slate-700"><b>Question answered:</b> Which genes distinguish this cluster from all other cells?</p><ul id="markers" class="mt-2 divide-y text-sm"></ul></div></div><div class="mt-4"><h3 class="font-semibold">Owner-supplied marker programmes</h3><p class="text-xs text-slate-500">TRANSCRIPT-DERIVED marker clues; presence is not proof of identity.</p><div id="programmes" class="mt-2 grid gap-2 sm:grid-cols-2"></div></div><div class="mt-4 rounded border p-3"><h3 class="font-semibold">Marker expression summary</h3><p class="text-xs text-slate-500">DATASET-DERIVED: dot colour shows mean log-normalised expression; dot size shows the percentage of cells expressing the gene. This is stable across cluster sizes.</p><div id="dotplot" class="h-[30rem] w-full"></div><div class="mt-3 overflow-x-auto"><table class="w-full text-left text-sm"><thead><tr class="border-b"><th class="p-2">Gene</th><th class="p-2">Programme</th><th class="p-2">Mean</th><th class="p-2">Median</th><th class="p-2">Cells positive</th></tr></thead><tbody id="expressionTable"></tbody></table></div></div><details class="mt-4 text-sm"><summary class="cursor-pointer font-semibold">Limitations and provenance</summary><div class="mt-2 space-y-2 text-slate-600"><p><b>Dataset-derived:</b> expression, quality, clusters, UMAP, and rankings.</p><p><b>Transcript-derived:</b> owner goals, supplied marker examples, and decision criteria.</p><p><b>Model knowledge:</b> biological annotations not explicitly in the transcript are not used as established evidence.</p><p>Clusters 6 and 7 contain only 13 and 10 cells. UMAP is a projection, not a measurement. There is no universal mitochondrial cutoff supplied here. Marker rankings do not prove identity.</p></div></details></section><section id="stabilityCard" class="col-span-1 mt-4 w-full rounded-xl bg-white p-4 shadow-sm lg:col-span-2"><h2 class="text-xl font-bold">Clustering stability</h2><p id="stabilityMethod" class="mt-1 text-xs text-slate-500"></p><div class="mt-3 overflow-x-auto"><table class="w-full text-left text-sm"><thead><tr class="border-b"><th class="p-2">Cluster</th><th class="p-2">Cells</th><th class="p-2">Median stability</th><th class="p-2">Range</th></tr></thead><tbody id="stabilityRows"></tbody></table></div><div class="mt-3 grid gap-2 sm:grid-cols-2"><div class="rounded bg-amber-50 p-3 text-sm"><b>Cluster 6:</b> Computationally stable; this does not establish identity or deletion safety.</div><div class="rounded bg-blue-50 p-3 text-sm"><b>Cluster 7:</b> Less stable; some cells are sensitive to settings, so treat its boundary cautiously.</div></div><p class="mt-3 text-xs text-slate-500">DATASET-DERIVED stability scores. Conclusions are computational observations, not cell-type annotations.</p></section><section id="annotationCard" class="col-span-1 mt-4 w-full rounded-xl bg-white p-4 shadow-sm lg:col-span-2"><h2 class="text-xl font-bold">Provisional cluster annotations</h2><p class="mt-1 text-xs text-slate-500">Dataset observations, transcript-supported interpretation, and model biological knowledge are kept in separate columns.</p><div class="mt-3 overflow-x-auto"><table class="w-full min-w-[70rem] text-left text-sm"><thead><tr class="border-b"><th class="p-2">Cluster</th><th class="p-2">Suggested label</th><th class="p-2">Confidence</th><th class="p-2">Dataset-derived</th><th class="p-2">Transcript-derived</th><th class="p-2">Model pretrained knowledge</th></tr></thead><tbody id="annotationRows"></tbody></table></div></section></div></main><script>
let data;const colors=['#2563eb','#dc2626','#16a34a','#9333ea','#ea580c','#0891b2','#db2777','#65a30d'];const cmap=Object.fromEntries(colors.map((x,i)=>[String(i),x]));
const $=id=>document.getElementById(id);async function get(url){const r=await fetch(url);if(!r.ok)throw Error((await r.json()).detail||r.statusText);return r.json()}
function plot(selected=null,values=null,label='Cluster'){const c=values||data.cluster.map(x=>selected===null?cmap[x]:x===selected?cmap[x]:'#cbd5e1');Plotly.react('plot',[{x:data.x,y:data.y,mode:'markers',type:'scattergl',text:data.cell_id,customdata:data.cluster,marker:{size:7,color:c,colorscale:'Viridis',showscale:!!values,colorbar:{title:label}},hovertemplate:'cell %{text}<br>cluster %{customdata}<extra></extra>'}],{margin:{l:45,r:15,t:15,b:45},xaxis:{title:'UMAP 1'},yaxis:{title:'UMAP 2'},dragmode:'pan'});}
function decisions(c){if(c==='6')return ['Do not delete yet.','Dataset quality is heterogeneous; review the actual expression before removal.','High confidence that size/median gene count alone is insufficient; low confidence for identity without further annotation.'];if(c==='7')return ['Do not prioritise yet.','The supplied lineage programmes overlap; no single requested identity is established.','Moderate confidence for the quality observation; low-to-moderate confidence for biological identity.'];return ['Exploratory view.','No owner decision is assigned to this cluster.','Not assessed.']}
async function select(c){plot(c);const [g,q,p,e]=await Promise.all([get(`/api/clusters/${c}/genes?n_genes=10`),get(`/api/clusters/${c}/quality`),get(`/api/clusters/${c}/programmes`),get(`/api/expression-summary/${c}`)]);$('card').classList.remove('hidden');$('title').textContent='Cluster '+c;$('quality').textContent=`${q.n_cells} cells · n_genes median ${q.quality.n_genes.median} · total counts median ${q.quality.total_counts.median} · pct_mito median ${q.quality.pct_mito.median.toFixed(2)}%`;$('average').innerHTML=g.average_expression.map(x=>`<li class="flex justify-between py-1"><span>${x.gene}</span><span>${x.value.toFixed(3)}</span></li>`).join('');$('markers').innerHTML=g.markers.map(x=>`<li class="flex justify-between py-1"><span>${x.gene}</span><span>${x.score.toFixed(3)}</span></li>`).join('');$('programmes').innerHTML=p.programmes.map(x=>`<div class="rounded bg-slate-50 p-2 text-sm"><b>${x.programme}</b>: ${x.positive_cells}/${x.n_cells} cells<br><span class="text-xs text-slate-500">${x.genes.join(', ')}</span></div>`).join('');const d=decisions(c);$('decision').textContent=d[0]+' '+d[1];$('decision').className='mt-3 rounded border-l-4 p-3 text-sm '+(c==='6'?'border-amber-500 bg-amber-50':'border-blue-500 bg-blue-50');$('confidence').textContent='Confidence: '+d[2];$('download').href=`/api/report/${c}.csv`;$('expressionTable').innerHTML=e.genes.map(x=>`<tr class="border-b"><td class="p-2 font-medium">${x.gene}</td><td class="p-2">${x.programme}</td><td class="p-2">${x.mean.toFixed(3)}</td><td class="p-2">${x.median.toFixed(3)}</td><td class="p-2">${x.pct_positive.toFixed(1)}%</td></tr>`).join('');const maxMean=Math.max(...e.genes.map(x=>x.mean),1);Plotly.react('dotplot',[{x:e.genes.map(x=>x.programme),y:e.genes.map(x=>x.gene),mode:'markers',type:'scatter',marker:{size:e.genes.map(x=>8+28*x.pct_positive/100),color:e.genes.map(x=>x.mean),colorscale:'Viridis',cmin:0,cmax:maxMean,colorbar:{title:'Mean expression'}},customdata:e.genes.map(x=>[x.mean,x.pct_positive]),hovertemplate:'%{y}<br>mean %{customdata[0]:.3f}<br>positive %{customdata[1]:.1f}%<extra></extra>'}],{margin:{l:90,r:20,t:10,b:50},xaxis:{title:'Owner-supplied programme'},yaxis:{title:'Gene',autorange:'reversed'},height:460})}
async function apply(){try{const mode=$('colourMode').value;if(mode==='cluster'){plot();return}if(mode==='gene'){const d=await get('/api/gene/'+encodeURIComponent($('gene').value.trim()));plot(null,d.values,$('gene').value.trim());return}const d=await get('/api/quality/'+mode);plot(null,d.values,mode)}catch(e){alert(e.message)}}
async function loadAnnotations(){const a=await get('/api/annotations');$('annotationRows').innerHTML=a.rows.map(x=>`<tr class="border-b align-top"><td class="p-2 font-semibold">${x.cluster}</td><td class="p-2">${x.label}</td><td class="p-2">${x.confidence}</td><td class="p-2">${x.dataset}</td><td class="p-2">${x.transcript}</td><td class="p-2">${x.model}</td></tr>`).join('')}
async function loadStability(){const s=await get('/api/stability');$('stabilityMethod').textContent=s.method;$('stabilityRows').innerHTML=s.clusters.filter(x=>['6','7'].includes(String(x.cluster))).map(x=>`<tr class="border-b"><td class="p-2">${x.cluster}</td><td class="p-2">${x.n_cells}</td><td class="p-2">${Number(x.median_same_original_cocluster_frequency).toFixed(3)}</td><td class="p-2">${Number(x.min_same_original_cocluster_frequency).toFixed(3)}–${Number(x.max_same_original_cocluster_frequency).toFixed(3)}</td></tr>`).join('')}
async function init(){data=await get('/api/umap');const clusters=[...new Set(data.cluster)].sort((a,b)=>+a-+b);clusters.forEach(c=>{const b=document.createElement('button');b.textContent='Cluster '+c;b.className='rounded px-2 py-2 text-sm font-medium';b.style.background=cmap[c];b.style.color='white';b.onclick=()=>select(c);$('clusterButtons').appendChild(b)});data.cell_id.slice(0,0);plot();await select('0')}
loadStability().catch(e=>console.warn('Stability summary unavailable',e));loadAnnotations().catch(e=>console.warn('Annotations unavailable',e));$('colourMode').onchange=()=>$('geneControls').classList.toggle('hidden',$('colourMode').value!=='gene');$('colourBtn').onclick=apply;$('resetBtn').onclick=()=>{plot();$('card').classList.add('hidden')};init().catch(e=>alert(e.message));
</script></body></html>'''


@app.get("/", response_class=HTMLResponse)
def index():
    return HTML
