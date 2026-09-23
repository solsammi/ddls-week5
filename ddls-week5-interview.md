# DDLS 2026 — Week 5 interview transcript

**Lab:** Interview a Single-Cell Data Owner
**Data owner:** Dr. Ravi Menon — Immunologist cleaning up a PBMC single-cell dataset before a figure
**Student:** Sammi Baudin <sammi.baudin@scilifelab.se>
**Exported:** 2026-09-23T11:49:18+00:00

> This is the record of your interview with the data owner. Hand it to your analyst agent as the raw source for `AGENTS.md` and `spec.md`.

---

**You**:

What are you trying to find out?

**Dr. Ravi Menon**:

I’m Dr. Ravi Menon. The facility returned my first PBMC single-cell dataset—**2,700 cells in clusters 0–7**—and I’m trying to decide what to do with the tiny cluster 6: **13 cells, about 350 genes each**. I’m about to bin it as dead-cell or empty-droplet junk; can you confirm it’s safe to delete?

**You**:

Some follow up questions to that. What is PBMC? Who collected the data, and how?

**Dr. Ravi Menon**:

PBMC means **peripheral blood mononuclear cells**—mainly lymphocytes and monocytes isolated from blood. I collected the blood sample in our university hospital lab; the core facility isolated the cells, ran the single-cell sequencing, and returned the processed file with the numbered groups and figure. Do you want the experimental provenance, or should we get back to whether cluster 6 is junk?

**You**:

What does the processed file look like? What are the rows and columns?

**Dr. Ravi Menon**:

I haven’t inspected that layout myself. It contains the **individual cells**, their gene readouts, the facility’s cluster number, and coordinates for the 2-D map; I don’t know whether cells are arranged as rows or columns in the file. Is that what you need for the analysis, or are you trying to decide what cluster 6 represents?

**You**:

What does one row represent?

**Dr. Ravi Menon**:

I can’t tell from the summary alone whether one row is a cell or a gene—the facility didn’t explain that to me. It does include the cell readouts and the cluster assignment, though. Can you inspect the column names or file description?

**You**:

What processing was done before the file was handed to you?

**Dr. Ravi Menon**:

I don’t know the exact steps. The core facility did the sequencing and sent back the cells already grouped into clusters 0–7 with a 2-D map and summary numbers; I wasn’t given a methods sheet. Are you asking about sample preparation or how they produced the groups?

**You**:

What are the columns?

**Dr. Ravi Menon**:

I don’t have the full schema in front of me. I know it includes the gene-expression values plus **leiden** (the cluster number), **n_genes**, **total_counts**, and **pct_mito**; there are also the 2-D map coordinates. The facility didn’t tell me whether those are stored as columns or attached cell metadata—can you see the header?

**You**:

Please give me every columns meaning, units and ranges. Do you know if there are any missing values, bias duplicates, or something that has already been cleaned? what do you trust the most or least in this data set?

**Dr. Ravi Menon**:

I can’t give you every column’s meaning because I haven’t been given the full schema. What I know is:

- **Gene-expression values:** one value per gene per cell; I don’t know whether these are raw counts or transformed values, so I can’t state units or ranges.
- **`leiden`:** the facility’s integer group label, 0–7.
- **`n_genes`:** number of genes detected per cell; the summary gives 350 for cluster 6 and 2363 for cluster 7.
- **`total_counts`:** total recorded expression counts per cell; units and range weren’t supplied.
- **`pct_mito`:** mitochondrial-expression percentage; for cluster 6 the summary says 1.6%.
- **Map coordinates:** two numbers per cell for the plotted position; I don’t know their scale or range.

I don’t know whether there are missing values, duplicate cells, batch effects, or what cleaning was done before delivery. I trust the facility’s reported cell and summary numbers as recorded; I trust the biological interpretation of the numbered groups least, because nobody has identified them for me. Can you pull the actual column list and a missing-value/duplicate summary from the file?

**You**:

What is "done" in your eyes?

**Dr. Ravi Menon**:

“Done” means I can give the facility a clear **yes-or-no** on whether to spend the next sequencing run on cluster 7, with a likely identity supported by the genes highlighted on the map and its top genes. I don’t want just “it has 10 cells” or a list of guesses. For cluster 6, I’d like to know whether deleting those **13 cells** is actually safe rather than relying only on the 350-gene figure.

**You**:

What does it mean to "spend the next sequencing run" on cluster 7? What do you mean by "likely identity"?

**Dr. Ravi Menon**:

I mean deciding whether to ask the facility to **prioritise cluster 7 for the next sequencing run**—spend time and budget getting more information from those 10 cells—or leave it out of the figure and use the run elsewhere.

By “likely identity,” I mean the blood-cell type it most resembles—such as a T cell, B cell, NK cell, or monocyte—based on which genes are switched on. Do you want me to provide the practical decision we need, or are you asking about how the facility would technically enrich that group?

**You**:

What is the total number of cluster? Are all cells a part of a cluster?

**Dr. Ravi Menon**:

There are **8 clusters**, numbered 0 through 7, covering all **2,700 cells** in the processed file. Yes, the facility assigned every cell to one of those groups. Are you checking the grouping itself, or deciding whether cluster 7 is worth the next run?

**You**:

Do you want to know something specific about the other clusters?

**Dr. Ravi Menon**:

Not particularly—the immediate decision is **cluster 7**, not a full audit of clusters 0–5. I’d want the genes highlighted for cluster 7 and its top genes, so we can decide whether those **10 cells** have a clear blood-cell identity or are too mixed to justify another run. Can you pull those out?

**You**:

What does top genes mean for you?

**Dr. Ravi Menon**:

The genes showing the strongest or most distinctive expression among the **10 cells in cluster 7**, compared with the other groups. I’d want their names and expression values, not just the cluster size or gene count. Can you extract that list and the genes highlighted on the map?

**You**:

What do you define as a "clear blood cell" identity? What would count as "too mixed"?

**Dr. Ravi Menon**:

A clear identity would be the **same blood-cell programme** across those 10 cells—for example, convincingly T-cell, B-cell, NK-cell, or monocyte-like, rather than just one isolated gene.

“Too mixed” would mean the cells show strong, incompatible programmes at once—say convincing T-cell and monocyte signals together—without one coherent identity. For cluster 7, do its highlighted and top genes point consistently to one lineage, or do they look like several cell types mixed together?

**You**:

do you know the markers of each of these cell types?

**Dr. Ravi Menon**:

Yes—at the bench I’d recognise the usual patterns: **CD3D/E and TRBC** for T cells, **MS4A1 and CD79A** for B cells, **NKG7 and GNLY** for NK cells, and **LYZ with LST1** for monocytes. But I don’t know which of those are present in cluster 7 until you show me the genes; that’s the part I need you to analyse.

**You**:

What would you need to make the yes/no decision about cluster 7? What evidence would point to dead cells in cluster 6?

**Dr. Ravi Menon**:

For **cluster 7**, I’d need the genes highlighted on the map and its top expressed genes, with their expression across the 10 cells. A consistent T-, B-, NK-, or monocyte-like pattern would support spending the run; incompatible patterns or mostly generic signals would argue against it.

For **cluster 6**, I realise I’ve only looked at the summary: 13 cells, about 350 genes each, and 1.6% mitochondrial reads. I hadn’t looked at what they’re actually expressing—I just assumed “hardly any genes” meant dead or empty. Are those cells blank, or are they consistently turning on something specific?

**You**:

Define A consistent T-, B-, NK-, or monocyte-like pattern, also define incompatible pattern and generic signals

**Dr. Ravi Menon**:

A **consistent pattern** would be several genes pointing to the same cell type across most of the 10 cells—for example, CD3D/E and TRBC together across the cells for T cells, or MS4A1 and CD79A together for B cells.

An **incompatible pattern** would be strong T-cell signals alongside strong monocyte signals in the same cells, rather than one coherent blood-cell identity. **Generic signals** would be mostly housekeeping, ribosomal, stress, or mitochondrial genes that don’t identify a blood-cell type. Can you see whether cluster 7 has a coherent set or just scattered signals?

**You**:

Can you list all marker genes, what importance they have and in what cell types they are marker genes?

**Dr. Ravi Menon**:

Not **all**—there isn’t one universal list, and their importance depends on the tissue and assay. The common clues I’d use at the bench are:

- **T cells:** `CD3D`, `CD3E`, `TRBC1/2`
- **B cells:** `MS4A1`, `CD79A`, `CD37`
- **NK cells:** `NKG7`, `GNLY`, `KLRD1`
- **Monocytes:** `LYZ`, `LST1`, `S100A8/A9` for more inflammatory monocytes
- **Dendritic-cell-like cells:** `FCER1A`, `CST3`
- **Platelet contamination:** `PPBP`, `PF4`

Those are clues, not proof, and I don’t know which are actually present in cluster 7. Do you want the facility’s highlighted genes interpreted cell by cell, or just the top-gene list compared with these clues?

**You**:

What markers/conditions would point to a cluster being dead cells or empty droplets?

**Dr. Ravi Menon**:

I’d expect dead or empty material to show **very little expression overall**, with mostly mitochondrial, ribosomal, or general housekeeping genes and no coherent immune-cell pattern. A dying cell might also show a conspicuously large mitochondrial fraction, whereas an empty droplet would have almost no convincing cell-specific signal.

For cluster 6, the summary gives 13 cells, about 350 genes each, and **1.6% mitochondrial reads**—but I haven’t seen the actual genes, so I can’t call it junk from those numbers alone. Can you show what is expressed across those 13 cells?

**You**:

What genes are mitochondrial, ribosomal and housekeeping genes? What conspicuously large mitochondrial fraction?counts as a coherent immune-cell pattern? What counts as a

**Dr. Ravi Menon**:

Examples:

- **Mitochondrial:** `MT-ND1`, `MT-ND2`, `MT-CO1`, `MT-CO2`, `MT-ATP6`, `MT-CYB`
- **Ribosomal:** `RPLP0`, `RPS3`, `RPL13`, `RPS18`
- **Housekeeping:** `GAPDH`, `ACTB`, `B2M`, `UBC`, `HSP90AA1`

There isn’t one fixed mitochondrial percentage that means “dead”; I’d call it conspicuous when it is much higher than the other cells and dominates the expression. A coherent immune pattern means several cell-type-specific genes occur together across the same cells—for example `CD3D/E` with `TRBC` for T cells—not just one isolated transcript. For cluster 6, I only know the summary says **1.6% mitochondrial reads**; I still haven’t seen its actual genes. Can you check whether those 13 cells share a specific immune pattern or mostly show the generic genes above?

**You**:

I will state your own goal and the task back to you as I understood it. Let me know what I miss. You want to know if cluster 7 consists or either Tcells, Bcells, monocytes or NK cells ONLY (no mix), if they do, this warrants another sequencing run. If they dont, another sequencing run is a no. For cluster 6, you want to know if these are dead cells, empty droplets, or something else.

**Dr. Ravi Menon**:

That’s mostly right. For **cluster 7**, I need a clear, coherent likely blood-cell identity—not absolute proof that every one of the 10 cells is identical. If it’s a convincing T-, B-, NK-, or monocyte-like group, I’d consider the next run; if the genes are mixed or uninformative, I’d leave it out.

For **cluster 6**, yes: determine whether the 13 cells are genuinely empty/dead-looking or expressing a consistent immune programme before I delete them. I had been treating the 350 genes per cell as enough evidence, but I haven’t checked the actual expression.

**You**:

What is the number of genes per cell in the other clusters?

**Dr. Ravi Menon**:

The summary lists:

- Cluster 0: **809 genes/cell**
- Cluster 1: **850**
- Cluster 2: **828**
- Cluster 3: **673**
- Cluster 4: **1,263**
- Cluster 5: **1,570**
- Cluster 6: **350**
- Cluster 7: **2,363**

Those are summary values for each cluster, not the individual cells. Are you comparing cluster 6’s expression depth with the others, or checking how unusual cluster 7 is?

**You**:

Is that an average?

**Dr. Ravi Menon**:

I don’t know what summary statistic the facility used—those figures are simply labelled **genes/cell** in the report. I can’t honestly call them averages without their documentation. Can you check the report header or methods note for whether it says mean, median, or something else?
