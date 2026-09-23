# AGENTS.md

## Environment

- Use `uv` for the Python environment: create it with `uv venv`.
- Run all Python with `uv run`; this works the same on every OS.
- All analysis code must be created as Python scripts so a human can review it by eye.
- Keep conventions and assumptions consistent across scripts.

## Data and outputs

- The supplied data archive is `ddls-week5-s1-junk-or-signal-dataset.zip`.
- It contains `data/pbmc3k.h5ad` and `data/ABOUT_THIS_FILE.txt`.
- The `.h5ad` file is a processed single-cell RNA-seq dataset of 2,700 human PBMCs and 13,714 genes. It contains log-normalised expression in `ad.X`, raw UMI counts in `ad.layers["counts"]`, gene symbols in `ad.var_names`, cluster labels in `ad.obs["leiden"]`, per-cell quality fields in `ad.obs`, and UMAP coordinates in `ad.obsm["X_umap"]`.
- Clusters and the map are already computed; do not recompute them unless the owner later asks for that.
- Write analysis outputs to `results/`.

## Version control

- This folder is now a git repository.
- Commit the current state before any big change.
- Commit again whenever something starts working.
- Use short, clear commit messages.
- Never commit `.env`, `.venv/`, `__pycache__/`, or `*.pyc`.

## Evidence and confidence

- Never report an answer about a cluster without first reporting the confidence that matches the claim.
- Separate what the files show from interpretation and uncertainty.
- Do not claim missing schema, processing, or biological facts that are not supported by the transcript or files.
