# ccc — differential cell–cell communication

`Differential-CCC.ipynb` — does ligand–receptor signalling between cell types differ between
PFI groups? Combines [LIANA](https://liana-py.readthedocs.io/) (per-sample LR scoring) with
[Tensor-cell2cell](https://earmingol.github.io/cell2cell/) (context-aware tensor
decomposition), then interprets the PFI-associated factors with PROGENy and GSEA.
Runs on a **GPU** (the scVI batch-correction step).

## Input

A spatial object **already carrying the cell2location deconvolution results** — i.e. the mapping output from [`cell-type_deconvolution/`](../cell-type_deconvolution/README.md) (`spatial.h5ad`, with `obsm["q05_cell_abundance_w_sf"]`, `obsm["means_cell_abundance_w_sf"]`, `uns["mod"]["factor_names"]`). The notebook expects it at `../cell-type_deconvolution/data/adata.h5ad` (`DATA_DIR` / `ADATA` variables) — copy or symlink the mapping output there.

## What the notebook does

1. **Niche assignment** — `cell2location.run_colocation` (NMF, `n_fact` 10–15) → dominant
   colocalisation compartment per spot (`major_cell_type`: Cancer, Fibroblast, Epithelial,
   Myeloid, B cell, T cell, ...).
2. **QC + normalisation** — spot filtering (`min_counts=250`, `min_genes=50`), then scVI
   batch correction across `patient` (`n_latent=30`, NB likelihood); `adata.X` ← batch-
   corrected, log- and library-normalised expression.
3. **LIANA** — `li.mt.rank_aggregate.by_sample`, grouped by `major_cell_type`, per
   `sample_context = class × PFI × patient`, `consensus` resource, `expr_prop=0.9`.
   → `LIANA_by_sample.csv`.
4. **Tensor-cell2cell** — `to_tensor_c2c` (score = `magnitude_rank`), metadata by PFI,
   `run_tensor_cell2cell_pipeline` across ranks; **rank 12** used for the final decomposition.
   → `tensor_factorized.pkl`, `tensor_metadata.pkl`.
5. **Factor ↔ PFI** — Spearman correlation of each context factor with PFI severity
   (long=1, medium=2, short=3): **Factor 1 ↔ short**, **Factor 7 ↔ long**. Context boxplots,
   loading clustermaps, per-factor CCC networks / circos plots.
6. **Interpretation** — PROGENy pathway footprint enrichment on the LR loadings
   (`decoupler.mt.ulm`); GSEA (`cell2cell.external.run_gsea`, GO-BP) on the LR loadings;
   Factor 1 vs Factor 7 contrast.

## Outputs

Written to `OUT_DIR` (`../differential_expression/figures/` as configured):
`processed.h5ad`, `LIANA_by_sample.csv`, `tensor_factorized.pkl`, `tensor_metadata.pkl`, `tensor_factors.pdf`, `context_boxplots.pdf`, `Clustermap-Contexts.pdf`, `comm_networks_pfi.pdf`, `PROGENy.pdf`, `GSEA-Dotplot.pdf`, and per-slide niche maps.

## Config variables to edit

| Variable | Meaning |
|---|---|
| `WORKING_DIR` | output root (currently `differential_expression/`) |
| `DATA_DIR` | folder holding the deconvolved `adata.h5ad` |
| `ADATA` | deconvolved spatial object (`spatial.h5ad` from the mapping step) |
| `OUT_DIR` | figure / results directory |
| `n_fact`, tensor `rank` | colocalisation factor count and tensor decomposition rank |
