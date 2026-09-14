# cell-type_deconvolution

Spot-level cell-type deconvolution of the primary-cohort Visium data with
[cell2location](https://cell2location.readthedocs.io/), using an external ovarian-cancer
scRNA-seq atlas (Zheng et al., *Nat Cancer* 2023) as the reference.

Outputs (per-spot cell-type abundances and colocalisation factors) are written into the spatial
object and are reused by the [`ccc/`](../ccc/README.md) analysis.

## Inputs

| Input | Path | Source |
|---|---|---|
| Primary-cohort spatial object | `../quality_control/primary-cohort/adata.h5ad` | [data deposition](../README.md#data-availability) |
| scRNA-seq reference | `data/OvC_adata_zheng.h5ad` | GSE180661 / doi:10.1038/s43018-023-00599-8 |

Create `cell-type_deconvolution/data/` and place the reference there (or point `DATA_DIR` /
`ADATA` at your own locations).

## Steps

| Order | File | Runs on | What it does |
|---|---|---|---|
| 0 (optional) | `load-reference.R` | R | Convert a Seurat `.rds` reference to `.h5ad` via `sceasy`, if your reference is distributed as a Seurat object. Skip if you already have `OvC_adata_zheng.h5ad`. |
| 1 | `prepare_reference_zheng.ipynb` | CPU | Subset the atlas to HGSOC patients and primary + metastatic tumour, relabel cell types, drop rare/ambiguous types (DC, HSC, Monocyte, Mesothelial). Writes `data/OvC_adata_zheng_preprocessed.h5ad`. |
| 2 | `c2l_training_zheng.py` (submit with `c2l_batch.sh`) | **GPU** | Train the cell2location `RegressionModel` (negative-binomial regression, `batch_key="Patients"`, `labels_key="cell_type"`, 300 epochs) to estimate reference signatures. Writes `spot-deconvolution_zheng/reference_signatures/` (`adata_sc.h5ad`, `adata_st.h5ad`, model, QC plots). |
| 3 | `c2l_spatial_mapping_zheng.py` (submit with `c2l_mapping_batch.sh`) | **GPU** | Fit the `Cell2location` spatial model (`N_cells_per_location=15`, `detection_alpha=20`, 30 000 epochs) and export the posterior. Writes `spot-deconvolution_zheng/<run_name>/spatial.h5ad` with `obsm["q05_cell_abundance_w_sf"]`, `obsm["means_cell_abundance_w_sf"]`. |
| 4 | `cell2location_spatial_mapping_zheng.ipynb` | CPU | Downstream: assign the dominant cell type per spot; cell-type composition by site / histology / PFI (heatmaps, stacked bars); Leiden region ("niche") clustering on the abundance matrix; UMAPs; cell-type colocalisation (correlation network + `cell2location.run_colocation` NMF factors). Figures → `visualizations/`. |

`cell2location_training_zheng.ipynb` is the interactive equivalent of step 2 (250 epochs);
`c2l_training_zheng.py` / `c2l_spatial_mapping_zheng.py` are the batch versions actually used.

> **Path caveat:** the mapping script writes to `.../cell2location_map_reg20/` while the
> downstream notebook reads `.../cell2location_map/`. Make the `run_name` in
> `c2l_spatial_mapping_zheng.py` and the notebook match, or rename the output directory.
> A second file, `cell2location_spatial_mapping_zheng.py`, is an earlier variant
> (`detection_alpha=200`) kept for reference — `c2l_spatial_mapping_zheng.py` is the current one.

## Config variables to edit

Every script/notebook has these near the top:

| Variable | Meaning |
|---|---|
| `WORKING_DIR` | this directory on your system |
| `DATA_DIR` | where `OvC_adata_zheng*.h5ad` live |
| `ADATA` | primary-cohort spatial object (`../quality_control/primary-cohort/adata.h5ad`) |
| `OUT_DIR` / `results_folder` | where model outputs go (`spot-deconvolution_zheng/`) |

`c2l_batch.sh` / `c2l_mapping_batch.sh`: set `#SBATCH --account` to your own allocation and adjust `module load` for your cluster. The environment is activated from `$ENV_PATH`, so export it before submitting:

```bash
ENV_PATH=/path/to/your/env/bin sbatch c2l_batch.sh
```

## Outputs

- `spot-deconvolution_zheng/reference_signatures/` — trained regression model + signatures
- `spot-deconvolution_zheng/<run>/spatial.h5ad` — spatial object with cell-type abundances
- `visualizations/` — UMAPs, composition bar charts, colocalization heatmap (not in git)
