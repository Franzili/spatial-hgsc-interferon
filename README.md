# spatial-hgsc-interferon

The repository contains everything needed to reproduce the figures and tables of the manuscript _Spatial transcriptomics reveals interferon pathway enrichment in early relapsing serous carcinoma_, starting from Space Ranger output and pathologist histology annotations. The [`quality_control/`](quality_control/README.md) notebooks build the processed `adata.h5ad` object that every downstream analysis uses. If you only want to rerun a specific downstream analysis, download the processed object instead (see [Data availability](#data-availability)).

---

## Study design at a glance

| Axis | Values | Where it lives |
|---|---|---|
| PFI group (`PFI`) | `short`, `medium`, `long` | `adata.obs["PFI"]` — several analyses collapse `medium`+`long` → `long` (`PFI_short_long`) |
| Sampling site (`class` / `class_overall`) | `OvaryL`/`OvaryR` → `Ovary` (adnexa), `Omentum` | `adata.obs["class"]`; `class_overall` merges left/right ovary |
| Tissue compartment (`histology`) | `Tumor Epithelium`, `Other` (stroma), `Mixed TumEpi+Other` | pathologist annotation per spot |
| Patient / sample | `patient`, `sample` (`Paxgene1`–`Paxgene4` TMA slides) | `adata.obs` |

`class == "Marker"` / `patient == "Marker"` spots are orientation markers on the TMA and are
dropped at the start of every downstream analysis.

Two cohorts are analyzed:

- **primary cohort** — the study's own Visium TMAs (`quality_control/primary-cohort/adata.h5ad`).
- **external cohort** — an independent Visium HGSC dataset used for a small case study (`quality_control/external-cohort/pre-processing/adata.h5ad`).

---

## Repository content

```
ovst-staging/
├── quality_control/              # Space Ranger → processed AnnData objects
│   ├── primary-cohort/
│   ├── external-cohort/
├── cell-type_deconvolution/      # cell2location: scRNA reference → spot-level cell-type abundance
├── copy-number-variations/       # inferCNV inference + downstream CNV / PFI association analyses
│   ├── primary-cohort/
│   ├── external-cohort/
│   ├── scripts/                  # CNV extraction scripts
│   └── files/                    # gene_order_file.tsv for inferCNV
├── differential_expression/      # pseudobulk DESeq2 + marker genes; primary & external cohort
├── gsea/                         # pre-ranked GSEA on the DESeq2 results + antigen-presentation analysis
├── ccc/                          # cell–cell communication analysis (LIANA + tensor-cell2cell)
├── src/
│   ├── cnv_inference/            # inferCNV batch-job helper + subcluster extraction (installable package)
│   └── differential_expression/  # run_deseq2.py + SLURM wrapper
├── pyproject.toml / uv.lock      # Python environment (uv)
└── requirements.txt              # pip-installable export of the same environment
```

Each analysis directory has its own `README.md` with inputs, ordered run steps, outputs and the
config variables to edit.

---

## Data availability

| Dataset | Accession | Place at |
|---|---|---|
| Primary cohort | [doi:10.5281/zenodo.22746547](https://doi.org/0.5281/zenodo.22746547) | `quality_control/primary-cohort/adata.h5ad` |
| External cohort (Laury et al., *Modern Pathology* 2024) | [doi:10.1016/j.modpat.2024.100508](https://doi.org/10.1016/j.modpat.2024.100508) | `quality_control/external-cohort/pre-processing/adata.h5ad` |
| scRNA-seq reference (Zheng et al., *Nature Cancer* 2023) | [GSE180661 / doi:10.1038/s43018-023-00599-8](https://doi.org/10.1038/s43018-023-00599-8) | `cell-type_deconvolution/data/OvC_adata_zheng.h5ad` |

> The processed `.h5ad` object is large (2–3 GB each) and is **not** stored in git. Download it from the accession above and put it at the paths shown. For the external cohort, no preprocessed `.h5ad` object is available, but it can be generated from the raw data under the accession above, through the `quality_control` steps for the external cohort in this repo.

### What's in the processed object

`adata` (spots × genes) with:

- `adata.X` / `adata.layers["counts"]` — raw filtered UMI counts
- `adata.layers["log-transformed"]` — `log1p` of library-normalized counts
- `adata.var["gene_ids"]` — Ensembl IDs (`var_names` are gene symbols)
- `adata.obs` — `patient`, `sample`, `class` (in primary cohort only), `PFI` (`outcome` in the external cohort), `histology` (see table above)
- `adata.obsm["spatial"]`, `adata.uns["spatial"]` — Visium coordinates and images per slide

See [`quality_control/README.md`](quality_control/README.md) for the full object specifications.

---

## Environment

**Python** (analysis notebooks and `src/` scripts) — managed with [uv](https://docs.astral.sh/uv/), Python ≥ 3.12:

```bash
uv sync                     # creates .venv/ from pyproject.toml + uv.lock
uv run jupyter lab          # or: source .venv/bin/activate
```

Key packages: `scanpy`, `squidpy`, `anndata`, `cell2location`, `scvi-tools`, `pydeseq2`,
`decoupler`, `gseapy`, `liana`, `cell2cell`, `esda`/`libpysal` (spatial statistics).
The `cell2location` / `scvi` training steps (deconvolution, CCC) need a **CUDA GPU**; every
other step runs on CPU.

**R** (inferCNV only) — a separate environment, not managed here. Needs `infercnv`, `rhdf5`, `optparse`, `jsonlite` (install via `BiocManager`).
`copy-number-variations/scripts/run_infercnv_sample.R` will `BiocManager::install()` anything missing on first run. The one-off Seurat → AnnData reference conversion (`cell-type_deconvolution/load-reference.R`) additionally needs `Seurat`, `SeuratDisk`, `sceasy`, `reticulate`.

**HPC** — the long-running jobs (cell2location training, inferCNV, DESeq2) were run on a SLURM cluster (CSC Roihu). The `*.sh` / `*_batch.sh` scripts and `config.json` files are templates.

---

## Analysis pipeline

Run in dependency order. Steps 3–5 are independent of each other once the object is in place.

| # | Analysis | Directory | Needs |
|---|---|---|---|
| 1 | Obtain data | — | Space Ranger output + annotation CSVs, or the processed object ([Data availability](#data-availability)) |
| 2 | Quality control / pre-processing | [`quality_control/`](quality_control/README.md) | Space Ranger output (skip if you downloaded the processed object) |
| 3 | Cell-type deconvolution | [`cell-type_deconvolution/`](cell-type_deconvolution/README.md) | primary object + scRNA reference |
| 4 | Copy-number variations | [`copy-number-variations/`](copy-number-variations/README.md) | cohort object (+ inferCNV R env) |
| 5 | Differential expression | [`differential_expression/`](differential_expression/README.md) | cohort object |
| 6 | GSEA | [`gsea/`](gsea/README.md) | DESeq2 result CSVs from step 5 |
| 7 | Cell–cell communication | [`ccc/`](ccc/README.md) | primary object + cell2location output from step 3 |
| 8 | Antigen-presentation analysis | [`gsea/`](gsea/README.md) (`antigen-presentation.ipynb`) | GSEA results from step 6 + primary object |

Notebooks and scripts carry a **configuration block near the top** (`ADATA`, `DATA_DIR`, `WORKING_DIR`,
`OUT_DIR`, `*_FILE`, ...). These are **repo-relative paths** and need no editing as long as you keep the
repository layout and place the downloaded objects where [Data availability](#data-availability) says.
Edit them only if your data lives elsewhere. Each sub-README lists the specific variables per file.

> **Run notebooks from the directory they live in** (the default for `jupyter lab` / VS Code) — the
> relative paths are resolved from there. The `src/` scripts and the `*.sh` wrappers instead resolve
> the repository root from their own location, so they can be launched from anywhere.

---

## License

[MIT](LICENSE).

## Citation

If you use this code, please cite the manuscript (see [`CITATION.cff`](CITATION.cff)) and the tools it builds on (cell2location, inferCNV, PyDESeq2, decoupler, LIANA, Tensor-cell2cell, GSEApy, scanpy/squidpy).
