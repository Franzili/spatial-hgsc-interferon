# differential_expression

Differential gene expression between PFI groups, histology compartments and sampling sites, in the primary and external cohorts. Two approaches:

- **pseudobulk DESeq2** (`src/differential_expression/`) — the quantitative pipeline; its CSVs feed [`gsea/`](../gsea/README.md).
- **marker genes** (`*/marker_genes.ipynb`) — scanpy `rank_genes_groups` for descriptive marker matrixplots.

```
differential_expression/
├── primary-cohort/
│   ├── deseq2.ipynb              figures from the DESeq2 result CSVs
│   ├── marker_genes.ipynb        scanpy marker genes (tumour epithelium; adnexa; omentum)
│   ├── tumor-stroma-border.ipynb C3 / IFI27 / BST2 per-patient boxplots by histology × PFI
│   ├── output/                   marker-gene matrixplots (not in git)
│   └── figures/                  deseq2.ipynb figures + result CSVs (not in git)
├── external-cohort/
│   ├── deseq2.ipynb              figures from the external-cohort DESeq2 results
│   ├── marker_genes.ipynb
│   └── DEGs-spatially.ipynb      exploratory pydeseq2 notebook (see below)
└── ../src/differential_expression/
    ├── run_deseq2.py             pseudobulk build + all DESeq2 contrasts
    └── run_deseq2.sh             SLURM wrapper
```

---

## 1. Pseudobulk DESeq2 — `src/differential_expression/run_deseq2.py`

Aggregates spot counts (`layers["counts"]`) to **pseudobulk samples** by `PFI_short_long × histology × class_overall × patient` (`medium`+`long` PFI merged to `long`; left/right ovary merged to `Ovary`; `Marker` spots dropped). Genes are filtered (`≥ min_count` in `≥ min_samples` pseudobulk samples, plus a low-CV filter), then [PyDESeq2](https://pydeseq2.readthedocs.io/) fits the following contrasts:

| Contrast name | Design | Question |
|---|---|---|
| `condition_short_vs_long` | `~ histology + condition` | overall PFI effect |
| `condition_in_histology_<h>` | `~ condition` (subset to one histology) | PFI effect within tumour epithelium / stroma / mixed |
| `condition_in_class_<c>` | `~ histology + condition` (subset to one site) | PFI effect within adnexa / omentum |
| `condition_in_class_<c>_tumor_epithelium` | `~ condition` (tumour epithelium only, one site) | site-specific PFI effect, tumour-only (unconfounded by tumour/stroma ratio) |
| `histology_<h1>_vs_<h2>` | `~ condition + histology` | pairwise histology |
| `class_<c1>_vs_<c2>` | `~ histology + class` | pairwise site |

Per contrast it writes `deseq2_results_batch/`:
`<name>_full.csv` (unshrunk — use its `padj` for significance),
`<name>_shrunk.csv` (apeglm-shrunk `log2FoldChange` — use for plotting),
`<name>_significant.csv` (`padj < 0.05` & `|log2FC| > 1`),
`<name>_dds.pkl` (fitted model), plus `summary.csv`, `pseudobulk_counts.csv`,
`pseudobulk_metadata.csv`.

### Run

```bash
# HPC:
sbatch src/differential_expression/run_deseq2.sh

# or locally / interactively:
uv run python src/differential_expression/run_deseq2.py \
  --anndata quality_control/primary-cohort/adata.h5ad \
  --outdir  differential_expression/primary-cohort/deseq2_results_batch \
  --min-count 10 --min-samples 2 --padj 0.05 --lfc 1.0 --ncpus 8
```

> Point `--outdir` at `differential_expression/primary-cohort/deseq2_results_batch/` — that is
> where `deseq2.ipynb` and `gsea/gsea.ipynb` (`DEG_DIR`) look for the CSVs.
> `run_deseq2.sh` resolves `REPO_ROOT` from its own location; set `#SBATCH --account` to your
> allocation, and export `ENV_PATH` if your cluster needs a specific Python environment
> (`ENV_PATH=/path/to/env/bin sbatch src/differential_expression/run_deseq2.sh`).

### 2. Figures — `primary-cohort/deseq2.ipynb`

Loads the CSVs above and produces: DEG-count summary bar, MA plots, volcano plots,
effect-comparison scatters (PFI effect: adnexa vs omentum; tumour vs stroma; all histology
pairs), top-DEG tables, per-contrast gene lists, and a pseudobulk clustermap.
Genes of interest highlighted throughout: **C3, IFI27, BST2**. Config at top: `OUT_DEGS`,
`FIG_DIR`, `PADJ`, `LFC`, `PFI_PALETTE`.

---

## 3. Marker genes — `*/marker_genes.ipynb`

`scanpy` `highly_variable_genes(flavor="seurat_v3")` + `rank_genes_groups(method="wilcoxon")`
on `PFI`, restricted to tumour-epithelium spots, then repeated for adnexa (`class_overall ==
"Ovary"`) and omentum. Produces marker matrixplots → `primary-cohort/output/*.{pdf,png}`.

| Notebook | Input object (`ADATA` variable) |
|---|---|
| `primary-cohort/marker_genes.ipynb` | `../../quality_control/primary-cohort/adata.h5ad` |
| `external-cohort/marker_genes.ipynb` | `../../quality_control/external-cohort/pre-processing/adata.h5ad` |

---

## 4. Supporting notebooks

- **`primary-cohort/tumor-stroma-border.ipynb`** — per-patient mean-expression boxplots of
  C3 / IFI27 / BST2 by histology × PFI (with per-site marker shapes), plus per-compartment `rank_genes_groups` stats.
- **`external-cohort/DEGs-spatially.ipynb`** — earlier PyDESeq2 exploration
  (includes a `class:condition` interaction design). **Superseded by `run_deseq2.py`**; kept
  for provenance.

---

## Config variables to edit

| File | Variables |
|---|---|
| `run_deseq2.py` | CLI args, or `ANNDATA_PATH` / `OUT_DEGS` env vars |
| `run_deseq2.sh` | `#SBATCH --account`; `ENV_PATH` (env var, optional); `OUT_DEGS` |
| `deseq2.ipynb` | `OUT_DEGS`, `FIG_DIR` |
| `marker_genes.ipynb` | `WORKING_DIR`, `OUT_DIR`, `ADATA` |
| `tumor-stroma-border.ipynb` | `WORKING_DIR`, `OUT_DIR`, `INPUT_H5AD` |
