# copy-number-variations

Copy-number-variation (CNV) inference from spatial transcriptomics with
[inferCNV](https://github.com/broadinstitute/inferCNV), followed by downstream analyses linking
CNV burden to platinum-free interval (PFI). Run for both the **primary** and the
**external** cohort (`primary-cohort/`, `external-cohort/`).

Pipeline: **prepare inputs → run inferCNV (SLURM) → downstream CNV/PFI analyses**.

```
primary-cohort/ , external-cohort/
├── notebooks/                prepare_annotations*.ipynb  +  downstream analysis notebooks
├── config.json              (primary)  inferCNV + SLURM parameters
├── cnv_inference/
│   ├── infercnv_input/       per-sample / combined counts (.h5) + annotation (.tsv)
│   │   └── */config.json    (external) inferCNV + SLURM parameters
│   ├── infercnv_runs/        inferCNV output (HMM predictions, heatmaps) — not in git
│   └── batch_*.sh           (external) ready-made SLURM array scripts
└── outputs/                  downstream results (gene lists, association tables) — not in git
scripts/                      standalone extraction scripts (both cohorts)
files/gene_order_file.tsv     genome-wide gene ordering for inferCNV (GRCh38)
../src/cnv_inference/          infercnv_helper.py (batch-script generator) + infercnv_subclusters.py
```

---

## 1. Prepare inferCNV inputs

`notebooks/prepare_annotations.ipynb` (primary cohort) and the external-cohort variants
(`prepare_annotations.ipynb`, `prepare_annotations_stroma.ipynb`,
`prepare_annotations_tumor-stroma-OVA09.ipynb`) take the corresponding `adata.h5ad` and write, per sample and/or across samples:

- `<sample>_counts.h5` — raw counts (`layers["counts"]`), genes keyed by **Ensembl ID**
- `<sample>_infercnv_annotations.tsv` — barcode → histology group

into `cnv_inference/infercnv_input/{per_sample,combined,...}/`. Stroma / `Other` spots are
used as the **normal reference**; necrotic-debris, blood and surface-epithelium spots are
excluded (external cohort). Edit `WORKING_DIR`, `OUT_DIR`, `ADATA` at the top of the notebook.

---

## 2. Run inferCNV

For each sample inferCNV is run by `scripts/run_infercnv_sample.R` (reads the `.h5` + `.tsv`,
builds the object, runs `infercnv::run(..., HMM=TRUE)`, exports the plotted matrix and HMM
predictions). Samples are submitted as a **SLURM array**, one task per sample.

### 2a. Primary cohort — generate the batch script with `infercnv_helper.py`

`config.json` (in `primary-cohort/`) is the main configuration for an infercnv run. Main sections:

- **`input_dir`** – directory with the per-sample `<sample>_counts.h5` / `<sample>_infercnv_annotations.tsv`
- **`gene_order_file`** – `files/gene_order_file.tsv`
- **`output_root`** – where inferCNV run directories are written
- **`infercnv`** – default parameters (`cutoff`, `denoise`, `cluster_by_groups`, `num_threads`)
- **`samples`** – per-sample overrides; each key must match a discovered input file stem and
  may set `HMM`, `ref_groups`, `cutoff`
- **`slurm`** – account, partition, resources, R module, optional `project_r_lib`

```bash
python ../src/cnv_inference/infercnv_helper.py \
  --config primary-cohort/config.json \
  --r-script copy-number-variations/scripts/run_infercnv_sample.R \
  --manifest primary-cohort/cnv_inference/infercnv_runs/manifest.tsv \
  --slurm-script primary-cohort/cnv_inference/submit_infercnv_array.sh

sbatch primary-cohort/cnv_inference/submit_infercnv_array.sh
```

This writes a tab-separated **manifest** (one sample per row: `sample`, `counts_h5`,
`annotations_tsv`, `out_dir`, `cutoff`, `hmm`, `ref_groups`) and the SLURM array script.
The manifest can also be written by hand.

### 2b. External cohort — pre-made scripts

`external-cohort/cnv_inference/batch_tumor-stroma.sh` and `batch_tumor-stroma-ova09.sh` are the equivalent array scripts. `MANIFEST`, `R_SCRIPT` and `GENE_ORDER` are derived from a `REPO_ROOT` that the script resolves from its own location, so only `#SBATCH --account` needs to be allocated. A cluster-specific R library can be passed optionally:

```bash
PROJECT_R_LIB=/path/to/your/R/library sbatch batch_tumor-stroma.sh
```

Their `config.json` files sit in `external-cohort/cnv_inference/infercnv_input/{combined,tumor-stroma-ova09}/`.

### inferCNV output files referenced downstream

| Cohort | Prefix of HMM prediction files |
|---|---|
| primary | `HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_{genes,regions}.dat` |
| external | `17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_{genes,regions}.dat` |

plus `infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observation_groupings.txt`
(barcode → subcluster).

---

## 3. Downstream analyses

All read the inferCNV HMM output + the cohort `adata.h5ad`. Set the `DATA_DIR` /
`*_FILE` / `ADATA_PATH` variables in the config cell at the top of each notebook.

### Scripts (`scripts/`)

| Script | Purpose |
|---|---|
| `extract_spots_with_cnv_gene.py` | Extract Visium spots carrying an amplification or deletion of a gene of interest (built-in table of ~40 HGSC/immune genes, or supply an Ensembl ID / coordinates). Example: `BST2` amplification → `primary-cohort/outputs/BST2_results/`. Usage examples are in the module docstring. |
| `extract_genes_altered.py` | Recurrently altered genes across patients (per-gene patient count, dominant direction/state, top genes per chromosome). Writes `outputs/recurrent_altered_genes.tsv`, `top_genes_per_chromosome.tsv`. Edit the `sample_to_file` dict and output paths at the bottom. |
| `../src/cnv_inference/infercnv_subclusters.py` | Map inferCNV subclusters (the horizontal bands in the heatmap) to their CNV fingerprints and to AnnData barcodes; plot per-subcluster CNV state maps. Edit `INFERCNV_DIR` / `ADATA_PATH` in the `__main__` block. |

### Notebooks (`{cohort}/notebooks/`)

| Notebook | Question |
|---|---|
| `genome_wide_cnv_scan.ipynb` | Which genomic regions (any chromosome, gains and losses) have CNV burden enriched in short-PFI vs long-PFI patients? Per-sample per-100 kb-bin gain/loss fraction, Mann–Whitney + BH correction. |
| `chr19_mcr_analysis.ipynb` | Minimal common region on chr19 consistently amplified in short-PFI patients. |
| `chr19_amplicon_pfi_analysis.ipynb` | Is the chr19 amplicon *as a whole* (not any single gene) associated with PFI category? Ordinal logistic + Kruskal/Mann–Whitney. |
| `cnv_pfi_category_analysis.ipynb` | Association between a specific CNV (e.g. `BST2` amplification, from `extract_spots_with_cnv_gene.py`) and PFI category, per sample. |
| `prepare_annotations*.ipynb` | (step 1 above) |
| `tumor-stroma.ipynb` (external) | Tumour vs stroma inferCNV comparison for the OVA09 sample. |

`primary-cohort/notebooks/chr19_mcr_analysis_backup.ipynb` is an earlier draft of
`chr19_mcr_analysis.ipynb` and can be ignored.

---

## Config / environment notes

- **R environment** for inferCNV: `infercnv`, `rhdf5`, `optparse`, `jsonlite` (via
  `BiocManager`). `run_infercnv_sample.R` installs missing packages on first run.
- `PFI_PALETTE = {'short': '#C7844A', 'medium': '#456EAE', 'long': '#538984'}` is shared
  across the downstream notebooks.
- SLURM scripts and `config.json` ship with the placeholder account `project_xxxxxxx` — replace it with your own allocation before submitting. `project_r_lib` in `config.json` is empty by default; set it if your cluster needs a project-specific R library path. All other paths are repo-relative.
