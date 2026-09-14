# quality_control

Builds the **processed spatial transcriptomics objects** that every downstream analysis in
this repository starts from, from Space Ranger output. One notebook per cohort.

```
quality_control/
├── primary-cohort/quality_control.ipynb          primary cohort → quality_control/primary-cohort/adata.h5ad
└── external-cohort/pre-processing/pre-processing.ipynb   external cohort → quality_control/external-cohort/pre-processing/adata.h5ad
```

## Inputs

| Input | Notes |
|---|---|
| Space Ranger `filtered_feature_bc_matrix.h5` + `spatial/` (per capture area) | aggregated run over the 4 TMA slides; from the raw-data deposition ([top-level README](../README.md#data-availability)) |
| Annotation CSVs (`DATA_DIR/annotations/aggr/`) | primary: `PAX_annotations.csv` (histology), `PAX_patients.csv` (patient), `PAX_class.csv` (anatomical site); external: `TissueType.csv`, `Patient-ids.csv`, `OutcomeGroups.csv` |

Barcode suffixes `-1..-4` map to slides `Paxgene1..4` (primary) / `HC-TMA1..4` (external).
PFI group (primary) / outcome (external) is joined per patient.

## What the notebooks do

1. Load the aggregated matrix, attach per-slide spatial images (`squidpy`), merge the
   annotation CSVs into `obs`.
2. QC metrics (`sc.pp.calculate_qc_metrics`; no mitochondrial metrics — the v1 Visium FFPE
   probe set has no MT probes). Diagnostic plots: joint `total_counts` vs `n_genes_by_counts`,
   violins, per-slide spatial maps.
3. **Spot filtering** — flag and drop spots with `total_counts < 150` (`qc_lib_size`) or
   `n_genes_by_counts < 100` (`qc_expressed_genes`); `discard = qc_lib_size | qc_expressed_genes`.
4. **Gene filtering** — `sc.pp.filter_genes(min_cells=3)`. Save `adata_filtered_raw.h5ad`
   (raw counts, filtered).
5. **Normalisation** — store `layers["counts"]` (raw), then add
   `layers["library-nomalized"]` (`normalize_total`, target 1e4, `exclude_highly_expressed`),
   `layers["log-transformed"]` (`log1p`), `layers["scaled"]` (`scale`, zero-centred),
   `layers["pearson"]` (`normalize_pearson_residuals`).
6. Write `adata.h5ad`.

## Output object structure

Spots × genes `AnnData`:

- `X` — raw filtered counts (identical to `layers["counts"]`)
- `layers["counts"]` — raw integer UMI counts (used by DESeq2 and inferCNV)
- `layers["log-transformed"]` — `log1p` of library-normalised counts (used by scanpy marker genes)
- `layers["library-nomalized"]`, `layers["scaled"]`, `layers["pearson"]` — alternative normalisations
- `var_names` — gene symbols; `var["gene_ids"]` — Ensembl IDs
- `obs` columns:
  | Column | Values | Notes |
  |---|---|---|
  | `patient` | patient ID (`H###`, primary / `OVA##`, external) | `"Marker"` for TMA orientation-marker spots |
  | `sample` | `Paxgene1`…`4` / `HC-TMA1`…`4` | TMA slide / Visium capture area |
  | `class` | `OvaryL`, `OvaryR`, `Omentum`, `Marker` | anatomical origin (primary only) |
  | `PFI` | `short`, `medium`, `long` | primary cohort |
  | `outcome` | `short`, `long` | external cohort (analogue of `PFI`) |
  | `histology` | `Tumor Epithelium`, `Other`, `Mixed TumEpi+Other` (primary); tissue-type labels (external) | per-spot pathologist annotation |
- `obsm["spatial"]`, `uns["spatial"]` — coordinates + H&E images, keyed by `sample`

## Derived columns added by downstream analyses

Most downstream notebooks re-create these on load:

```python
adata.obs["class"]          = adata.obs["class"].astype(str)
adata.obs["class_overall"]  = adata.obs["class"].replace({"OvaryR": "Ovary", "OvaryL": "Ovary"})
adata.obs["PFI_short_long"] = adata.obs["PFI"].replace({"medium": "long", "long": "long"})
adata = adata[~adata.obs["class"].isin(["Marker"])].copy()   # drop orientation markers
```

## Config variables to edit

All paths are **repo-relative**; run each notebook from its own directory.

| Variable | Default | Meaning |
|---|---|---|
| `WORKING_DIR` | `"."` | output root — where `adata.h5ad` is written |
| `DATA_DIR` | `"data"` | Space Ranger output + `annotations/` (not in git — see [Data availability](../README.md#data-availability)) |
| `OUT_DIR` | `"figures"` | subdirectory for QC plots |
