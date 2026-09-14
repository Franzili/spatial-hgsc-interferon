# gsea

Pathway-level interpretation of the PFI differential expression signal.

- `gsea.ipynb` — pre-ranked GSEA on the DESeq2 contrasts.
- `antigen-presentation.ipynb` — focused analysis of the MHC-I / MHC-II antigen-presentation
  module within the GSEA leading edges.

Both need internet access on first run to fetch gene-set libraries from
[Enrichr](https://maayanlab.cloud/Enrichr/) (cached locally afterwards).

---

## `gsea.ipynb` — pre-ranked GSEA

**Input:** the DESeq2 result CSVs from
[`differential_expression/`](../differential_expression/README.md), via
`DEG_DIR = "../differential_expression/primary-cohort/deseq2_results_batch"` (and `summary.csv`
in that folder). Run the DESeq2 pipeline first — including the
`condition_in_class_<site>_tumor_epithelium` contrasts, which family 2b depends on.

**Method:** genes are ranked per contrast (signed by `log2FoldChange`, weighted by
significance), then `gseapy.prerank` is run against **MSigDB Hallmark 2020** (1000 permutations,
gene-set size 15–500). FDR < 0.25 is the significance cut for the figures.

**Three analysis families**, run and plotted independently:

| Family | Contrasts | Question |
|---|---|---|
| 1 — overall | `condition_short_vs_long` | Does the PFI effect (short vs medium+long) show up pooling all samples? |
| 2 — by site | `condition_in_class_Ovary` / `_Omentum` | Is the PFI effect consistent between adnexa and omentum? |
| 2b — by site, tumour-epithelium only | `condition_in_class_*_tumor_epithelium` | …controlling for differing tumour/stroma proportions between sites. |
| 3 — by histology | `condition_in_histology_Tumor Epithelium` / `_Other` | Is the PFI effect consistent between tumour epithelium and stroma? |

**Outputs** (`GSEA_DIR = "gsea_results"`, `FIG_DIR = "figures"`):
`gene_ranking_overall.csv`, `gsea_pfi_overall_MSigDB_Hallmark_2020.csv` (+ per-family GSEA
tables), the downloaded `.gmt`, and NES heatmaps / bar charts / divergence scatters.

---

## `antigen-presentation.ipynb` — MHC module analysis

**Inputs:**
- `gsea_results/gsea_pfi_overall_MSigDB_Hallmark_2020.csv` and
  `gsea_results/gene_ranking_overall.csv` (from `gsea.ipynb`)
- `../quality_control/primary-cohort/adata.h5ad` (primary object; see the
  [quality_control README](../quality_control/README.md) path note)
- `MSigDB_Hallmark_2020.tsv` — auto-downloaded from Enrichr and cached

**Method:** a curated MHC-I / peptide-loading and MHC-II gene module is intersected with the
Hallmark Interferon-γ / Interferon-α leading edges; a **detection-bounded depletion test**
(Fisher exact, accounting for which module genes are on the Visium probe set) asks whether
antigen-presentation genes are under-represented in the leading edge, with a sensitivity
analysis around `HLA-B`.

**Outputs:** leading-edge composition and test tables (`suppl_T1_leading_edge_composition.csv`
and printed summaries).

---

## Config variables to edit

| File | Variables |
|---|---|
| `gsea.ipynb` | `DEG_DIR`, `GSEA_DIR`, `FIG_DIR`, `FIG_FORMAT`, `FIG_DPI`, `GENE_SETS` |
| `antigen-presentation.ipynb` | `GSEA_CSV`, `RANKING`, `LIB_CACHE`, the `adata` path in the load cell |
