import pandas as pd
from pathlib import Path
import scanpy as sc

def load_all_hmm_pred_genes(sample_to_file):
    """
    sample_to_file: dict mapping patient/sample ID -> pred_cnv_genes.dat path
    """
    dfs = []
    for sample, path in sample_to_file.items():
        df = pd.read_csv(path, sep="\t")
        df["sample"] = sample
        dfs.append(df)
    out = pd.concat(dfs, ignore_index=True)

    out["gene"] = out["gene"].astype(str)
    out["sample"] = out["sample"].astype(str)
    out["state"] = pd.to_numeric(out["state"])

    return out

def get_recurrent_altered_genes(
    hmm_df,
    neutral_state=3,
):
    """
    Count genes altered in how many patients.
    A gene counts once per patient if it is altered in any subcluster there.
    """
    altered = hmm_df[hmm_df["state"] != neutral_state].copy()

    # count each gene once per patient
    gene_patient = altered[["gene", "sample"]].drop_duplicates()

    recurrence = (
        gene_patient
        .groupby("gene")
        .agg(
            n_patients_altered=("sample", "nunique"),
            patients=("sample", lambda x: ",".join(sorted(set(x))))
        )
        .reset_index()
    )

    n_total_patients = hmm_df["sample"].nunique()
    recurrence["fraction_patients_altered"] = (
        recurrence["n_patients_altered"] / n_total_patients
    )

    recurrence = recurrence.sort_values(
        ["n_patients_altered", "fraction_patients_altered", "gene"],
        ascending=[False, False, True]
    ).reset_index(drop=True)

    return recurrence


def get_recurrent_altered_genes_with_meta(
    hmm_df,
    neutral_state=3,
):
    """
    Compute recurrent altered genes across patients, including:
      - chromosome
      - dominant state
      - direction (loss/gain)
      - state distribution
    """

    df = hmm_df.copy()
    df["gene"] = df["gene"].astype(str)
    df["sample"] = df["sample"].astype(str)
    df["state"] = pd.to_numeric(df["state"])

    # keep only altered
    df = df[df["state"] != neutral_state].copy()

    # direction
    df["direction"] = "neutral"
    df.loc[df["state"] < neutral_state, "direction"] = "loss"
    df.loc[df["state"] > neutral_state, "direction"] = "gain"

    # ensure gene-chr consistency (take first non-null)
    gene_chr = (
        df[["gene", "chr"]]
        .dropna()
        .drop_duplicates()
        .groupby("gene")["chr"]
        .agg(lambda x: x.iloc[0])
    )

    # gene-patient uniqueness
    gene_patient = df[["gene", "sample"]].drop_duplicates()

    recurrence = (
        gene_patient
        .groupby("gene")
        .agg(
            n_patients_altered=("sample", "nunique"),
            patients=("sample", lambda x: ",".join(sorted(set(x))))
        )
    )

    # state distribution
    state_counts = (
        df.groupby(["gene", "state"])
        .size()
        .rename("count")
        .reset_index()
    )

    state_summary = (
        state_counts
        .groupby("gene")
        .agg(state_counts=("state", lambda s: dict(zip(s, state_counts.loc[s.index, "count"]))))
    )

    # dominant state
    dominant_state = (
        state_counts.sort_values(["gene", "count"], ascending=[True, False])
        .drop_duplicates("gene")
        .set_index("gene")["state"]
        .rename("dominant_state")
    )

    # dominant direction
    direction_counts = (
        df.groupby(["gene", "direction"])
        .size()
        .rename("count")
        .reset_index()
    )

    dominant_direction = (
        direction_counts.sort_values(["gene", "count"], ascending=[True, False])
        .drop_duplicates("gene")
        .set_index("gene")["direction"]
        .rename("dominant_direction")
    )

    # combine everything
    result = (
        recurrence
        .join(gene_chr.rename("chr"), how="left")
        .join(dominant_state, how="left")
        .join(dominant_direction, how="left")
        .join(state_summary, how="left")
        .reset_index()
    )

    # fraction
    n_total_patients = df["sample"].nunique()
    result["fraction_patients_altered"] = (
        result["n_patients_altered"] / n_total_patients
    )

    # final sorting
    result = result.sort_values(
        ["n_patients_altered", "fraction_patients_altered", "gene"],
        ascending=[False, False, True]
    ).reset_index(drop=True)

    return result


import pandas as pd


def extract_top_genes_per_chromosome(
    df: pd.DataFrame,
    top_n: int = 3,
    chr_col: str = "chr",
    gene_col: str = "gene",
    gene_symbol_col: str = "gene_symbol",
    sort_by: list[str] | None = None,
    ascending: list[bool] | None = None,
    drop_na_chr: bool = True,
    exclude_chromosomes: list[str] | None = None,
) -> pd.DataFrame:
    """
    Extract the top N genes per chromosome from a recurrence/results dataframe.

    Parameters
    ----------
    df
        Input dataframe with at least chromosome and gene columns.
    top_n
        Number of top genes to keep per chromosome.
    chr_col
        Column containing chromosome labels.
    gene_col
        Column containing gene IDs (e.g. Ensembl IDs).
    gene_symbol_col
        Column containing gene symbols. If missing, output will still work.
    sort_by
        Columns to sort genes by within chromosome.
        Default:
            ["n_patients_altered", "fraction_patients_altered"]
    ascending
        Sort direction for `sort_by`.
        Default:
            [False, False]
    drop_na_chr
        Whether to drop rows with missing chromosome.
    exclude_chromosomes
        Optional list of chromosome labels to exclude.

    Returns
    -------
    DataFrame
        Top genes per chromosome with an added rank column.
    """
    out = df.copy()

    if drop_na_chr:
        out = out[out[chr_col].notna()].copy()

    if exclude_chromosomes is not None:
        out = out[~out[chr_col].isin(exclude_chromosomes)].copy()

    if sort_by is None:
        sort_by = ["n_patients_altered", "fraction_patients_altered"]

    if ascending is None:
        ascending = [False] * len(sort_by)

    missing = [c for c in [chr_col, gene_col] + sort_by if c not in out.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # sort globally first
    out = out.sort_values(
        [chr_col] + sort_by + [gene_col],
        ascending=[True] + ascending + [True]
    ).copy()

    # take top N per chromosome
    out["rank_within_chr"] = (
        out.groupby(chr_col).cumcount() + 1
    )

    out = out[out["rank_within_chr"] <= top_n].copy()

    # reorder columns nicely if available
    preferred_cols = [
        chr_col,
        "rank_within_chr",
        gene_symbol_col,
        gene_col,
        "n_patients_altered",
        "fraction_patients_altered",
        "dominant_direction",
        "dominant_state",
        "state_counts",
        "patients",
    ]
    cols = [c for c in preferred_cols if c in out.columns] + [
        c for c in out.columns if c not in preferred_cols
    ]

    return out[cols].reset_index(drop=True)


def collapse_top_genes_per_chromosome(
    top_df: pd.DataFrame,
    chr_col: str = "chr",
    gene_symbol_col: str = "gene_symbol",
    fallback_gene_col: str = "gene",
) -> pd.DataFrame:
    label_col = gene_symbol_col if gene_symbol_col in top_df.columns else fallback_gene_col

    out = (
        top_df.groupby(chr_col)[label_col]
        .apply(lambda x: ", ".join(x.astype(str)))
        .reset_index(name="top_genes")
    )
    return out


def print_genes_detailed(
    top_df,
    chr_col="chr",
    gene_symbol_col="gene_symbol",
):
    for chrom, sub in top_df.groupby(chr_col):
        print(f"\n{chrom}")
        for _, row in sub.sort_values("rank_within_chr").iterrows():
            gene = row.get(gene_symbol_col, row["gene"])
            direction = row.get("dominant_direction", "")
            n = row.get("n_patients_altered", "")
            print(f"  - {gene} ({direction}, n={n})")


def build_gene_mapping_from_adata(adata, ensembl_col='gene_ids', symbol_col=None):
    """
    Build mapping from Ensembl IDs to gene symbols from AnnData.
    """

    var = adata.var.copy()

    # determine columns
    if ensembl_col is None:
        # often var_names are Ensembl IDs
        var["gene_ids"] = var.index.astype(str)
    else:
        var["ensembl"] = var[ensembl_col].astype(str)

    var["symbol"] = var.index.astype(str)

    mapping = dict(zip(var["ensembl"], var["symbol"]))

    return mapping


SCRIPT_DIR = Path(__file__).resolve().parent
CNV_DIR    = SCRIPT_DIR.parent          # copy-number-variations/
REPO_ROOT  = CNV_DIR.parent
OUT_DIR    = CNV_DIR / "primary-cohort" / "outputs"

sample_to_file = {
    "combined_stroma_ref": str(
        CNV_DIR / "primary-cohort" / "cnv_inference/infercnv_runs/combined_stroma_ref"
        / "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat"
    ),
}

adata = sc.read_h5ad(REPO_ROOT / "quality_control" / "primary-cohort" / "adata.h5ad")

hmm_df = load_all_hmm_pred_genes(sample_to_file)
recurrent_genes = get_recurrent_altered_genes_with_meta(hmm_df, neutral_state=3)

gene_map = build_gene_mapping_from_adata(adata)
recurrent_genes["gene_symbol"] = recurrent_genes["gene"].map(gene_map)

recurrent_genes.to_csv(
    OUT_DIR / "recurrent_altered_genes.tsv",
    sep="\t",
    index=False
)
top_per_chr = extract_top_genes_per_chromosome(
    recurrent_genes,
    top_n=3
)

top_labels = collapse_top_genes_per_chromosome(top_per_chr.sort_values("n_patients_altered", ascending=False))
top_labels.to_csv(
    OUT_DIR / "top_genes_per_chromosome.tsv",
    sep="\t",
    index=False
)

print_genes_detailed(top_per_chr.sort_values(["n_patients_altered", "fraction_patients_altered"], ascending=[False, False]))
