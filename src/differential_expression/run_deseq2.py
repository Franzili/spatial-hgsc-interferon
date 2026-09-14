"""
run_deseq2.py
─────────────
Pseudobulk DESeq2 differential expression analysis for Visium FFPE data.
Designed to be run on Puhti (or any SLURM HPC) via run_deseq2.sh.

Inputs:
  ANNDATA_PATH   : path to the .h5ad AnnData file
  OUT_DEGS       : output directory for results CSVs and pkl files

Usage (via SLURM):
  sbatch run_deseq2.sh
"""

import os
import sys
import pickle
import argparse

import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference


# ── CONFIG ────────────────────────────────────────────────────────────────────

DEFAULT_ANNDATA = os.environ.get("ANNDATA_PATH", "quality_control/adata.h5ad")
DEFAULT_OUTDIR  = os.environ.get("OUT_DEGS",     "deseq2_results")

# DESeq2 significance thresholds
PADJ_THRESH = 0.05
LFC_THRESH  = 1.0

# Gene pre-filtering: keep genes with >= MIN_COUNT in >= MIN_SAMPLES samples
MIN_COUNT   = 10
MIN_SAMPLES = 2

# Number of CPUs for pydeseq2 inference
N_CPUS = int(os.environ.get("SLURM_CPUS_PER_TASK", 24))


# ── ARGUMENT PARSING ─────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(description="Pseudobulk DESeq2 analysis")
    parser.add_argument("--anndata", default=DEFAULT_ANNDATA,
                        help="Path to input .h5ad AnnData file")
    parser.add_argument("--outdir",  default=DEFAULT_OUTDIR,
                        help="Output directory for results")
    parser.add_argument("--min-count",   type=int, default=MIN_COUNT)
    parser.add_argument("--min-samples", type=int, default=MIN_SAMPLES)
    parser.add_argument("--padj",  type=float, default=PADJ_THRESH)
    parser.add_argument("--lfc",   type=float, default=LFC_THRESH)
    parser.add_argument("--ncpus", type=int,   default=N_CPUS)
    return parser.parse_args()


# ── HELPER FUNCTIONS ─────────────────────────────────────────────────────────

def filter_genes_for_pseudobulk(counts_df, min_count, min_samples):
    """
    Remove lowly expressed genes before passing to DESeq2.

    Keeps only genes with at least min_count counts in at least
    min_samples pseudobulk samples. After aggregation across many spots,
    nearly every gene passes a simple total-count filter; this stricter
    criterion removes near-zero genes that inflate runtime and produce
    unreliable dispersion estimates.
    """
    numeric_cols = counts_df.select_dtypes(include=[np.number]).columns
    counts_df    = counts_df[numeric_cols]

    # Filter minimum expression
    expr_mask = (counts_df >= min_count).sum(axis=0) >= min_samples
    counts_df = counts_df.loc[:, expr_mask]
    print(f"Genes after expression filter : {counts_df.shape[1]:,}")

    # Filter remove near-constant genes (variance ~ 0 across samples)
    gene_means = counts_df.mean(axis=0)
    gene_stds  = counts_df.std(axis=0)

    cv         = gene_stds / gene_means.replace(0, np.nan)
    var_mask   = cv > 0.1
    counts_df  = counts_df.loc[:, var_mask]
    print(f"Genes after variance filter   : {counts_df.shape[1]:,}")

    return counts_df


def fit_deseq_model(counts_df, metadata, design, name, outdir, inference):
    """
    Fit a DESeq2 model and immediately pickle the fitted dds object.

    Saving immediately means the model is available for reuse even if a
    later contrast or the rest of the script fails mid-run.
    """
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design=design,
        refit_cooks=True,
        inference=inference,
    )
    dds.deseq2()

    name_clean = name.replace(" ", "_").replace(":", "_")
    pkl_path   = os.path.join(outdir, f"{name_clean}_dds.pkl")
    with open(pkl_path, "wb") as f:
        pickle.dump(dds, f)
    print(f"  Saved dds: {pkl_path}")

    return dds


def extract_and_save_contrast(dds, contrast, name, outdir, inference,
                               padj=PADJ_THRESH, lfc=LFC_THRESH):
    """
    Extract one contrast from a fitted model, apply LFC shrinkage, and
    save full (unshrunk), shrunk and significant results to CSV immediately.

    Significance calls use unshrunk padj; effect sizes for plotting use
    shrunk log2FoldChange. Both are saved separately so downstream code
    can pick the appropriate column for each purpose.
    """
    stat_res = DeseqStats(dds, contrast=contrast, inference=inference)
    stat_res.summary()

    # Unshrunk results — padj is computed here and must not be overwritten
    results_unshrunk = stat_res.results_df.copy()

    name_clean = name.replace(" ", "_").replace(":", "_")

    # Full unshrunk — use padj column from this file for significance
    full_path = os.path.join(outdir, f"{name_clean}_full.csv")
    results_unshrunk.to_csv(full_path)

    available_coeffs = dds.varm["LFC"].columns.tolist()
    factor, level_num, level_ref = contrast

    matching = [
        c for c in available_coeffs
        if factor in c and level_num in c
    ]

    # LFC shrinkage
    results_shrunk = results_unshrunk  # fallback if shrinkage fails
    if len(matching) == 1:
        coeff = matching[0]
        try:
            stat_res.lfc_shrink(coeff=coeff)
            results_shrunk = stat_res.results_df.copy()
            print(f"  LFC shrinkage applied (coeff: {coeff})")
        except Exception as e:
            print(f"  LFC shrinkage failed ({e}) — using unshrunk LFCs")
    elif len(matching) == 0:
        # This happens when level_num is the reference level — DESeq2 only
        # creates a coefficient for the non-reference levels, so the contrast
        # is estimable but there is no direct coefficient to shrink against.
        # In this case unshrunk LFCs are used for plotting.
        print(f"  No coefficient found for contrast {contrast} "
            f"('{level_num}' may be the reference level) — "
            f"skipping shrinkage, using unshrunk LFCs")
    else:
        print(f"  Ambiguous coefficient match for contrast {contrast}: "
            f"{matching} — skipping shrinkage, using unshrunk LFCs")

    # Shrunk — use log2FoldChange column from this file for plotting
    shrunk_path = os.path.join(outdir, f"{name_clean}_shrunk.csv")
    results_shrunk.to_csv(shrunk_path)

    # Significant subset — padj from unshrunk, LFC from shrunk
    sig_df   = get_sig(
        results_unshrunk.assign(log2FoldChange=results_shrunk["log2FoldChange"]),
        padj=padj, lfc=lfc,
    )
    sig_path = os.path.join(outdir, f"{name_clean}_significant.csv")
    sig_df.to_csv(sig_path)

    print(f"  Significant DEGs : {len(sig_df)}")
    print(f"  Saved: {full_path}")
    print(f"  Saved: {shrunk_path}")
    print(f"  Saved: {sig_path}")

    # Return unshrunk for summary table (padj needed there);
    # caller can load the shrunk CSV when plotting is needed
    return results_unshrunk, results_shrunk


def get_sig(results_df, padj=PADJ_THRESH, lfc=LFC_THRESH):
    """Filter DESeq2 results to significant DEGs."""
    return results_df[
        (results_df["padj"] < padj) &
        (results_df["log2FoldChange"].abs() > lfc)
    ].sort_values("padj")


def collect_summary(results_df, name, summary_rows,
                    padj=PADJ_THRESH, lfc=LFC_THRESH):
    """Append one row to the running summary list for a given contrast."""
    sig = get_sig(results_df, padj=padj, lfc=lfc)
    summary_rows.append({
        "comparison":      name,
        "n_tested":        results_df["padj"].notna().sum(),
        "n_significant":   len(sig),
        "n_upregulated":   (sig["log2FoldChange"] > 0).sum(),
        "n_downregulated": (sig["log2FoldChange"] < 0).sum(),
        "padj_threshold":  padj,
        "lfc_threshold":   lfc,
    })


# ── PSEUDOBULK PREPARATION ────────────────────────────────────────────────────

def build_pseudobulk(adata):
    import scipy.sparse

    counts = adata.layers["counts"]
    if scipy.sparse.issparse(counts):
        counts = counts.toarray()

    # Build counts DataFrame from gene columns only
    gene_names = adata.var_names.tolist()
    counts_df  = pd.DataFrame(
        counts.astype("int32"),
        index=adata.obs_names,
        columns=gene_names,
    )

    # Attach grouping columns separately
    counts_df["condition"] = adata.obs["PFI_short_long"].values
    counts_df["histology"] = adata.obs["histology"].values
    counts_df["class"]     = adata.obs["class_overall"].values  # use class_overall
    counts_df["patient"]   = adata.obs["patient"].values

    # Verify no gene name collides with the groupby keys before aggregating
    groupby_keys = ["condition", "histology", "class", "patient"]
    collisions   = [k for k in groupby_keys if k in gene_names]
    if collisions:
        raise ValueError(
            f"Gene name(s) {collisions} collide with groupby column names. "
            f"Rename the groupby keys before proceeding."
        )

    pseudobulk_df = (
        counts_df
        .groupby(groupby_keys, observed=True)
        .sum()
    )

    metadata = pseudobulk_df.index.to_frame(index=False)
    metadata.columns = groupby_keys
    for col in metadata.columns:
        metadata[col] = metadata[col].astype("category")

    sample_ids = [
        f"{c}_{h}_{s}_{p}"
        for c, h, s, p in zip(
            metadata["condition"], metadata["histology"],
            metadata["class"],     metadata["patient"],
        )
    ]
    pseudobulk_df.index = sample_ids
    metadata.index      = sample_ids

    pseudobulk_df = pseudobulk_df.drop(columns=groupby_keys, errors="ignore")

    non_numeric = pseudobulk_df.select_dtypes(exclude=[np.number]).columns.tolist()
    if non_numeric:
        raise ValueError(
            f"Non-numeric columns remain after pseudobulking: {non_numeric}. "
            f"This likely means a gene name collides with an obs column name."
        )

    print(f"  Pseudobulk matrix : {pseudobulk_df.shape[0]} samples × "
          f"{pseudobulk_df.shape[1]} genes")
    
    print()
    print(f"Total pseudobulk samples: {len(metadata)}")
    print(f"Patients: {metadata['patient'].nunique()}")

    print("\nSamples per condition:")
    print(metadata["condition"].value_counts())

    print("\nSamples per histology:")
    print(metadata["histology"].value_counts())

    print("\nSamples per class:")
    print(metadata["class"].value_counts())

    print("\nSamples per patient:")
    print(metadata["patient"].value_counts().sort_index())

    print("\nCondition x class cross-tabulation:")
    print(pd.crosstab(metadata["condition"], metadata["class"]))

    print("\nCondition x histology cross-tabulation:")
    print(pd.crosstab(metadata["condition"], metadata["histology"]))

    print("\nPatient x condition cross-tabulation:")
    print(pd.crosstab(metadata["patient"], metadata["condition"]))

    print("\nSamples per patient x condition x histology:")
    print(metadata.groupby(["patient", "condition", "histology"], observed=True)
        .size()
        .unstack(fill_value=0))
    print()

    return pseudobulk_df, metadata


# ── MAIN ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Output directory : {args.outdir}")
    print(f"AnnData path     : {args.anndata}")
    print(f"CPUs             : {args.ncpus}")
    print(f"Gene filter      : min_count={args.min_count}, "
          f"min_samples={args.min_samples}")
    print(f"Significance     : padj<{args.padj}, |LFC|>{args.lfc}\n")

    # ── Load data ─────────────────────────────────────────────────────────────
    import scanpy as sc
    print("Loading AnnData...")
    adata = sc.read_h5ad(args.anndata)
    print(f"  {adata.n_obs} spots, {adata.n_vars} genes")

    adata.obs['class'] = adata.obs['class'].astype(str)
    adata.obs['class_overall'] = adata.obs['class'].replace({'OvaryR': 'Ovary', 'OvaryL': 'Ovary'})
    adata.obs['PFI'] = adata.obs.PFI.astype(str)
    adata.obs['PFI_short_long'] = adata.obs['PFI'].replace({'short': 'short', 'medium': 'long', 'long': 'long'})
    adata.obs['patient'] = adata.obs['patient'].astype(str)
    adata.obs['anno'] = adata.obs['class'] + '_' + adata.obs['patient'] + '_PFI-' + adata.obs['PFI']
    adata = adata[~adata.obs['patient'].isin(['Marker'])].copy()

    # ── Build pseudobulk ─────────────────────────────────────────────────────
    print("\nBuilding pseudobulk matrix...")
    counts_df, metadata = build_pseudobulk(adata)
    print(f"  {counts_df.shape[0]} pseudobulk samples, "
          f"{counts_df.shape[1]} genes before filtering")
    print(f"\nMetadata preview:\n{metadata.head()}\n")
    print(f"Sample counts per group:\n"
          f"{metadata.groupby(['condition','histology','class']).size()}\n")

    # ── Inference engine ──────────────────────────────────────────────────────
    inference = DefaultInference(n_cpus=args.ncpus)

    # ── Pre-filter genes once on the full matrix ──────────────────────────────
    counts_df = filter_genes_for_pseudobulk(
        counts_df, args.min_count, args.min_samples
    )

    summary_rows = []

    # ── Condition: short PFI vs long PFI ───────────────────────────────────
    print("\n── Condition contrast ──")
    dds_cond = fit_deseq_model(
        counts_df.loc[metadata.index],
        metadata,
        design="~ histology + condition",
        name="model_condition",
        outdir=args.outdir,
        inference=inference,
    )
    res_condition_unshrunk, _ = extract_and_save_contrast(
        dds_cond,
        contrast=["condition", "short", "long"],
        name="condition_short_vs_long",
        outdir=args.outdir,
        inference=inference,
        padj=args.padj,
        lfc=args.lfc,
    )
    collect_summary(res_condition_unshrunk, "condition_short_vs_long",
                    summary_rows, padj=args.padj, lfc=args.lfc)

    # ── Histology: pairwise comparisons between tissue domains ─────────────
    print("\n── Histology contrasts ──")

    hist_counts = metadata["histology"].value_counts()
    keep_hist   = hist_counts[hist_counts >= 2].index
    meta_hist   = metadata[metadata["histology"].isin(keep_hist)]
    counts_hist = counts_df.loc[meta_hist.index]
    hist_levels = (
        meta_hist["histology"]
        .cat.remove_unused_categories()
        .cat.categories
        .tolist()
    )
    print(f"Histology levels with sufficient replication: {hist_levels}")

    # Fit once, extract all pairwise contrasts from the same model
    dds_hist = fit_deseq_model(
        counts_hist, meta_hist,
        design="~ condition + histology",
        name="model_histology",
        outdir=args.outdir,
        inference=inference,
    )
    for i, h1 in enumerate(hist_levels):
        for h2 in hist_levels[i + 1:]:
            print(f"  Histology: {h1} vs {h2}")
            key = f"{h1}_vs_{h2}"
            res_unshrunk, _ = extract_and_save_contrast(
                dds_hist,
                contrast=["histology", h1, h2],
                name=f"histology_{key}",
                outdir=args.outdir,
                inference=inference,
                padj=args.padj,
                lfc=args.lfc,
            )
            collect_summary(res_unshrunk, f"histology_{key}",
                            summary_rows, padj=args.padj, lfc=args.lfc)
            
    # ── Condition effect stratified by histology ──────────────────────────────────
    print("\n── Condition within each histology ──")

    for hist in metadata["histology"].cat.categories:
        meta_sub        = metadata[metadata["histology"] == hist].copy()
        cond_counts_sub = meta_sub["condition"].value_counts()

        print(f"\n  Histology: {hist}")
        print(f"  Condition counts: {cond_counts_sub.to_dict()}")

        if (cond_counts_sub >= 2).all():
            counts_sub = counts_df.loc[meta_sub.index]
            
            design = "~ condition"
            print(f"  Design: {design}")

            dds_sub = fit_deseq_model(
                counts_sub,
                meta_sub,
                design=design,
                name=f"model_condition_in_histology_{hist}",
                outdir=args.outdir,
                inference=inference,
            )
            res_unshrunk, _ = extract_and_save_contrast(
                dds_sub,
                contrast=["condition", "short", "long"],
                name=f"condition_in_histology_{hist}",
                outdir=args.outdir,
                inference=inference,
                padj=args.padj,
                lfc=args.lfc,
            )
            collect_summary(
                res_unshrunk,
                f"condition_in_histology_{hist}",
                summary_rows,
                padj=args.padj,
                lfc=args.lfc,
            )
        else:
            print(f"  Skipping: need >= 2 samples per condition "
                f"(need >= 2 per level)")

    # ── Class: pairwise comparisons between sampling locations ─────────────
    print("\n── Class contrasts ──")

    class_counts_n = metadata["class"].value_counts()
    keep_class     = class_counts_n[class_counts_n >= 2].index
    meta_class     = metadata[metadata["class"].isin(keep_class)]
    counts_class   = counts_df.loc[meta_class.index]
    class_levels   = (
        meta_class["class"]
        .cat.remove_unused_categories()
        .cat.categories
        .tolist()
    )
    print(f"Class levels with sufficient replication: {class_levels}")

    dds_class = fit_deseq_model(
        counts_class, meta_class,
        design="~ histology + class",
        name="model_class",
        outdir=args.outdir,
        inference=inference,
    )
    for i, c1 in enumerate(class_levels):
        for c2 in class_levels[i + 1:]:
            print(f"  Class: {c1} vs {c2}")
            key = f"{c1}_vs_{c2}"
            res_unshrunk, _ = extract_and_save_contrast(
                dds_class,
                contrast=["class", c1, c2],
                name=f"class_{key}",
                outdir=args.outdir,
                inference=inference,
                padj=args.padj,
                lfc=args.lfc,
            )
            collect_summary(res_unshrunk, f"class_{key}",
                            summary_rows, padj=args.padj, lfc=args.lfc)

    # ── 4. Condition effect stratified by class ───────────────────────────────
    # Subset to each class and test condition within it, rather than fitting
    # a class:condition interaction term which requires dense replication
    print("\n── Condition within each class ──")

    for cls in metadata["class"].cat.categories:
        meta_sub       = metadata[metadata["class"] == cls].copy()
        cond_counts_sub = meta_sub["condition"].value_counts()

        if (cond_counts_sub >= 2).all():
            counts_sub = counts_df.loc[meta_sub.index]
            print(f"  Fitting condition model within class: {cls}")
            dds_sub = fit_deseq_model(
                counts_sub, meta_sub,
                design="~ histology + condition",
                name=f"model_condition_in_class_{cls}",
                outdir=args.outdir,
                inference=inference,
            )
            res_unshrunk, _ = extract_and_save_contrast(
                dds_sub,
                contrast=["condition", "short", "long"],
                name=f"condition_in_class_{cls}",
                outdir=args.outdir,
                inference=inference,
                padj=args.padj,
                lfc=args.lfc,
            )
            collect_summary(res_unshrunk, f"condition_in_class_{cls}",
                            summary_rows, padj=args.padj, lfc=args.lfc)
        else:
            print(f"  Skipping class '{cls}': condition counts = "
                  f"{cond_counts_sub.to_dict()} (need >= 2 per level)")

    # ── 5. Condition effect stratified by class, within Tumor Epithelium ──────
    # Same idea as #4, but subset to a single histology first so the site
    # comparison isn't confounded by differing tumor/stroma proportions
    # between sites. Histology is constant within each subset, so it's
    # dropped from the design (no variation left to control for).
    print("\n── Condition within each class, Tumor Epithelium only ──")

    TARGET_HISTOLOGY = "Tumor Epithelium"
    meta_tumor = metadata[metadata["histology"] == TARGET_HISTOLOGY]

    for cls in metadata["class"].cat.categories:
        meta_sub        = meta_tumor[meta_tumor["class"] == cls].copy()
        cond_counts_sub = meta_sub["condition"].value_counts()

        if (cond_counts_sub >= 2).all():
            counts_sub = counts_df.loc[meta_sub.index]
            print(f"  Fitting condition model within class: {cls} "
                  f"({TARGET_HISTOLOGY} only)")
            dds_sub = fit_deseq_model(
                counts_sub, meta_sub,
                design="~ condition",
                name=f"model_condition_in_class_{cls}_tumor_epithelium",
                outdir=args.outdir,
                inference=inference,
            )
            res_unshrunk, _ = extract_and_save_contrast(
                dds_sub,
                contrast=["condition", "short", "long"],
                name=f"condition_in_class_{cls}_tumor_epithelium",
                outdir=args.outdir,
                inference=inference,
                padj=args.padj,
                lfc=args.lfc,
            )
            collect_summary(res_unshrunk, f"condition_in_class_{cls}_tumor_epithelium",
                            summary_rows, padj=args.padj, lfc=args.lfc)
        else:
            print(f"  Skipping class '{cls}' ({TARGET_HISTOLOGY} only): "
                  f"condition counts = {cond_counts_sub.to_dict()} "
                  f"(need >= 2 per level)")


    # ── Summary table ─────────────────────────────────────────────────────────
    print("\n── Writing summary ──")
    summary_df   = pd.DataFrame(summary_rows)
    summary_path = os.path.join(args.outdir, "summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"Summary saved to {summary_path}")
    print(summary_df.to_string(index=False))
    print("\nDone.")

    counts_df.to_csv(os.path.join(args.outdir, "pseudobulk_counts.csv"))
    metadata.to_csv(os.path.join(args.outdir, "pseudobulk_metadata.csv"))
    print("Saved pseudobulk counts and metadata")


if __name__ == "__main__":
    main()
