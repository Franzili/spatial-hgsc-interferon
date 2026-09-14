#!/usr/bin/env python3
"""
extract_cnv_spots.py
====================
Extract Visium spots (or single cells) with a copy-number alteration —
amplification OR deletion — for any gene of interest, from inferCNV HMM output.

Input files (standard inferCNV HMM output):
  • pred_cnv_genes.dat          – per-gene HMM state per subcluster
  • pred_cnv_regions.dat        – per-region HMM state per subcluster
  • observation_groupings.txt   – barcode → subcluster mapping

HMM state model (6-state inferCNV):
  1 = complete loss (0 copies)
  2 = loss of one copy (1 copy, heterozygous loss)
  3 = neutral (2 copies)  ← reference
  4 = gain (3 copies)
  5 = amplification (4+ copies)
  6 = high-level amplification

Gene query options (use any one):
  --gene  ENSG00000130303      Ensembl gene ID   (most reliable)
  --gene  BST2                 Gene symbol       (resolved via built-in table
                                                  or coordinate scan)
  --chr chr19 --start 17402939 --end 17405630    Genomic coordinates only
                                                  (GRCh38)

CNV type options:
  --cnv_type amplification     states ≥ --min_state (default 4)
  --cnv_type deletion          states ≤ --max_state (default 2)
  --cnv_type both              any non-neutral state

Usage examples:
  # BST2 amplification (default)
  python extract_cnv_spots.py --genes_file pred_cnv_genes.dat \\
      --regions_file pred_cnv_regions.dat --groupings obs_groupings.txt \\
      --gene BST2 --cnv_type amplification

  # TP53 deletion
  python extract_cnv_spots.py ... --gene TP53 --cnv_type deletion

  # CCNE1 amplification, strict threshold (only state 5+)
  python extract_cnv_spots.py ... --gene CCNE1 --cnv_type amplification --min_state 5

  # Any CNV at a locus by coordinates
  python extract_cnv_spots.py ... --chr chr17 --start 7661779 --end 7687550 \\
      --gene_label TP53 --cnv_type deletion

  # All CNV types for MYC
  python extract_cnv_spots.py ... --gene MYC --cnv_type both
"""

import argparse
import sys
import textwrap
from pathlib import Path

import pandas as pd

GENE_TABLE = {
    # Ovarian / HGSC drivers
    "TP53":    ("ENSG00000141510", "chr17", 7_661_779,  7_687_550),
    "BRCA1":   ("ENSG00000012048", "chr17", 43_044_295, 43_125_483),
    "BRCA2":   ("ENSG00000139618", "chr13", 32_315_086, 32_400_268),
    "CCNE1":   ("ENSG00000105173", "chr19", 28_817_812, 28_833_147),
    "MYC":     ("ENSG00000136997", "chr8",  127_735_434,127_742_951),
    "MECOM":   ("ENSG00000085djm", "chr3",  168_804_498,169_054_396),
    "PIK3CA":  ("ENSG00000121879", "chr3",  179_148_114,179_240_085),
    "PTEN":    ("ENSG00000171862", "chr10", 87_863_438, 87_971_930),
    "RB1":     ("ENSG00000139687", "chr13", 48_303_751, 48_481_890),
    "NF1":     ("ENSG00000196712", "chr17", 31_094_927, 31_377_677),
    "CDK12":   ("ENSG00000167258", "chr17", 37_673_802, 37_755_678),
    "KRAS":    ("ENSG00000133703", "chr12", 25_205_246, 25_250_936),
    "ERBB2":   ("ENSG00000141736", "chr17", 39_688_094, 39_728_658),
    "AKT1":    ("ENSG00000142208", "chr14", 104_769_349,104_795_748),
    "AKT2":    ("ENSG00000105221", "chr19", 40_230_317, 40_299_871),
    "AKT3":    ("ENSG00000117020", "chr1",  243_736_478,244_006_509),
    "NOTCH3":  ("ENSG00000074181", "chr19", 15_160_098, 15_200_995),
    "PALB2":   ("ENSG00000083093", "chr16", 23_603_160, 23_641_181),
    "CHEK2":   ("ENSG00000183765", "chr22", 28_687_742, 28_742_192),
    "RAD51C":  ("ENSG00000108384", "chr17", 56_769_984, 56_811_946),
    # Immune / BST2 context
    "BST2":    ("ENSG00000130303", "chr19", 17_402_939, 17_405_630),
    "CD274":   ("ENSG00000120217", "chr9",  5_450_503,  5_470_566),   # PD-L1
    "PDCD1":   ("ENSG00000188389", "chr2",  241_849_872,241_858_908), # PD-1
    "CD8A":    ("ENSG00000153563", "chr2",  86_945_987, 86_956_517),
    "TIGIT":   ("ENSG00000181492", "chr3",  114_239_168,114_256_744),
    # Other common cancer genes
    "EGFR":    ("ENSG00000146648", "chr7",  55_019_017, 55_211_628),
    "CDH1":    ("ENSG00000039068", "chr16", 68_737_291, 68_835_516),
    "CDKN2A":  ("ENSG00000147889", "chr9",  21_967_752, 21_995_324),
    "MDM2":    ("ENSG00000135679", "chr12", 69_202_998, 69_239_214),
    "VEGFA":   ("ENSG00000112715", "chr6",  43_737_901, 43_754_224),
    "FGFR1":   ("ENSG00000077782", "chr8",  38_411_132, 38_468_833),
    "FGFR2":   ("ENSG00000066468", "chr10", 121_478_397,121_598_238),
    "FGFR3":   ("ENSG00000068078", "chr4",  1_793_294,  1_808_872),
    "BRAF":    ("ENSG00000157764", "chr7",  140_719_327,140_924_764),
    "STK11":   ("ENSG00000118046", "chr19", 1_177_559,  1_228_434),
    "VHL":     ("ENSG00000134086", "chr3",  10_141_538, 10_153_959),
    "ARID1A":  ("ENSG00000117713", "chr1",  26_696_558, 26_883_524),
    "ATM":     ("ENSG00000149311", "chr11", 108_222_832,108_369_102),
    "MLH1":    ("ENSG00000076242", "chr3",  36_993_330, 37_050_846),
    "MSH2":    ("ENSG00000095002", "chr2",  47_402_913, 47_559_913),
    "DICER1":  ("ENSG00000100697", "chr14", 95_085_985, 95_164_534),
    "FOXL2":   ("ENSG00000183770", "chr3",  138_944_754,138_947_486),
}

HMM_STATE_LABELS = {
    1: "complete loss (0 copies)",
    2: "loss, 1 copy",
    3: "neutral (2 copies)",
    4: "gain (3 copies)",
    5: "amplification (4+ copies)",
    6: "high-level amplification",
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Extract spots/cells with a CNV (amp or del) for any gene from inferCNV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            HMM state guide (6-state model):
              1 = complete loss   2 = loss (1 copy)   3 = neutral (2 copies)
              4 = gain (3 copies) 5 = amplification   6 = high-level amp

            Gene symbol lookup uses a built-in table (~40 common genes).
            For any other gene, supply --gene as an Ensembl ID (ENSG...)
            or use --chr / --start / --end directly.
        """)
    )
    # Required files
    p.add_argument("--genes_file",   required=True,
                   help="pred_cnv_genes.dat from inferCNV")
    p.add_argument("--regions_file", required=True,
                   help="pred_cnv_regions.dat from inferCNV")
    p.add_argument("--groupings",    required=True,
                   help="observation_groupings.txt from inferCNV")

    # Gene specification (at least one required)
    gq = p.add_argument_group("Gene / locus specification (supply --gene OR coordinates)")
    gq.add_argument("--gene", default=None,
                    help="Ensembl ID (ENSG...) or gene symbol. "
                         "If a symbol, resolved via built-in table or coordinate scan.")
    gq.add_argument("--chr",   default=None, help="Chromosome, e.g. chr19")
    gq.add_argument("--start", type=int, default=None, help="Genomic start (GRCh38)")
    gq.add_argument("--end",   type=int, default=None, help="Genomic end (GRCh38)")
    gq.add_argument("--gene_label", default=None,
                    help="Label used in output filenames when querying by coordinates only")

    # CNV type
    p.add_argument("--cnv_type", default="amplification",
                   choices=["amplification", "deletion", "both"],
                   help="Type of CNV to extract. [default: amplification]")

    # Thresholds
    p.add_argument("--min_state", type=int, default=4, choices=range(1, 7),
                   help="For amplification: min HMM state (inclusive). [default: 4]")
    p.add_argument("--max_state", type=int, default=2, choices=range(1, 7),
                   help="For deletion: max HMM state (inclusive). [default: 2]")

    # Output
    p.add_argument("--out_dir", default="cnv_results",
                   help="Output directory. [default: cnv_results]")

    return p.parse_args()


# ---------------------------------------------------------------------------
# Gene resolution
# ---------------------------------------------------------------------------
def resolve_gene(gene_arg, chr_arg, start_arg, end_arg, gene_label_arg, genes_df):
    """
    Returns (ensg_id_or_None, chr, start, end, label) to use for querying.
    Resolution order:
      1. Explicit coordinates (--chr/--start/--end) override everything
      2. ENSG ID → look up coordinates from genes_df
      3. Known symbol → built-in GENE_TABLE
      4. Unknown symbol → scan genes_df for a gene at approximate coordinates
         (requires --chr/--start/--end to be provided alongside --gene)
    """
    ensg  = None
    chrom = chr_arg
    start = start_arg
    end   = end_arg
    label = gene_label_arg

    # Case 1: user gave full coordinates without a gene name
    if gene_arg is None:
        if chrom is None or start is None or end is None:
            sys.exit("ERROR: Supply either --gene or all of --chr/--start/--end.")
        label = label or f"{chrom}:{start}-{end}"
        return ensg, chrom, start, end, label

    # Case 2: Ensembl ID
    if gene_arg.upper().startswith("ENSG"):
        ensg = gene_arg.upper()
        # Try to get coordinates from genes_df
        rows = genes_df[genes_df["gene"].str.upper() == ensg]
        if len(rows):
            chrom = chrom or rows["chr"].iloc[0]
            start = start or int(rows["start"].min())
            end   = end   or int(rows["end"].max())
        elif chrom is None or start is None or end is None:
            sys.exit(f"ERROR: {ensg} not found in genes file and no coordinates supplied.\n"
                     f"       Add --chr/--start/--end manually.")
        label = label or ensg
        return ensg, chrom, start, end, label

    # Case 3: symbol in built-in table
    symbol_upper = gene_arg.upper()
    if symbol_upper in GENE_TABLE:
        entry = GENE_TABLE[symbol_upper]
        ensg  = entry[0]
        chrom = chrom or entry[1]
        start = start or entry[2]
        end   = end   or entry[3]
        label = label or gene_arg
        # Warn if the ENSG is not in the file (may have been filtered by inferCNV)
        if ensg not in genes_df["gene"].values:
            print(f"  NOTE: {gene_arg} ({ensg}) is not in the genes file "
                  f"(filtered out by inferCNV expression cutoff).")
            print(f"        Falling back to coordinate-based query: "
                  f"{chrom}:{start:,}-{end:,}")
            ensg = None   # will query by coords only
        return ensg, chrom, start, end, label

    # Case 4: unknown symbol – try coordinate scan if coords given
    if chrom and start and end:
        print(f"  NOTE: Symbol '{gene_arg}' not in built-in table. "
              f"Querying by supplied coordinates.")
        label = label or gene_arg
        return ensg, chrom, start, end, label

    sys.exit(
        f"ERROR: Gene '{gene_arg}' not recognised as an Ensembl ID and not in the\n"
        f"       built-in symbol table ({len(GENE_TABLE)} genes).\n"
        f"       Options:\n"
        f"         1. Use the Ensembl ID directly: --gene ENSG...\n"
        f"         2. Add --chr/--start/--end to query by coordinates\n"
        f"         3. Check the built-in table at the top of this script and add your gene\n"
        f"       Built-in symbols: {', '.join(sorted(GENE_TABLE.keys()))}"
    )


# ---------------------------------------------------------------------------
# State filter
# ---------------------------------------------------------------------------
def state_mask(series, cnv_type, min_state, max_state):
    """Return boolean mask for rows matching the requested CNV type."""
    if cnv_type == "amplification":
        return series >= min_state
    elif cnv_type == "deletion":
        return series <= max_state
    else:  # both
        return (series >= min_state) | (series <= max_state)


def cnv_type_label(cnv_type, min_state, max_state):
    if cnv_type == "amplification":
        return f"amplification (state ≥ {min_state}, {HMM_STATE_LABELS.get(min_state,'')})"
    elif cnv_type == "deletion":
        return f"deletion (state ≤ {max_state}, {HMM_STATE_LABELS.get(max_state,'')})"
    else:
        return (f"any CNV (amp: state ≥ {min_state}  |  del: state ≤ {max_state})")


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def load_genes(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()
    df["subcluster"] = df["cell_group_name"].str.split(".", n=1).str[-1]
    return df


def load_regions(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()
    df["subcluster"] = df["cell_group_name"].str.split(".", n=1).str[-1]
    return df


def load_groupings(path):
    df = pd.read_csv(path, sep=" ", quotechar='"', index_col=0)
    df.columns = df.columns.str.strip()
    df.index = df.index.str.strip()
    df = df.rename(columns={
        "Dendrogram Group": "subcluster",
        "Annotation Group": "annotation_group",
    })
    return df


# ---------------------------------------------------------------------------
# Core extraction
# ---------------------------------------------------------------------------
def find_gene_subclusters(genes_df, ensg, chrom, start, end,
                          cnv_type, min_state, max_state):
    """Find subclusters where the target gene has the requested CNV state."""
    if ensg:
        mask_gene = genes_df["gene"].str.upper() == ensg.upper()
    else:
        mask_gene = (
            (genes_df["chr"] == chrom) &
            (genes_df["start"] <= end) &
            (genes_df["end"] >= start)
        )
    hit = genes_df[mask_gene].copy()
    amp = hit[state_mask(hit["state"], cnv_type, min_state, max_state)].copy()
    return hit, amp   # return both so caller can see total rows for this gene


def find_region_subclusters(regions_df, chrom, start, end,
                             cnv_type, min_state, max_state):
    """Find subclusters where a CNV region overlaps the gene and has the CNV state."""
    coord_mask = (
        (regions_df["chr"] == chrom) &
        (regions_df["start"] <= end) &
        (regions_df["end"] >= start)
    )
    hit = regions_df[coord_mask].copy()
    amp = hit[state_mask(hit["state"], cnv_type, min_state, max_state)].copy()
    return hit, amp


def map_to_spots(subclusters, groupings):
    mask = groupings["subcluster"].isin(subclusters)
    df = groupings[mask].copy().reset_index().rename(columns={"index": "barcode"})
    return df[["barcode", "subcluster", "annotation_group"]]


# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load files
    # ------------------------------------------------------------------
    print("=" * 62)
    print("  InferCNV CNV spot extractor")
    print("=" * 62)
    print("[1/4] Loading inferCNV files ...")
    genes_df   = load_genes(args.genes_file)
    regions_df = load_regions(args.regions_file)
    groupings  = load_groupings(args.groupings)
    n_spots    = len(groupings)
    print(f"  Spots         : {n_spots:,}")
    print(f"  Subclusters   : {groupings['subcluster'].nunique()}")
    print(f"  Genes in file : {genes_df['gene'].nunique():,}  (unique Ensembl IDs)")

    # ------------------------------------------------------------------
    # Resolve gene
    # ------------------------------------------------------------------
    print("\n[2/4] Resolving gene query ...")
    ensg, chrom, start, end, label = resolve_gene(
        args.gene, args.chr, args.start, args.end, args.gene_label, genes_df
    )
    cnv_desc = cnv_type_label(args.cnv_type, args.min_state, args.max_state)

    print(f"  Gene label    : {label}")
    if ensg:
        print(f"  Ensembl ID    : {ensg}")
    print(f"  Locus         : {chrom}:{start:,}-{end:,}")
    print(f"  CNV type      : {cnv_desc}")

    # Safe filename prefix
    safe_label = label.replace(" ", "_").replace("/", "-").replace(":", "-")

    # ------------------------------------------------------------------
    # Find CNV subclusters
    # ------------------------------------------------------------------
    print("\n[3/4] Identifying CNV subclusters ...")

    gene_all_df, gene_cnv_df = find_gene_subclusters(
        genes_df, ensg, chrom, start, end,
        args.cnv_type, args.min_state, args.max_state
    )
    region_all_df, region_cnv_df = find_region_subclusters(
        regions_df, chrom, start, end,
        args.cnv_type, args.min_state, args.max_state
    )

    if len(gene_all_df) == 0 and len(region_all_df) == 0:
        print(f"\n  WARNING: No entries found for {label} in either file.")
        print(f"  The gene may have been excluded from the inferCNV run (expression cutoff),")
        print(f"  or the coordinates may not match the genome build used.")
        sys.exit(0)

    gene_cnv_subs   = set(gene_cnv_df["subcluster"].unique())
    region_cnv_subs = set(region_cnv_df["subcluster"].unique())
    any_subs  = gene_cnv_subs | region_cnv_subs
    both_subs = gene_cnv_subs & region_cnv_subs

    # Subclusters with data but not the queried CNV type (i.e. neutral or opposite)
    gene_neutral_subs = set(gene_all_df["subcluster"].unique()) - gene_cnv_subs

    print(f"  Gene-level   : {len(gene_cnv_subs):>4} subclusters  "
          f"(of {gene_all_df['subcluster'].nunique()} with any entry for this gene)")
    print(f"  Region-level : {len(region_cnv_subs):>4} subclusters  "
          f"(of {region_all_df['subcluster'].nunique()} with overlapping region)")
    print(f"  Union        : {len(any_subs):>4} subclusters")
    print(f"  Intersection : {len(both_subs):>4} subclusters")

    # ------------------------------------------------------------------
    # Map to spots
    # ------------------------------------------------------------------
    spots_gene   = map_to_spots(gene_cnv_subs,  groupings)
    spots_region = map_to_spots(region_cnv_subs, groupings)
    spots_any    = map_to_spots(any_subs,        groupings)
    spots_both   = map_to_spots(both_subs,       groupings)

    # Warn about subclusters present in CNV files but absent from groupings
    missing = any_subs - set(groupings["subcluster"].unique())
    if missing:
        print(f"\n  NOTE: {len(missing)} CNV-positive subclusters absent from groupings")
        print(f"  (reference 'Other' cells — excluded from spot output):")
        for s in sorted(missing)[:5]:
            print(f"    {s}")
        if len(missing) > 5:
            print(f"    ... and {len(missing)-5} more")

    print(f"\n  Spots — gene-level evidence   : {len(spots_gene):,}  "
          f"({100*len(spots_gene)/n_spots:.1f}%)")
    print(f"  Spots — region-level evidence : {len(spots_region):,}  "
          f"({100*len(spots_region)/n_spots:.1f}%)")
    print(f"  Spots — any evidence          : {len(spots_any):,}  "
          f"({100*len(spots_any)/n_spots:.1f}%)")
    print(f"  Spots — both evidence         : {len(spots_both):,}  "
          f"({100*len(spots_both)/n_spots:.1f}%)")

    # ------------------------------------------------------------------
    # Build annotated spot table
    # ------------------------------------------------------------------
    # Per-subcluster state summary
    def state_summary(cnv_df, prefix):
        if len(cnv_df) == 0:
            return pd.DataFrame(columns=["subcluster",
                                          f"{prefix}_state_max",
                                          f"{prefix}_state_min"])
        return (cnv_df.groupby("subcluster")["state"]
                .agg(**{f"{prefix}_state_max": "max",
                        f"{prefix}_state_min": "min"})
                .reset_index())

    gene_state_df   = state_summary(gene_cnv_df,   "gene")
    region_state_df = state_summary(region_cnv_df, "region")

    spot_table = (
        spots_any
        .assign(
            gene_cnv   = spots_any["subcluster"].isin(gene_cnv_subs),
            region_cnv = spots_any["subcluster"].isin(region_cnv_subs),
        )
        .merge(gene_state_df,   on="subcluster", how="left")
        .merge(region_state_df, on="subcluster", how="left")
    )

    # Subcluster-level summary (all subclusters in groupings)
    sub_summary = (
        groupings.groupby("subcluster")
        .size().reset_index(name="n_spots")
        .merge(gene_state_df,   on="subcluster", how="left")
        .merge(region_state_df, on="subcluster", how="left")
    )
    sub_summary["gene_cnv"]   = sub_summary["subcluster"].isin(gene_cnv_subs)
    sub_summary["region_cnv"] = sub_summary["subcluster"].isin(region_cnv_subs)
    sub_summary["any_cnv"]    = sub_summary["subcluster"].isin(any_subs)

    # State distribution among CNV+ subclusters
    if len(gene_cnv_df):
        state_dist = (
            gene_cnv_df.groupby("state")
            .agg(n_subclusters=("subcluster", "nunique"))
            .reset_index()
            .assign(state_label=lambda d: d["state"].map(HMM_STATE_LABELS))
        )
    else:
        state_dist = pd.DataFrame()

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    print("\n[4/4] Writing outputs ...")

    for ev_label, df in [
        ("gene_evidence",   spots_gene),
        ("region_evidence", spots_region),
        ("any_evidence",    spots_any),
        ("both_evidence",   spots_both),
    ]:
        if len(df):
            path = out_dir / f"{safe_label}_{args.cnv_type}_{ev_label}_barcodes.txt"
            df["barcode"].to_csv(path, index=False, header=False)
            print(f"  Written: {path}  ({len(df):,} barcodes)")

    spot_path = out_dir / f"{safe_label}_{args.cnv_type}_spots_annotated.tsv"
    spot_table.to_csv(spot_path, sep="\t", index=False)
    print(f"  Written: {spot_path}")

    sub_path = out_dir / f"{safe_label}_{args.cnv_type}_subcluster_summary.tsv"
    sub_summary.to_csv(sub_path, sep="\t", index=False)
    print(f"  Written: {sub_path}")

    if len(state_dist):
        sd_path = out_dir / f"{safe_label}_{args.cnv_type}_state_distribution.tsv"
        state_dist.to_csv(sd_path, sep="\t", index=False)
        print(f"  Written: {sd_path}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 62)
    print("  SUMMARY")
    print("=" * 62)
    print(f"  Gene           : {label}")
    print(f"  Locus          : {chrom}:{start:,}-{end:,}")
    print(f"  CNV type       : {cnv_desc}")
    print(f"  Total spots    : {n_spots:,}")
    print(f"  Affected spots : {len(spots_any):,}  ({100*len(spots_any)/n_spots:.1f}%)")

    if len(state_dist):
        print(f"\n  HMM state breakdown (gene-level):")
        for _, row in state_dist.iterrows():
            print(f"    state {int(row['state'])}  {row['state_label']:<35}"
                  f"  {int(row['n_subclusters'])} subclusters")

    recommended = out_dir / f"{safe_label}_{args.cnv_type}_gene_evidence_barcodes.txt"
    if not recommended.exists():
        recommended = out_dir / f"{safe_label}_{args.cnv_type}_any_evidence_barcodes.txt"


if __name__ == "__main__":
    main()


# Example usage (BST2 amplification):
# Run from the copy-number-variations/ directory:
#
# python scripts/extract_spots_with_cnv_gene.py \
#     --genes_file "primary-cohort/cnv_inference/infercnv_runs/combined_stroma_ref/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_genes.dat" \
#     --regions_file "primary-cohort/cnv_inference/infercnv_runs/combined_stroma_ref/HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat" \
#     --groupings "primary-cohort/cnv_inference/infercnv_runs/combined_stroma_ref/infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observation_groupings.txt" \
#     --out_dir "primary-cohort/outputs/BST2_results" \
#     --cnv_type amplification \
#     --gene BST2
