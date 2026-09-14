"""
infercnv_subclusters.py
=======================
Extract inferCNV subclusters and their CNV fingerprints from HMM output.

The subcluster IDs here correspond 1-to-1 with the horizontal colour bands
in the inferCNV heatmap, so you can always trace a subcluster back to its
visual position in the plot.

Expected inferCNV output files (all in your --output-path directory)
---------------------------------------------------------------------
observation_groupings.txt
    Tab-separated, no header.
    Columns: barcode  subcluster_id
    Example:
        AAACAAGTATCTCCCA-1    Observation_1_1
        AAACACCAATAACTGC-1    Observation_1_1
        AAACAGAGCGACTCCT-1    Observation_1_2

HMM_CNV_predictions.HMMi6.leiden.hmm_post_marginals.pred_cnv_regions.dat
    Tab-separated, header row.
    Rows = genomic regions (chr:start-end), columns = cell barcodes.
    Values = integer CNV state (1=deep del, 2=del, 3=neutral, 4=dup, 5=amp,
             or 0/6 depending on inferCNV version).
    This is the *regions* file — one row per HMM-called region.

HMM_CNV_predictions.HMMi6.leiden.hmm_post_marginals.pred_cnv_genes.dat
    Same format but one row per gene.
    Use this when you want gene-level resolution.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns


# CNV state colour map

CNV_STATE_COLORS = {
    1: "#2166AC",
    2: "#92C5DE",
    3: "#F7F7F7",
    4: "#F4A582",
    5: "#D6604D",
    6: "#B2182B",
}

CNV_STATE_LABELS = {
    1: "Deletion",
    2: "Loss",
    3: "Neutral",
    4: "Gain",
    5: "Amplification",
    6: "High amplification",
}


# File loading

def load_observation_groupings(infercnv_dir: str | Path) -> pd.Series:
    path = Path(infercnv_dir) / "infercnv.17_HMM_predHMMi6.leiden.hmm_mode-subclusters.observation_groupings.txt"
    df = pd.read_csv(path, sep=r"\s+", header=0, quotechar='"')
    df.index = df.index.astype(str).str.strip('"')
    return df.iloc[:, 0].str.strip('"').rename("subcluster")


def load_hmm_predictions(
    infercnv_dir: str | Path,
    level: str = "regions",
) -> pd.DataFrame:
    infercnv_dir = Path(infercnv_dir)
    stem = f"pred_cnv_{level}.dat"
    matches = list(infercnv_dir.glob(f"*{stem}"))
    if not matches:
        raise FileNotFoundError(f"No file matching '*{stem}' found in {infercnv_dir}.")
    path = matches[0]
    print(f"Loading HMM predictions from: {path.name}")

    raw = pd.read_csv(path, sep="\t")

    cell_col   = "cell_group_name"
    region_col = "cnv_name" if level == "regions" else "gene_region_name"

    # For regions: each cell×region combination is already unique — no aggregation needed
    # For genes: multiple genes per region, take modal state
    if level == "regions":
        wide = raw.pivot(index=cell_col, columns=region_col, values="state")
    else:
        wide = (
            raw.groupby([cell_col, region_col])["state"]
            .agg(lambda x: x.mode().iloc[0])
            .unstack(level=region_col)
        )
    return wide


# extraction

def extract_subclusters(
    infercnv_dir: str | Path,
    level: str = "regions",
    reference_pattern: str = "Stroma",
) -> dict[str, pd.DataFrame]:
    infercnv_dir = Path(infercnv_dir)

    groupings = load_observation_groupings(infercnv_dir)   # barcode → subcluster
    cnv       = load_hmm_predictions(infercnv_dir, level=level)  # subcluster → regions

    # cnv index = full cell_group_name like "OVA09 - Tumor epithelium.OVA09 - Tumor epithelium_s1"
    # groupings values = subcluster like "OVA09 - Tumor epithelium_s1"
    # They share the suffix after the dot — build a mapping from subcluster → cnv row
    cnv_index_map = {idx.split(".")[-1]: idx for idx in cnv.index}
    groupings_mapped = groupings.map(cnv_index_map)  # barcode → full cnv row key

    missing_cnv      = groupings[groupings_mapped.isna()]
    missing_grouping = set(cnv.index) - set(groupings_mapped.dropna())

    if len(missing_grouping) > 0:
        print(f"[warn] {len(missing_grouping)} CNV rows have no matching grouping entry.")
    if len(missing_cnv) > 0:
        print(
            f"[info] {len(missing_cnv)} barcodes in groupings have no CNV data "
            f"— expected for reference spots (inferCNV skips HMM on reference)."
        )

    groupings_clean  = groupings[groupings_mapped.notna()].astype(str)
    groupings_mapped = groupings_mapped.dropna()

    # ── Barcode ↔ subcluster table ────────────────────────────────────────
    barcodes_df = pd.DataFrame({
        "barcode":      groupings_clean.index,
        "subcluster":   groupings_clean.values,
        "is_reference": groupings_clean.str.contains(reference_pattern, regex=True).values,
    })

    # ── CNV fingerprints: one row per subcluster (already aggregated) ─────
    # cnv rows are already per-subcluster; just rename index to short subcluster name
    fingerprints = cnv.copy()
    fingerprints.index = fingerprints.index.str.split(".").str[-1]
    fingerprints = fingerprints.fillna(3).astype(float)
    fingerprints = fingerprints.astype(float)

    # ── Per-subcluster summary ────────────────────────────────────────────
    n_barcodes   = groupings_clean.value_counts().rename("n_barcodes")
    is_reference = barcodes_df.groupby("subcluster")["is_reference"].first()

    region_cols  = list(fingerprints.columns)
    frac_neutral = (fingerprints[region_cols] == 3).mean(axis=1).rename("fraction_neutral")

    def _dominant_states(row):
        non_neutral = row[row != 3]
        if non_neutral.empty:
            return "all neutral"
        counts = non_neutral.value_counts()
        return ", ".join(
            f"state {int(s)} ({CNV_STATE_LABELS.get(int(s), '?')}): {n} regions"
            for s, n in counts.items()
        )

    dominant = fingerprints[region_cols].apply(_dominant_states, axis=1).rename("dominant_cnv_states")
    summary  = pd.concat([n_barcodes, is_reference, frac_neutral, dominant], axis=1)

    # fill NaN is_reference using the subcluster name
    summary["is_reference"] = summary["is_reference"].fillna(
        pd.Series(
            summary.index.str.contains(reference_pattern, regex=True),
            index=summary.index
        )
    )
    summary["n_barcodes"] = summary["n_barcodes"].fillna(0).astype(int)
    summary  = summary.sort_values(["is_reference", "fraction_neutral"], ascending=[True, False])

    return {
        "barcodes":           barcodes_df,
        "cnv_fingerprints":   fingerprints,
        "subcluster_summary": summary,
    }


# Add subclusters to AnnData

def add_subclusters_to_adata(adata, result: dict, obs_key: str = "infercnv_subcluster"):
    """
    Join the subcluster labels from extract_subclusters() onto an AnnData object.

    The join is done on obs_names (barcodes).  Barcodes not present in the
    inferCNV output are set to NaN.

    Parameters
    ----------
    adata
        AnnData whose obs_names are the same barcodes used in inferCNV.
    result
        The dict returned by extract_subclusters().
    obs_key
        Column name to add to adata.obs.
    """
    mapping = result["barcodes"].set_index("barcode")["subcluster"]
    adata.obs[obs_key] = adata.obs_names.map(mapping)

    n_mapped = adata.obs[obs_key].notna().sum()
    print(f"Mapped {n_mapped}/{adata.n_obs} barcodes to inferCNV subclusters.")
    print(adata.obs[obs_key].value_counts())
    return adata


# Visualisation

def load_region_coordinates(infercnv_dir: str | Path) -> pd.DataFrame:
    """
    Load genomic coordinates for each CNV region from the regions pred file.
    Returns a DataFrame indexed by region name with columns: chr, start, end, mid_mb.
    """
    infercnv_dir = Path(infercnv_dir)
    matches = list(infercnv_dir.glob("*pred_cnv_regions.dat"))
    if not matches:
        raise FileNotFoundError("No pred_cnv_regions.dat found.")

    raw = pd.read_csv(matches[0], sep="\t")
    coords = (
        raw[["cnv_name", "chr", "start", "end"]]
        .drop_duplicates("cnv_name")
        .set_index("cnv_name")
    )
    coords["mid_mb"] = (coords["start"] + coords["end"]) / 2 / 1e6

    def _sort_key(r):
        chrom_id = r.split("-")[0].replace("chr", "")
        chrom_num = int(chrom_id) if chrom_id.isdigit() else {"X": 23, "Y": 24, "M": 25}.get(chrom_id, 99)
        return (chrom_num, coords.loc[r, "start"])

    coords = coords.loc[sorted(coords.index, key=_sort_key)]
    return coords

def plot_cnv_fingerprints(
    result: dict,
    infercnv_dir: str | Path | None = None,
    show_reference: bool = False,
    figsize: tuple = (16, 6),
    save_path: str | None = None,
):
    fingerprints = result["cnv_fingerprints"].copy()
    summary      = result["subcluster_summary"]

    if not show_reference:
        obs_subs     = summary[~summary["is_reference"].astype(bool)].index
        fingerprints = fingerprints.loc[fingerprints.index.isin(obs_subs)]

    def _region_sort_key(r):
        chrom_id = r.split("-")[0].replace("chr", "")
        chrom_num = int(chrom_id) if chrom_id.isdigit() else {"X": 23, "Y": 24, "M": 25}.get(chrom_id, 99)
        region_num = int(r.split("region_")[-1]) if "region_" in r else 0
        return (chrom_num, region_num)

    region_cols = sorted(fingerprints.columns, key=_region_sort_key)
    data        = fingerprints[region_cols]

    # Load coordinates
    coords = None
    if infercnv_dir is not None:
        try:
            coords = load_region_coordinates(infercnv_dir)
            coords = coords.loc[coords.index.isin(region_cols)]
            coords = coords.loc[sorted(coords.index, key=_region_sort_key)]
            region_cols = list(coords.index)
            data = data[region_cols]
        except Exception as e:
            print(f"[warn] Could not load region coordinates: {e}")
            coords = None

    states = sorted(CNV_STATE_COLORS.keys())
    cmap   = mcolors.ListedColormap([CNV_STATE_COLORS[s] for s in states])
    bounds = [s - 0.5 for s in states] + [states[-1] + 0.5]
    norm   = mcolors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(figsize=figsize)
    n_rows  = len(data)

    if coords is not None:
        def _chrom_sort_key(c):
            cid = c.replace("chr", "")
            return int(cid) if cid.isdigit() else {"X": 23, "Y": 24, "M": 25}.get(cid, 99)

        unique_chroms = sorted(coords["chr"].unique(), key=_chrom_sort_key)

        chrom_max_pos = {
            chrom: coords[coords["chr"] == chrom]["end"].max()
            for chrom in unique_chroms
        }

        def _x_pos(chrom, pos_bp):
            """Map a bp position to a plot x coordinate."""
            chrom_idx  = unique_chroms.index(chrom)
            normalised = pos_bp / chrom_max_pos[chrom]   # 0–1 within chromosome
            return chrom_idx + normalised                  # offset by chromosome index

        # ── Draw neutral background for every chromosome ──────────────────
        for row_idx in range(n_rows):
            for chrom_idx, chrom in enumerate(unique_chroms):
                ax.add_patch(plt.Rectangle(
                    (chrom_idx, row_idx - 0.5), 1.0, 1.0,
                    color=CNV_STATE_COLORS[3], linewidth=0, zorder=0,
                ))

        # ── Draw CNV regions on top ───────────────────────────────────────
        for region in region_cols:
            row      = coords.loc[region]
            chrom    = row["chr"]
            x_start  = _x_pos(chrom, row["start"])
            x_end    = _x_pos(chrom, row["end"])
            width    = max(x_end - x_start, 0.005)

            for row_idx, subcluster in enumerate(data.index):
                state = int(data.loc[subcluster, region])
                if state == 3:
                    continue   # already drawn as background
                color = CNV_STATE_COLORS.get(state, CNV_STATE_COLORS[3])
                ax.add_patch(plt.Rectangle(
                    (x_start, row_idx - 0.5), width, 1.0,
                    color=color, linewidth=0, zorder=1,
                ))

        total_width = len(unique_chroms)
        ax.set_xlim(0, total_width)
        ax.set_ylim(-0.5, n_rows - 0.5)

        # Chromosome boundary lines and labels at the centre of each slot
        for chrom_idx, chrom in enumerate(unique_chroms):
            ax.axvline(chrom_idx, color="black", linewidth=0.5, alpha=0.5, zorder=2)
            ax.text(chrom_idx + 0.5, -0.6, chrom.replace("chr", ""),
                    ha="center", va="top", fontsize=7)
            
        ax.set_xticks([])
        ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.set_xlabel("Chromosome", fontsize=10, labelpad=15)

    else:
        # evenly spaced imshow
        im = ax.imshow(data.values, aspect="auto", cmap=cmap, norm=norm,
                       interpolation="nearest")
        chrom_of_region = [r.split("-")[0] for r in region_cols]
        chrom_ticks = {}
        for i, chrom in enumerate(chrom_of_region):
            if chrom not in chrom_ticks:
                chrom_ticks[chrom] = i
        ax.set_xticks(list(chrom_ticks.values()))
        ax.set_xticklabels(
            [c.replace("chr", "") for c in chrom_ticks.keys()],
            rotation=90, fontsize=7,
        )
        for tick_pos in chrom_ticks.values():
            ax.axvline(tick_pos - 0.5, color="black", linewidth=0.5, alpha=0.4)

        ax.set_xticks([])
        ax.xaxis.set_minor_locator(plt.NullLocator())
        ax.set_xlabel("Chromosome", fontsize=10, labelpad=15)

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(data.index, fontsize=8)
    ax.set_ylabel("Subcluster", fontsize=10)
    ax.set_title("InferCNV CNV fingerprints per subcluster", fontsize=12)
    ax.tick_params(axis="x", bottom=False, labelbottom=True)

    # Colourbar
    import matplotlib.cm as mcm
    sm   = mcm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, ax=ax, ticks=states, shrink=0.6, pad=0.02)
    cbar.ax.set_yticklabels([CNV_STATE_LABELS.get(s, str(s)) for s in states], fontsize=8)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()


def plot_subcluster_sizes(
    result: dict,
    show_reference: bool = False,
    figsize: tuple = (8, 5),
    save_path: str | None = None,
):
    """Bar chart of how many barcodes belong to each subcluster."""
    summary = result["subcluster_summary"].copy()
    if not show_reference:
        summary = summary[~summary["is_reference"]]

    fig, ax = plt.subplots(figsize=figsize)
    ax.barh(summary.index, summary["n_barcodes"], color="#4C72B0")
    ax.set_xlabel("Number of barcodes", fontsize=10)
    ax.set_ylabel("Subcluster", fontsize=10)
    ax.set_title("Spots per inferCNV subcluster", fontsize=12)
    ax.invert_yaxis()
    sns.despine()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    import anndata as ad

    from pathlib import Path

    REPO_ROOT = Path(__file__).resolve().parents[2]
    INFERCNV_DIR = str(
        REPO_ROOT
        / "copy-number-variations/external-cohort/cnv_inference"
        / "infercnv_runs/tumor-stroma-ova09"
    )
    ADATA_PATH = str(
        REPO_ROOT / "quality_control/external-cohort/pre-processing/adata.h5ad"
    )
    OUT_DIR = "outputs/infercnv_subcluster_outputs"

    # ── extract everything ────────────────────────────────────────
    result = extract_subclusters(
        infercnv_dir=INFERCNV_DIR,
        level="regions",
        reference_pattern="Ovarian stroma",
    )

    print("\n=== Subcluster summary ===")
    print(result["subcluster_summary"])

    # Which barcodes belong to which heatmap band?
    print("\n=== Barcode → subcluster (first 10) ===")
    print(result["barcodes"].head(10))

    # What is the CNV profile of each band?
    print("\n=== CNV fingerprints (first 5 regions) ===")
    print(result["cnv_fingerprints"].iloc[:, :5])

    # attach to AnnData
    adata = ad.read_h5ad(ADATA_PATH)
    adata = add_subclusters_to_adata(adata, result, obs_key="infercnv_subcluster")

    # visualise
    import os
    os.makedirs(OUT_DIR, exist_ok=True)

    plot_cnv_fingerprints(
        result,
        infercnv_dir=INFERCNV_DIR,
        show_reference=False,
        save_path=f"{OUT_DIR}/cnv_fingerprints.png",
    )
    plot_subcluster_sizes(
        result,
        show_reference=False,
        save_path=f"{OUT_DIR}/subcluster_sizes.png",
    )

    # save tables
    result["barcodes"].to_csv(f"{OUT_DIR}/barcode_to_subcluster.tsv", sep="\t", index=False)
    result["subcluster_summary"].to_csv(f"{OUT_DIR}/subcluster_summary.tsv", sep="\t")
    result["cnv_fingerprints"].to_csv(f"{OUT_DIR}/cnv_fingerprints.tsv", sep="\t")
