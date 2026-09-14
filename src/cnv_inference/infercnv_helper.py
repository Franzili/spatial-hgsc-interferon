from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class SampleRun:
    sample: str
    counts_h5: Path
    annotations_tsv: Path
    out_dir: Path
    ref_groups: list[str]
    cutoff: float
    hmm: bool


def discover_sample_inputs(input_dir: Path) -> dict[str, dict[str, Path]]:
    """
    Expect files named like:
      <sample>_counts.h5
      <sample>_infercnv_annotations.tsv
    """
    found: dict[str, dict[str, Path]] = {}

    for h5 in input_dir.glob("*_counts.h5"):
        sample = h5.name.removesuffix("_counts.h5")
        found.setdefault(sample, {})["counts_h5"] = h5.resolve()

    for ann in input_dir.glob("*_infercnv_annotations.tsv"):
        sample = ann.name.removesuffix("_infercnv_annotations.tsv")
        found.setdefault(sample, {})["annotations_tsv"] = ann.resolve()

    return found


def load_config(config_path: Path) -> dict[str, Any]:
    with open(config_path) as f:
        return json.load(f)


def parse_bool(value: Any, field_name: str) -> bool:
    """
    Parse booleans safely from bool or common string forms.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        val = value.strip().lower()
        if val in {"true", "t", "1", "yes", "y"}:
            return True
        if val in {"false", "f", "0", "no", "n"}:
            return False
    raise ValueError(f"'{field_name}' must be a boolean or boolean-like string, got: {value!r}")


def build_runs(config: dict[str, Any]) -> list[SampleRun]:
    input_dir = Path(config["input_dir"]).resolve()
    output_root = Path(config["output_root"]).resolve()
    discovered = discover_sample_inputs(input_dir)

    infercnv_defaults = config.get("infercnv", {})
    default_cutoff = float(infercnv_defaults.get("cutoff", 0.1))
    default_refs = infercnv_defaults.get("ref_groups", [])
    default_hmm = parse_bool(
        infercnv_defaults.get("HMM", infercnv_defaults.get("hmm", False)),
        "infercnv.HMM",
    )

    sample_cfg = config.get("samples", {})
    runs: list[SampleRun] = []

    for sample, cfg in sample_cfg.items():
        if sample not in discovered:
            raise FileNotFoundError(
                f"Sample '{sample}' is in config but no matching input files were found in {input_dir}"
            )

        sample_files = discovered[sample]
        if "counts_h5" not in sample_files:
            raise FileNotFoundError(f"Missing counts H5 for sample '{sample}'")
        if "annotations_tsv" not in sample_files:
            raise FileNotFoundError(f"Missing annotations TSV for sample '{sample}'")

        ref_groups = cfg.get("ref_groups", default_refs)
        if not isinstance(ref_groups, list) or not all(isinstance(x, str) for x in ref_groups):
            raise ValueError(f"'ref_groups' for sample '{sample}' must be a list of strings")

        cutoff = float(cfg.get("cutoff", default_cutoff))
        hmm = parse_bool(cfg.get("HMM", cfg.get("hmm", default_hmm)), f"samples.{sample}.HMM")

        runs.append(
            SampleRun(
                sample=sample,
                counts_h5=sample_files["counts_h5"],
                annotations_tsv=sample_files["annotations_tsv"],
                out_dir=(output_root / sample),
                ref_groups=ref_groups,
                cutoff=cutoff,
                hmm=hmm,
            )
        )

    return runs


def write_manifest(runs: list[SampleRun], manifest_path: Path) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    with open(manifest_path, "w") as f:
        # header
        f.write("\t".join([
            "sample",
            "counts_h5",
            "annotations_tsv",
            "out_dir",
            "cutoff",
            "hmm",
            "ref_groups",
        ]) + "\n")

        for run in runs:
            f.write("\t".join([
                run.sample,
                str(run.counts_h5),
                str(run.annotations_tsv),
                str(run.out_dir),
                str(run.cutoff),
                str(run.hmm).upper(),
                ",".join(g.replace(" ", "__SPACE__") for g in run.ref_groups) if run.ref_groups else "",
            ]) + "\n")


def write_slurm_array_script(
    config: dict[str, Any],
    runs: list[SampleRun],
    manifest_path: Path,
    r_script_path: Path,
    slurm_script_path: Path,
) -> None:
    slurm = config.get("slurm", {})
    infercnv_cfg = config.get("infercnv", {})

    account = slurm["account"]
    partition = slurm.get("partition", "small")
    time = slurm.get("time", "08:00:00")
    cpus = int(slurm.get("cpus_per_task", 8))
    mem = slurm.get("mem", "64G")
    r_module = slurm.get("r_module", "r-env")
    job_name = slurm.get("job_name", "infercnv")
    array_max = slurm.get("array_max_concurrent", 8)
    project_r_lib = slurm.get("project_r_lib", "")

    gene_order_file = Path(config["gene_order_file"]).resolve()
    cluster_by_groups = str(
        parse_bool(infercnv_cfg.get("cluster_by_groups", True), "infercnv.cluster_by_groups")
    ).upper()
    denoise = str(
        parse_bool(infercnv_cfg.get("denoise", True), "infercnv.denoise")
    ).upper()
    num_threads = int(infercnv_cfg.get("num_threads", cpus))
    num_jobs = len(runs)

    script = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --account={account}
#SBATCH --partition={partition}
#SBATCH --time={time}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --array=1-{num_jobs}%{array_max}
#SBATCH --output=cnv_inference/logs/%x_%A_%a.out
#SBATCH --error=cnv_inference/logs/%x_%A_%a.err

set -euo pipefail

module load {r_module}

MANIFEST="{manifest_path}"
R_SCRIPT="{r_script_path}"
GENE_ORDER="{gene_order_file}"

mkdir -p logs
mkdir -p "$(dirname "$MANIFEST")"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$MANIFEST")

IFS=$'\\t' read -r SAMPLE COUNTS_H5 ANNOTATIONS_TSV OUT_DIR CUTOFF HMM REF_GROUPS <<< "$LINE"

mkdir -p "$OUT_DIR"
"""

    if project_r_lib:
        script += f"""
export R_LIBS_USER="{project_r_lib}:$R_LIBS_USER"
"""

    script += f"""
echo "Running sample: $SAMPLE"
echo "Counts: $COUNTS_H5"
echo "Annotations: $ANNOTATIONS_TSV"
echo "Out dir: $OUT_DIR"
echo "Reference groups: $REF_GROUPS"
echo "HMM: $HMM"

cmd=(
  Rscript "$R_SCRIPT"
  --counts_h5 "$COUNTS_H5"
  --annotations_tsv "$ANNOTATIONS_TSV"
  --gene_order_file "$GENE_ORDER"
  --out_dir "$OUT_DIR"
  --cutoff "$CUTOFF"
  --cluster_by_groups TRUE
  --denoise TRUE
  --num_threads 8
)

# ---------- optional ref_groups ----------
if [[ -n "${{REF_GROUPS:-}}" ]]; then
  cmd+=("--ref_groups=$REF_GROUPS")
fi

# ---------- optional HMM ----------
if [[ -n "${{HMM:-}}" ]]; then
  case "$HMM" in
    TRUE|FALSE)
      cmd+=(--HMM "$HMM")
      ;;
    *)
      echo "Error: HMM must be TRUE, FALSE, or empty in manifest, got: '$HMM'" >&2
      exit 1
      ;;
  esac
fi

printf 'Command arg: <%s>\\n' "${{cmd[@]}}"

"${{cmd[@]}}"
"""

    slurm_script_path.parent.mkdir(parents=True, exist_ok=True)
    slurm_script_path.write_text(script)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to JSON config")
    parser.add_argument(
        "--r-script",
        default="scripts/run_infercnv_sample.R",
        help="Path to the R runner script",
    )
    parser.add_argument(
        "--manifest",
        default=None,
        help="Optional output path for manifest.tsv",
    )
    parser.add_argument(
        "--slurm-script",
        default=None,
        help="Optional output path for submit_infercnv_array.sh",
    )
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    config = load_config(config_path)
    runs = build_runs(config)

    output_root = Path(config["output_root"]).resolve()
    manifest_path = Path(args.manifest).resolve() if args.manifest else output_root / "manifest.tsv"
    slurm_script_path = (
        Path(args.slurm_script).resolve()
        if args.slurm_script
        else output_root / "submit_infercnv_array.sh"
    )
    r_script_path = Path(args.r_script).resolve()

    output_root.mkdir(parents=True, exist_ok=True)
    for run in runs:
        run.out_dir.mkdir(parents=True, exist_ok=True)

    write_manifest(runs, manifest_path)
    write_slurm_array_script(config, runs, manifest_path, r_script_path, slurm_script_path)

    print(f"Wrote manifest: {manifest_path}")
    print(f"Wrote Slurm script: {slurm_script_path}")
    print("Samples:")
    for run in runs:
        print(f"  - {run.sample}: refs={run.ref_groups}, cutoff={run.cutoff}, HMM={run.hmm}")


if __name__ == "__main__":
    main()


    