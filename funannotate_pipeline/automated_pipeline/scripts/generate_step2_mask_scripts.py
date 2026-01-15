#!/usr/bin/env python3
"""
generate_step2_mask_scripts.py

Purpose
-------
Generate (and optionally submit) SLURM job scripts for Step 2:
  - Build RepeatModeler database
  - Run RepeatModeler
  - Copy key RepeatModeler outputs back into each genome's prep folder
  - Softmask the sorted assembly using funannotate mask

Step 2 expects Step 1 outputs:
  needs_annotation/working_files/<OME>/prep/<OME>.sort.fasta

Step 2 produces:
  needs_annotation/working_files/<OME>/prep/<OME>.masked.fasta
  needs_annotation/working_files/<OME>/prep/<OME>-families.fa
  needs_annotation/working_files/<OME>/prep/<OME>-families.stk
  needs_annotation/working_files/<OME>/prep/<OME>-rmod.log

RepeatModeler working directory (large outputs) goes to:
  /90daydata/arsef/bulk_genome_annotation/repmod/<OME>/

Typical usage
-------------
Run for all genomes ready for Step 2 (based on genomes marked as step1-complete and step2-empty in progress TSV):
    python generate_step2_mask_scripts.py

Run only for a specific set of OMEs (one per line):
    python generate_step2_mask_scripts.py \
      --ome_list /project/arsef/projects/bulk_genome_annotation/commands/fusarium_input_omes.txt

Generate scripts but do NOT submit:
    python generate_step2_mask_scripts.py --no_submit

Include genomes even if step2_done is already set:
    python generate_step2_mask_scripts.py --include_done

Notes
-----
- This script uses the master progress TSV:
    progress/annotation_master_progress.tsv
- By default, it submits jobs via sbatch.
- It will skip a genome if BOTH the repeat library and masked fasta exist
  (unless you use --force).
"""

from __future__ import annotations

import argparse
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd


DEFAULT_PROGRESS_COLUMNS = [
    "OMEcode", "genus",
    "step1_job", "step1_done",
    "step2_job", "step2_done",
    "step3_job", "step3_done",
    "step4_job", "step4_done",
    "note",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate (and optionally submit) SLURM scripts for RepeatModeler + funannotate mask (Step 2)."
    )

    p.add_argument(
        "--project_dir",
        default="/project/arsef/projects/bulk_genome_annotation",
        help="Base project directory (default: /project/arsef/projects/bulk_genome_annotation).",
    )

    # Optional targeted selection
    p.add_argument(
        "--ome_list",
        default=None,
        help="Path to a text file with one OME code per line. If omitted, uses all genomes ready per progress TSV.",
    )

    # Behavior toggles
    p.add_argument(
        "--include_done",
        action="store_true",
        help="Include genomes even if step2_done is already recorded.",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Do not skip genomes even if masked + repeat library files already exist.",
    )
    p.add_argument(
        "--no_submit",
        action="store_true",
        help="Generate scripts but do NOT submit them via sbatch.",
    )

    # RepeatModeler scratch base
    p.add_argument(
        "--repmod_base",
        default="/90daydata/arsef/bulk_genome_annotation/repmod",
        help="Base directory for RepeatModeler working dirs (default: /90daydata/arsef/bulk_genome_annotation/repmod).",
    )

    # SLURM options (match your existing defaults)
    p.add_argument("--time", default="168:00:00", help="SLURM time (default: 168:00:00).")
    p.add_argument("--nodes", type=int, default=1, help="SLURM nodes (default: 1).")
    p.add_argument("--ntasks_per_node", type=int, default=48, help="Tasks per node (default: 48).")
    p.add_argument("--mem_per_cpu", default="5000MB", help="Memory per CPU (default: 5000MB).")
    p.add_argument("--partition", default="ceres", help="SLURM partition (default: ceres).")
    p.add_argument("--account", default="arsef", help="SLURM account (default: arsef).")

    # Environment / command options
    p.add_argument(
        "--module_load",
        default="module load miniconda/24.7.1-2",
        help='Module load line to include in SLURM script (default: "module load miniconda/24.7.1-2").',
    )
    p.add_argument(
        "--conda_env",
        default="/project/arsef/environments/funannotate",
        help="Conda environment path to activate (default: /project/arsef/environments/funannotate).",
    )
    p.add_argument(
        "--blast_usage_report",
        default="false",
        choices=["true", "false"],
        help="Set BLAST_USAGE_REPORT (default: false).",
    )

    return p.parse_args()


def load_or_init_progress(progress_file: Path) -> pd.DataFrame:
    if progress_file.exists():
        df = pd.read_csv(progress_file, sep="\t", dtype=str).fillna("")
        for col in DEFAULT_PROGRESS_COLUMNS:
            if col not in df.columns:
                df[col] = ""
        ordered = DEFAULT_PROGRESS_COLUMNS + [c for c in df.columns if c not in DEFAULT_PROGRESS_COLUMNS]
        return df[ordered]
    return pd.DataFrame(columns=DEFAULT_PROGRESS_COLUMNS)


def read_ome_list_file(fp: Path) -> List[str]:
    omes: List[str] = []
    with fp.open("r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            omes.append(line)
    return omes


def update_progress_row(df: pd.DataFrame, ome: str, col_val_dict: Dict[str, str]) -> pd.DataFrame:
    if (df["OMEcode"] == ome).any():
        for col, val in col_val_dict.items():
            df.loc[df["OMEcode"] == ome, col] = val
    else:
        new_row = {c: "" for c in df.columns}
        new_row["OMEcode"] = ome
        for col, val in col_val_dict.items():
            new_row[col] = val
        df.loc[len(df)] = new_row
    return df


def sbatch_submit(script_path: Path) -> Tuple[bool, str]:
    try:
        result = subprocess.run(
            ["sbatch", str(script_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        stdout = (result.stdout or "").strip()
        job_id = stdout.split()[-1] if stdout else ""
        if not re.fullmatch(r"\d+", job_id):
            return True, stdout
        return True, job_id
    except subprocess.CalledProcessError as e:
        err = (e.stderr or "").strip() or (e.stdout or "").strip() or str(e)
        return False, err


def is_ready_for_step2(row: pd.Series) -> bool:
    # ready means step1_done is set (non-empty)
    return str(row.get("step1_done", "")).strip() != ""


def step2_already_done(row: pd.Series) -> bool:
    return str(row.get("step2_done", "")).strip() != ""


def build_slurm_script(
    ome: str,
    project_dir: Path,
    repmod_dir: Path,
    prep_dir: Path,
    sort_input: Path,
    slurm_out: Path,
    slurm_err: Path,
    args: argparse.Namespace,
) -> str:
    repeat_fa = prep_dir / f"{ome}-families.fa"

    return f"""#!/bin/bash
#SBATCH --time={args.time}
#SBATCH --nodes={args.nodes}
#SBATCH --ntasks-per-node={args.ntasks_per_node}
#SBATCH --mem-per-cpu={args.mem_per_cpu}
#SBATCH --partition={args.partition}
#SBATCH --account={args.account}
#SBATCH --job-name=mask_{ome}
#SBATCH --output={slurm_out}
#SBATCH --error={slurm_err}

set -euo pipefail

{args.module_load}
source activate {args.conda_env}

# Disable/enable NCBI BLAST usage reporting to avoid delays
export BLAST_USAGE_REPORT={args.blast_usage_report}

mkdir -p "{repmod_dir}"
cd "{repmod_dir}"

# Build RepeatModeler database
BuildDatabase -name "{ome}" "{sort_input}"

# Run RepeatModeler
RepeatModeler -threads {args.ntasks_per_node} -database "{ome}" -engine ncbi

# Copy key outputs back to prep dir
cp "{repmod_dir}/{ome}-families.fa" "{repmod_dir}/{ome}-families.stk" "{repmod_dir}/{ome}-rmod.log" "{prep_dir}/"

# Softmask with funannotate (uses the RepeatModeler library)
cd "{prep_dir}"
funannotate mask -i "{sort_input}" -m repeatmodeler -l "{repeat_fa}" -o "{ome}.masked.fasta" --cpus {args.ntasks_per_node}
"""


def main() -> None:
    args = parse_args()

    project_dir = Path(args.project_dir).resolve()
    working_base = project_dir / "needs_annotation" / "working_files"
    repmod_base = Path(args.repmod_base).resolve()

    script_output_dir = project_dir / "scripts" / "step2_mask"
    log_output_dir = project_dir / "logs" / "step2_mask"
    progress_file = project_dir / "progress" / "annotation_master_progress.tsv"

    script_output_dir.mkdir(parents=True, exist_ok=True)
    log_output_dir.mkdir(parents=True, exist_ok=True)
    progress_file.parent.mkdir(parents=True, exist_ok=True)

    progress_df = load_or_init_progress(progress_file)
    progress_df.fillna("", inplace=True)

    # Determine target OMEs
    if args.ome_list:
        ome_list = read_ome_list_file(Path(args.ome_list).resolve())
        # Ensure those OMEs exist in progress, but don’t require they already be present
        for ome in ome_list:
            if not (progress_df["OMEcode"] == ome).any():
                progress_df = update_progress_row(progress_df, ome, {"note": "Added for step2 selection"})
    else:
        # Use all ready-for-step2 based on progress TSV
        ome_list = progress_df.loc[progress_df.apply(is_ready_for_step2, axis=1), "OMEcode"].astype(str).tolist()

    # Filter out already-done step2 unless requested
    if not args.include_done:
        ome_list = [
            ome for ome in ome_list
            if not step2_already_done(progress_df.loc[progress_df["OMEcode"] == ome].iloc[0])
        ]

    if not ome_list:
        print("No genomes to process for Step 2 (after filtering).")
        print(f"Progress file: {progress_file}")
        return

    for ome in ome_list:
        row = progress_df.loc[progress_df["OMEcode"] == ome].iloc[0]

        if not is_ready_for_step2(row):
            print(f"Skipping {ome}: not ready for Step 2 (step1_done is empty).")
            continue

        genome_dir = working_base / ome
        prep_dir = genome_dir / "prep"
        repmod_dir = repmod_base / ome

        sort_input = prep_dir / f"{ome}.sort.fasta"
        masked_output = prep_dir / f"{ome}.masked.fasta"
        repeat_fa = prep_dir / f"{ome}-families.fa"

        if not sort_input.exists():
            print(f"Skipping {ome}: missing Step 1 output: {sort_input}")
            progress_df = update_progress_row(progress_df, ome, {"note": "Missing Step 1 sorted fasta for Step 2"})
            continue

        if (not args.force) and repeat_fa.exists() and masked_output.exists():
            print(f"[✓] Skipping {ome} — repeat + mask already exist.")
            continue

        slurm_script_path = script_output_dir / f"{ome}_step2_mask.sh"
        slurm_out = log_output_dir / f"{ome}_%j.out"
        slurm_err = log_output_dir / f"{ome}_%j.err"

        script_text = build_slurm_script(
            ome=ome,
            project_dir=project_dir,
            repmod_dir=repmod_dir,
            prep_dir=prep_dir,
            sort_input=sort_input,
            slurm_out=slurm_out,
            slurm_err=slurm_err,
            args=args,
        )

        slurm_script_path.write_text(script_text)

        if args.no_submit:
            progress_df = update_progress_row(
                progress_df,
                ome,
                {"step2_job": "", "step2_done": "", "note": "Step2 script generated (not submitted)"},
            )
            print(f"Generated script (not submitted): {ome} -> {slurm_script_path}")
            continue

        ok, job_or_err = sbatch_submit(slurm_script_path)
        if ok:
            progress_df = update_progress_row(progress_df, ome, {"step2_job": job_or_err, "note": ""})
            print(f"[✓] Submitted step2 for {ome}: Job {job_or_err}")
        else:
            progress_df = update_progress_row(
                progress_df,
                ome,
                {"step2_job": "FAILED", "step2_done": "FAILED", "note": f"Submission failed: {job_or_err}"},
            )
            print(f"[✗] Failed to submit step2 for {ome}: {job_or_err}")

    progress_df.to_csv(progress_file, sep="\t", index=False)
    print(f"\nUpdated: {progress_file}")


if __name__ == "__main__":
    main()
