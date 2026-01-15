#!/usr/bin/env python3
"""
generate_step1_sort_scripts.py

Purpose
-------
Generate (and optionally submit) SLURM job scripts to run:
    funannotate sort --minlen 1000
for each genome FASTA in your input set.

This is Step 1 of the bulk genome annotation pipeline:
    Input:  needs_annotation/fna_input/<OME>.fna
    Output: needs_annotation/working_files/<OME>/prep/<OME>.sort.fasta

What this script does
---------------------
1) Determines which genomes to process by either:
   - Scanning an input directory for *.fna files (default), OR
   - Reading a user-provided list of OME codes (one per line), OR
   - Reading a user-provided list of genome file paths (one per line)

2) Skips genomes already marked complete for Step 1 in:
      progress/annotation_master_progress.tsv
   unless you use --include_done.

3) Writes one SLURM script per genome into:
      scripts/step1_sort/

4) Optionally submits the scripts with sbatch (default: submit).
   You can disable submission with --no_submit.

Typical usage
-------------
Default behavior (scan default input dir and submit jobs):
    python generate_step1_sort_scripts.py

Scan a custom input directory:
    python generate_step1_sort_scripts.py --fna_input_dir /path/to/fna_input

Provide an explicit list of OME codes (one per line):
    python generate_step1_sort_scripts.py --ome_list omes_to_add.txt

Provide an explicit list of genome file paths (one per line):
    python generate_step1_sort_scripts.py --genome_paths genome_paths.txt

Generate scripts but DO NOT submit:
    python generate_step1_sort_scripts.py --no_submit

Generate scripts and submit, but do not filter out already-done genomes:
    python generate_step1_sort_scripts.py --include_done

Notes / assumptions
-------------------
- By default, genomes are expected to be named <OME>.fna in the fna_input_dir.
- This script updates/creates a progress TSV:
      progress/annotation_master_progress.tsv
- The "genus" column is preserved but not filled by this step.
- You can tune SLURM resources via CLI flags (time/cpus/mem/account).

"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
        description="Generate (and optionally submit) SLURM scripts for funannotate sort (Step 1)."
    )

    # Core project structure
    p.add_argument(
        "--project_dir",
        default="/project/arsef/projects/bulk_genome_annotation",
        help="Base project directory (default: /project/arsef/projects/bulk_genome_annotation).",
    )
    p.add_argument(
        "--fna_input_dir",
        default=None,
        help="Directory containing input .fna files. Defaults to <project_dir>/needs_annotation/fna_input",
    )

    # Input selection options
    group = p.add_mutually_exclusive_group()
    group.add_argument(
        "--ome_list",
        default=None,
        help="Path to a text file with one OME code per line (e.g., OME_00001).",
    )
    group.add_argument(
        "--genome_paths",
        default=None,
        help="Path to a text file with one genome file path per line (must exist).",
    )

    # Behavior toggles
    p.add_argument(
        "--include_done",
        action="store_true",
        help="Include genomes even if step1_done is already recorded in the progress file.",
    )
    p.add_argument(
        "--no_submit",
        action="store_true",
        help="Generate scripts but do NOT submit them via sbatch.",
    )

    # SLURM/resource options
    p.add_argument("--time", default="01:00:00", help="SLURM time (default: 01:00:00).")
    p.add_argument("--cpus", type=int, default=2, help="CPUs per task (default: 2).")
    p.add_argument("--mem", default="4000MB", help="Memory per job (default: 4000MB).")
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
        "--minlen",
        type=int,
        default=1000,
        help="Minimum contig length for funannotate sort (default: 1000).",
    )

    return p.parse_args()


def load_or_init_progress(progress_file: Path) -> pd.DataFrame:
    if progress_file.exists():
        df = pd.read_csv(progress_file, sep="\t", dtype=str).fillna("")
        # Ensure required columns exist (in case the schema evolved)
        for col in DEFAULT_PROGRESS_COLUMNS:
            if col not in df.columns:
                df[col] = ""
        # Keep only known columns first, then any extras
        ordered = DEFAULT_PROGRESS_COLUMNS + [c for c in df.columns if c not in DEFAULT_PROGRESS_COLUMNS]
        return df[ordered]
    else:
        return pd.DataFrame(columns=DEFAULT_PROGRESS_COLUMNS)


def normalize_ome_from_filename(path: Path) -> str:
    """
    Extract OME code from filename (expects <OME>.fna), but will also
    accept other extensions and strip them.
    """
    name = path.name
    # Strip common extensions
    for ext in [".fna", ".fa", ".fasta", ".fas"]:
        if name.endswith(ext):
            return name[: -len(ext)]
    return path.stem


def read_ome_list_file(fp: Path) -> List[str]:
    omes: List[str] = []
    with fp.open("r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            omes.append(line)
    return omes


def read_genome_paths_file(fp: Path) -> List[Path]:
    paths: List[Path] = []
    with fp.open("r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            paths.append(Path(line))
    return paths


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


def get_done_omes(df: pd.DataFrame) -> List[str]:
    done = df.loc[df["step1_done"] != "", "OMEcode"].astype(str).tolist()
    return [x for x in done if x]



def build_slurm_script(
    ome: str,
    project_dir: Path,
    input_fna: Path,
    slurm_out: Path,
    slurm_err: Path,
    time: str,
    cpus: int,
    mem: str,
    account: str,
    module_load: str,
    conda_env: str,
    minlen: int,
) -> str:
    work_prep_dir = project_dir / "needs_annotation" / "working_files" / ome / "prep"
    out_fasta = work_prep_dir / f"{ome}.sort.fasta"

    # Use absolute paths in the SLURM script to reduce ambiguity
    return f"""#!/bin/bash
#SBATCH --time={time}
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --account={account}
#SBATCH --job-name=sort_{ome}
#SBATCH --output={slurm_out}
#SBATCH --error={slurm_err}

set -euo pipefail

{module_load}
source activate {conda_env}

mkdir -p "{work_prep_dir}"

funannotate sort -i "{input_fna}" \\
  -o "{out_fasta}" \\
  --minlen {minlen}
"""


def sbatch_submit(script_path: Path) -> Tuple[bool, str]:
    """
    Returns (success, job_id_or_error).
    """
    try:
        result = subprocess.run(
            ["sbatch", str(script_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        stdout = (result.stdout or "").strip()
        # Typical sbatch output: "Submitted batch job 123456"
        job_id = stdout.split()[-1] if stdout else ""
        if not re.fullmatch(r"\d+", job_id):
            # If parsing fails, still return stdout for debugging
            return True, stdout
        return True, job_id
    except subprocess.CalledProcessError as e:
        err = (e.stderr or "").strip() or (e.stdout or "").strip() or str(e)
        return False, err


def main() -> None:
    args = parse_args()

    project_dir = Path(args.project_dir).resolve()
    fna_input_dir = Path(args.fna_input_dir).resolve() if args.fna_input_dir else (project_dir / "needs_annotation" / "fna_input")

    script_output_dir = project_dir / "scripts" / "step1_sort"
    log_output_dir = project_dir / "logs" / "step1_sort"
    progress_file = project_dir / "progress" / "annotation_master_progress.tsv"

    script_output_dir.mkdir(parents=True, exist_ok=True)
    log_output_dir.mkdir(parents=True, exist_ok=True)
    progress_file.parent.mkdir(parents=True, exist_ok=True)

    progress_df = load_or_init_progress(progress_file)

    # Determine targets
    ome_to_input: Dict[str, Path] = {}

    if args.genome_paths:
        paths_file = Path(args.genome_paths).resolve()
        genome_paths = read_genome_paths_file(paths_file)
        for gp in genome_paths:
            if not gp.exists():
                raise FileNotFoundError(f"Genome path does not exist: {gp}")
            ome = normalize_ome_from_filename(gp)
            ome_to_input[ome] = gp.resolve()

    elif args.ome_list:
        omes_file = Path(args.ome_list).resolve()
        omes = read_ome_list_file(omes_file)
        for ome in omes:
            fp = (fna_input_dir / f"{ome}.fna").resolve()
            if not fp.exists():
                raise FileNotFoundError(
                    f"Expected input genome not found for OME '{ome}': {fp}\n"
                    f"Either place it in the input dir or use --genome_paths."
                )
            ome_to_input[ome] = fp

    else:
        if not fna_input_dir.exists():
            raise FileNotFoundError(f"Input directory does not exist: {fna_input_dir}")
        for fp in sorted(fna_input_dir.glob("*.fna")):
            ome = normalize_ome_from_filename(fp)
            ome_to_input[ome] = fp.resolve()

    ome_list = sorted(ome_to_input.keys())

    # Filter out already-done unless requested
    if not args.include_done:
        done_omes = set(get_done_omes(progress_df))
        ome_list = [ome for ome in ome_list if ome not in done_omes]

    if not ome_list:
        print("No genomes to process (after filtering).")
        print(f"Progress file: {progress_file}")
        return

    # Generate scripts (and optionally submit)
    for ome in ome_list:
        input_fna = ome_to_input[ome]
        slurm_script_path = script_output_dir / f"{ome}_sort.sh"
        slurm_out = log_output_dir / f"{ome}_%j.out"
        slurm_err = log_output_dir / f"{ome}_%j.err"

        script_content = build_slurm_script(
            ome=ome,
            project_dir=project_dir,
            input_fna=input_fna,
            slurm_out=slurm_out,
            slurm_err=slurm_err,
            time=args.time,
            cpus=args.cpus,
            mem=args.mem,
            account=args.account,
            module_load=args.module_load,
            conda_env=args.conda_env,
            minlen=args.minlen,
        )

        slurm_script_path.write_text(script_content)

        if args.no_submit:
            # Record that scripts exist, but no job was submitted
            progress_df = update_progress_row(
                progress_df,
                ome,
                {
                    "step1_job": "",
                    "step1_done": "",
                    "note": "Step1 script generated (not submitted)",
                },
            )
            print(f"Generated script (not submitted): {ome} -> {slurm_script_path}")
            continue

        ok, job_or_err = sbatch_submit(slurm_script_path)
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        if ok:
            # If we got a numeric job id, store it; otherwise store stdout text
            job_id = job_or_err
            progress_df = update_progress_row(
                progress_df,
                ome,
                {
                    "step1_job": job_id,
                    "step1_done": "",
                    "note": "",
                },
            )
            print(f"Submitted {ome}: {job_id}")
        else:
            progress_df = update_progress_row(
                progress_df,
                ome,
                {
                    "step1_job": "FAILED",
                    "step1_done": "FAILED",
                    "note": f"Submission failed: {job_or_err}",
                },
            )
            print(f"Failed to submit {ome}: {job_or_err}")

    progress_df.to_csv(progress_file, sep="\t", index=False)
    print(f"\nProgress tracker updated: {progress_file}")


if __name__ == "__main__":
    main()
