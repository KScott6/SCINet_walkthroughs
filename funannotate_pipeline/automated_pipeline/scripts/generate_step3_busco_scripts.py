#!/usr/bin/env python3
"""
generate_step3_busco_scripts.py

Purpose
-------
Generate (and optionally submit) SLURM job scripts for Step 3:
  Run BUSCO in genome mode on the Step 1 sorted assembly, and (optionally)
  copy the trained AUGUSTUS species model into the shared AUGUSTUS config
  inside your funannotate environment.

Step 3 expects Step 1 outputs:
  needs_annotation/working_files/<OME>/prep/<OME>.sort.fasta

Step 3 produces:
  needs_annotation/working_files/<OME>/busco/<OME>_prelim/...

Typical usage
-------------
Run for all genomes ready for Step 3 (based on progress TSV):
    python generate_step3_busco_scripts.py

Run only for a specific set of OMEs (one per line):
    python generate_step3_busco_scripts.py \
      --ome_list /project/arsef/projects/bulk_genome_annotation/commands/fusarium_input_omes.txt

Generate scripts but do NOT submit:
    python generate_step3_busco_scripts.py --no_submit

Include genomes even if step3_done is already set:
    python generate_step3_busco_scripts.py --include_done

Notes / cautions
----------------
- Step 3 is intentionally NOT dependent on Step 2 (repeat masking).
- The default behavior copies trained AUGUSTUS species models into the
  shared AUGUSTUS config inside the conda environment. This can be risky
  if many BUSCO jobs run simultaneously. If you want BUSCO QC only (no
  shared env modification), use:
      --no_augustus_copy
"""

from __future__ import annotations

import argparse
import re
import subprocess
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
        description="Generate (and optionally submit) SLURM scripts for BUSCO + AUGUSTUS training (Step 3)."
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

    # BUSCO settings
    p.add_argument(
        "--busco_lineage",
        default="hypocreales_odb10",
        help="BUSCO lineage name (default: hypocreales_odb10).",
    )
    p.add_argument(
        "--busco_db_base",
        default="/project/arsef/databases/funannotate_databases",
        help="Base folder containing BUSCO lineage DBs (default: /project/arsef/databases/funannotate_databases).",
    )

    # Behavior toggles
    p.add_argument(
        "--include_done",
        action="store_true",
        help="Include genomes even if step3_done is already recorded (useful for re-runs).",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Do not skip genomes even if BUSCO output directory already exists.",
    )
    p.add_argument(
        "--no_submit",
        action="store_true",
        help="Generate scripts but do NOT submit them via sbatch.",
    )
    p.add_argument(
        "--no_augustus_copy",
        action="store_true",
        help="Do not copy the trained AUGUSTUS species model into the shared AUGUSTUS config.",
    )

    # SLURM options (match your original defaults)
    p.add_argument("--time", default="48:00:00", help="SLURM time (default: 48:00:00).")
    p.add_argument("--nodes", type=int, default=1, help="SLURM nodes (default: 1).")
    p.add_argument("--ntasks_per_node", type=int, default=16, help="Tasks per node (default: 16).")
    p.add_argument("--mem_per_cpu", default="3000MB", help="Memory per CPU (default: 3000MB).")
    p.add_argument("--partition", default="ceres", help="SLURM partition (default: ceres).")
    p.add_argument("--account", default="arsef", help="SLURM account (default: arsef).")

    # Environment
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


def is_ready_for_step3(row: pd.Series) -> bool:
    # Step 3 depends only on Step 1 completion
    return str(row.get("step1_done", "")).strip() != ""


def step3_already_done(row: pd.Series) -> bool:
    val = str(row.get("step3_done", "")).strip()
    return val != "" and val != "FAILED"


def build_slurm_script(
    ome: str,
    sort_input: Path,
    busco_output_dir: Path,
    busco_lineage: str,
    busco_db_path: Path,
    slurm_out: Path,
    slurm_err: Path,
    args: argparse.Namespace,
) -> str:
    # BUSCO output name and retraining path
    out_name = f"{ome}_prelim"
    busco_output_path = busco_output_dir / out_name

    # This matches your previous layout
    retrain_path = busco_output_path / f"run_{busco_lineage}" / "augustus_output" / "retraining_parameters" / f"BUSCO_{ome}_prelim"

    # Optional AUGUSTUS copy block
    if args.no_augustus_copy:
        augustus_block = f"""
echo "[i] --no_augustus_copy set: skipping AUGUSTUS species copy for {ome}"
"""
    else:
        augustus_block = f"""
AUGUSTUS_CONFIG_PATH=$(dirname $(dirname $(which augustus)))/config/

# Copy trained species model into shared AUGUSTUS config
if [ -d "{retrain_path}" ]; then
  cp -R "{retrain_path}" "${{AUGUSTUS_CONFIG_PATH}}/species/"
  funannotate setup -u -w
  chmod -R 770 "${{AUGUSTUS_CONFIG_PATH}}/species/BUSCO_{ome}_prelim"
else
  echo "[✗] BUSCO training output not found for {ome}" >&2
  exit 1
fi
"""

    return f"""#!/bin/bash
#SBATCH --time={args.time}
#SBATCH --nodes={args.nodes}
#SBATCH --ntasks-per-node={args.ntasks_per_node}
#SBATCH --mem-per-cpu={args.mem_per_cpu}
#SBATCH --partition={args.partition}
#SBATCH --account={args.account}
#SBATCH --job-name=busco_{ome}
#SBATCH --output={slurm_out}
#SBATCH --error={slurm_err}

set -euo pipefail

{args.module_load}
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate {args.conda_env}

mkdir -p "{busco_output_dir}"
cd "{busco_output_dir}"

busco --long --augustus -i "{sort_input}" \\
  -o "{out_name}" -l "{busco_db_path}" -m genome -c {args.ntasks_per_node}

{augustus_block}
"""


def main() -> None:
    args = parse_args()

    project_dir = Path(args.project_dir).resolve()
    working_base = project_dir / "needs_annotation" / "working_files"
    script_output_dir = project_dir / "scripts" / "step3_busco"
    log_output_dir = project_dir / "logs" / "step3_busco"
    progress_file = project_dir / "progress" / "annotation_master_progress.tsv"

    busco_db_path = Path(args.busco_db_base).resolve() / args.busco_lineage

    script_output_dir.mkdir(parents=True, exist_ok=True)
    log_output_dir.mkdir(parents=True, exist_ok=True)
    progress_file.parent.mkdir(parents=True, exist_ok=True)

    df = load_or_init_progress(progress_file)

    # Determine target OMEs
    if args.ome_list:
        ome_list = read_ome_list_file(Path(args.ome_list).resolve())
        for ome in ome_list:
            if not (df["OMEcode"] == ome).any():
                df = update_progress_row(df, ome, {"note": "Added for step3 selection"})
    else:
        ome_list = df.loc[df.apply(is_ready_for_step3, axis=1), "OMEcode"].astype(str).tolist()

    # Filter out already-done unless requested
    if not args.include_done:
        filtered: List[str] = []
        for ome in ome_list:
            row = df.loc[df["OMEcode"] == ome].iloc[0]
            if not step3_already_done(row):
                filtered.append(ome)
        ome_list = filtered

    if not ome_list:
        print("No genomes to process for Step 3 (after filtering).")
        print(f"Progress file: {progress_file}")
        return

    # Validate BUSCO DB exists (fail early)
    if not busco_db_path.exists():
        raise FileNotFoundError(f"BUSCO DB path does not exist: {busco_db_path}")

    for ome in ome_list:
        row = df.loc[df["OMEcode"] == ome].iloc[0]

        if not is_ready_for_step3(row):
            print(f"Skipping {ome}: not ready for Step 3 (step1_done is empty).")
            continue

        genome_dir = working_base / ome
        prep_dir = genome_dir / "prep"
        sort_input = prep_dir / f"{ome}.sort.fasta"
        busco_output_dir = genome_dir / "busco"

        if not sort_input.exists():
            print(f"[✗] Missing sorted genome for {ome}, skipping: {sort_input}")
            df = update_progress_row(df, ome, {"note": "Missing Step 1 sorted fasta for Step 3"})
            continue

        # Optional skip if output already exists
        expected_dir = busco_output_dir / f"{ome}_prelim"
        if (not args.force) and expected_dir.exists():
            print(f"[✓] Skipping {ome} — BUSCO output already exists: {expected_dir}")
            continue

        slurm_script_path = script_output_dir / f"{ome}_step3_busco.sh"
        slurm_out = log_output_dir / f"{ome}_%j.out"
        slurm_err = log_output_dir / f"{ome}_%j.err"

        script_text = build_slurm_script(
            ome=ome,
            sort_input=sort_input,
            busco_output_dir=busco_output_dir,
            busco_lineage=args.busco_lineage,
            busco_db_path=busco_db_path,
            slurm_out=slurm_out,
            slurm_err=slurm_err,
            args=args,
        )
        slurm_script_path.write_text(script_text)

        if args.no_submit:
            df = update_progress_row(df, ome, {"step3_job": "", "step3_done": "", "note": "Step3 script generated (not submitted)"})
            print(f"Generated script (not submitted): {ome} -> {slurm_script_path}")
            continue

        ok, job_or_err = sbatch_submit(slurm_script_path)
        if ok:
            df = update_progress_row(df, ome, {"step3_job": job_or_err, "step3_done": "", "note": ""})
            print(f"[✓] Submitted BUSCO for {ome}: Job {job_or_err}")
        else:
            df = update_progress_row(df, ome, {"step3_job": "FAILED", "step3_done": "FAILED", "note": f"step3 submission failed: {job_or_err}"})
            print(f"[✗] Failed to submit BUSCO for {ome}: {job_or_err}")

    df.to_csv(progress_file, sep="\t", index=False)
    print(f"\nProgress tracker updated: {progress_file}")


if __name__ == "__main__":
    main()
