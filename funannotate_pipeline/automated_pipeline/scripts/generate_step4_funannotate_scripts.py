#!/usr/bin/env python3
"""
generate_step4_funannotate_scripts.py

Purpose
-------
Generate (and optionally submit) SLURM job scripts for Step 4 (Funannotate predict).

Step 4 expects Step 2 outputs:
  needs_annotation/working_files/<OME>/prep/<OME>.masked.fasta

Step 4 produces (in an intermediate output dir on 90daydata):
  <out_base>/<OME>_output/   (or <out_base>/<OME>/ if you choose)

and then copies key results to final locations:
  final_funannotate_results/by_ome/<OME>/
  final_funannotate_results/busco_short_summaries/
  final_funannotate_results/gff3/

Selection logic
---------------
Default (no --ome_list):
  - Uses the progress TSV to find genomes that are ready for Step 4:
      step2_done is non-empty  AND  step4_done is empty
  - Also requires that the masked fasta exists.

With --ome_list:
  - Restricts to only OMEs in that file, then applies the same readiness checks.

Submission behavior
-------------------
By default this script only GENERATES job scripts (does NOT submit).
To submit jobs automatically, add:
  --submit

Typical usage
-------------
Generate scripts for all ready genomes (do not submit):
  python generate_step4_funannotate_scripts.py

Generate scripts for a specific batch list (do not submit):
  python generate_step4_funannotate_scripts.py \
    --ome_list /project/arsef/projects/bulk_genome_annotation/commands/fusarium_input_omes.txt

Generate scripts and submit:
  python generate_step4_funannotate_scripts.py --submit

Generate scripts + submit only for fusarium batch:
  python generate_step4_funannotate_scripts.py \
    --ome_list /project/arsef/projects/bulk_genome_annotation/commands/fusarium_input_omes.txt \
    --submit

Notes
-----
- UniProt evidence is always included as protein evidence.
- Genus-specific evidence is pulled from:
    <evidence_base>/transcript/<genus>/
    <evidence_base>/protein/<genus>/
- This version creates a unique TMPDIR per job under --tmpdir_base.
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
        description="Generate (and optionally submit) SLURM scripts for Funannotate predict (Step 4)."
    )

    p.add_argument(
        "--project_dir",
        default="/project/arsef/projects/bulk_genome_annotation",
        help="Base project directory.",
    )
    p.add_argument(
        "--progress_file",
        default=None,
        help="Path to annotation_master_progress.tsv (default: <project_dir>/progress/annotation_master_progress.tsv).",
    )

    # Optional selection restriction
    p.add_argument(
        "--ome_list",
        default=None,
        help="Path to a text file with one OME code per line. If omitted, selects all ready genomes from progress TSV.",
    )

    # Output and evidence
    p.add_argument(
        "--evidence_base",
        default=None,
        help="Evidence base directory (default: <project_dir>/evidence).",
    )
    p.add_argument(
        "--uniprot_path",
        default="/project/arsef/databases/funannotate_databases/uniprot_sprot.fasta",
        help="UniProt Swiss-Prot fasta (always included as protein evidence).",
    )

    p.add_argument(
        "--working_base",
        default=None,
        help="Working files base (default: <project_dir>/needs_annotation/working_files).",
    )

    p.add_argument(
        "--out_base",
        default="/90daydata/arsef/output_funannotate",
        help="Base directory for funannotate output (default: /90daydata/arsef/output_funannotate).",
    )
    p.add_argument(
        "--tmpdir_base",
        default="/90daydata/arsef/tmp_funannotate",
        help="Base directory under which per-job TMPDIRs will be created (default: /90daydata/arsef/tmp_funannotate).",
    )

    # Final destinations
    p.add_argument(
        "--final_by_ome_dir",
        default=None,
        help="Final results base folder (default: <project_dir>/needs_annotation/final_funannotate_results/by_ome).",
    )
    p.add_argument(
        "--short_summary_dest",
        default=None,
        help="Destination folder for BUSCO short_summary files (default: <project_dir>/needs_annotation/final_funannotate_results/busco_short_summaries).",
    )
    p.add_argument(
        "--gff_dest",
        default=None,
        help="Destination folder for GFF3 files (default: <project_dir>/needs_annotation/final_funannotate_results/gff3).",
    )

    # Funannotate settings
    p.add_argument(
        "--busco_db",
        default="hypocreales_odb10",
        help="BUSCO lineage name for funannotate predict (default: hypocreales_odb10).",
    )
    p.add_argument(
        "--run_number",
        type=int,
        default=1,
        help="Run number used in sample name (default: 1).",
    )
    p.add_argument(
        "--cpus",
        type=int,
        default=16,
        help="CPUs for funannotate predict (default: 16).",
    )

    # SLURM settings
    p.add_argument("--time", default="72:00:00", help="SLURM time (default: 72:00:00).")
    p.add_argument("--nodes", type=int, default=1, help="SLURM nodes (default: 1).")
    p.add_argument("--ntasks_per_node", type=int, default=16, help="Tasks per node (default: 16).")
    p.add_argument("--mem_per_cpu", default="3000MB", help="Memory per CPU (default: 3000MB).")
    p.add_argument("--partition", default="ceres", help="SLURM partition (default: ceres).")
    p.add_argument("--account", default="arsef", help="SLURM account (default: arsef).")

    # Env / modules (keeps your current conda-based predict)
    p.add_argument(
        "--module_load",
        default="module load miniconda",
        help='Module load line for the SLURM script (default: "module load miniconda").',
    )
    p.add_argument(
        "--conda_env",
        default="/project/arsef/environments/funannotate_working/",
        help="Conda env to activate for funannotate predict (default: /project/arsef/environments/funannotate_working/).",
    )
    p.add_argument(
        "--module_unload",
        default="module unload miniconda",
        help='Optional module unload line (default: "module unload miniconda").',
    )

    # Funannotate env vars (as in your script)
    p.add_argument(
        "--genemark_path",
        default="/project/arsef/environments/__external_software_funannotate/gmes_linux_64_4",
        help="GENEMARK_PATH export value.",
    )
    p.add_argument(
        "--funannotate_db",
        default="/project/arsef/databases/funannotate_databases",
        help="FUNANNOTATE_DB export value.",
    )

    # Behavior toggles
    p.add_argument(
        "--include_done",
        action="store_true",
        help="Include genomes even if step4_done is already set in progress TSV.",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Pass --force to funannotate predict (and allow overwriting existing output dir).",
    )
    p.add_argument(
        "--submit",
        action="store_true",
        help="Submit generated scripts via sbatch. (Default: generate only).",
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


def is_ready_for_step4(row: pd.Series) -> bool:
    # Step 4 depends on Step 2 having completed (masked genome exists)
    return str(row.get("step2_done", "")).strip() != ""


def step4_already_done(row: pd.Series) -> bool:
    val = str(row.get("step4_done", "")).strip()
    return val != "" and val != "FAILED"


def gather_evidence_files(evidence_base: Path, genus: str, uniprot_path: Path) -> Tuple[List[Path], List[Path]]:
    """
    Returns (transcripts, proteins) as Path lists.
    UniProt is NOT included here; caller appends it to proteins.
    """
    transcripts: List[Path] = []
    proteins: List[Path] = []

    if genus:
        transcript_dir = evidence_base / "transcript" / genus
        protein_dir = evidence_base / "protein" / genus

        if transcript_dir.is_dir():
            exts = ["*.fasta", "*.fa", "*.fna", "*.nt.fasta"]
            for ext in exts:
                transcripts.extend([Path(p) for p in transcript_dir.glob(ext)])
        if protein_dir.is_dir():
            exts = ["*.fasta", "*.fa", "*.aa.fasta", "*.faa"]
            for ext in exts:
                proteins.extend([Path(p) for p in protein_dir.glob(ext)])

    # De-duplicate + sort
    transcripts = sorted({p.resolve() for p in transcripts})
    proteins = sorted({p.resolve() for p in proteins})

    return transcripts, proteins


def build_slurm_script(
    ome: str,
    masked_fasta: Path,
    genus: str,
    transcripts: List[Path],
    proteins: List[Path],
    uniprot_path: Path,
    out_dir: Path,
    final_out_dir: Path,
    short_summary_dest: Path,
    gff_dest: Path,
    slurm_out: Path,
    slurm_err: Path,
    args: argparse.Namespace,
) -> str:
    # Build evidence strings (funannotate accepts multiple paths)
    transcript_block = ""
    if transcripts:
        transcript_block = "--transcript_evidence \\\n  " + " \\\n  ".join(str(p) for p in transcripts) + " \\\n  "

    # Always include UniProt in protein evidence
    all_proteins = proteins + [uniprot_path]
    protein_block = "--protein_evidence \\\n  " + " \\\n  ".join(str(p) for p in all_proteins) + " \\\n  "

    force_flag = "--force" if args.force else ""

    # Sample name
    sample_name = f"{ome}_run_{args.run_number}"

    # BUSCO short summary path (funannotate predict_misc layout varies; keep your original target but guarded)
    # Your original: predict_misc/busco/run_<sample>/short_summary_<sample>.txt
    short_summary_src = out_dir / "predict_misc" / "busco" / f"run_{sample_name}" / f"short_summary_{sample_name}.txt"
    gff_src = out_dir / "predict_results" / f"{sample_name}.gff3"

    # Unique tmpdir created per job under tmpdir_base
    tmpdir_job = f'{args.tmpdir_base}/{ome}_${{SLURM_JOB_ID}}'

    return f"""#!/bin/bash
#SBATCH --time={args.time}
#SBATCH --nodes={args.nodes}
#SBATCH --ntasks-per-node={args.ntasks_per_node}
#SBATCH --partition={args.partition}
#SBATCH --account={args.account}
#SBATCH --mem-per-cpu={args.mem_per_cpu}
#SBATCH --job-name={ome}_fun
#SBATCH --output={slurm_out}
#SBATCH --error={slurm_err}

set -euo pipefail

{args.module_load}
source activate {args.conda_env}
{args.module_unload}

export GENEMARK_PATH="{args.genemark_path}"
export FUNANNOTATE_DB="{args.funannotate_db}"

OMEcode="{ome}"
GENUS="{genus}"
BUSCO_LINEAGE="{args.busco_db}"
RUN_NUMBER="{args.run_number}"
CPUS="{args.cpus}"

# Create unique TMPDIR for this job
export TMPDIR="{tmpdir_job}"
mkdir -p "$TMPDIR"
echo "Using TMPDIR: $TMPDIR"

# Prepare output directory
rm -rf "{out_dir}"
mkdir -p "{out_dir}"
cd "{out_dir}"

funannotate predict \\
  -i "{masked_fasta}" \\
  -s "{sample_name}" \\
  {transcript_block}{protein_block}--cpus {args.cpus} \\
  --optimize_augustus \\
  --busco_db "{args.busco_db}" \\
  -o "{out_dir}" \\
  --tmpdir "$TMPDIR" \\
  {force_flag}

# Copy results to final locations
mkdir -p "{final_out_dir}"
mkdir -p "{short_summary_dest}"
mkdir -p "{gff_dest}"

# Copy predict_results contents into per-OME final folder
if [ -d "{out_dir}/predict_results" ]; then
  cp -r "{out_dir}/predict_results/"* "{final_out_dir}/" || true
else
  echo "No predict_results directory found for {ome}"
fi

# Copy BUSCO short summary (guarded)
if [ -f "{short_summary_src}" ]; then
  cp "{short_summary_src}" "{short_summary_dest}/" || true
else
  echo "Missing BUSCO short summary for {ome}: {short_summary_src}"
fi

# Copy GFF3 (guarded)
if [ -f "{gff_src}" ]; then
  cp "{gff_src}" "{gff_dest}/" || true
else
  echo "Missing GFF3 for {ome}: {gff_src}"
fi

# Optional: cleanup tmpdir
rm -rf "$TMPDIR" || true
"""


def main() -> None:
    args = parse_args()

    project_dir = Path(args.project_dir).resolve()
    progress_file = Path(args.progress_file).resolve() if args.progress_file else (project_dir / "progress" / "annotation_master_progress.tsv")

    script_output_dir = project_dir / "scripts" / "step4_funannotate"
    log_output_dir = project_dir / "logs" / "step4_funannotate"

    evidence_base = Path(args.evidence_base).resolve() if args.evidence_base else (project_dir / "evidence")
    working_base = Path(args.working_base).resolve() if args.working_base else (project_dir / "needs_annotation" / "working_files")
    out_base = Path(args.out_base).resolve()

    final_by_ome_dir = Path(args.final_by_ome_dir).resolve() if args.final_by_ome_dir else (project_dir / "needs_annotation" / "final_funannotate_results" / "by_ome")
    short_summary_dest = Path(args.short_summary_dest).resolve() if args.short_summary_dest else (project_dir / "needs_annotation" / "final_funannotate_results" / "busco_short_summaries")
    gff_dest = Path(args.gff_dest).resolve() if args.gff_dest else (project_dir / "needs_annotation" / "final_funannotate_results" / "gff3")

    uniprot_path = Path(args.uniprot_path).resolve()

    script_output_dir.mkdir(parents=True, exist_ok=True)
    log_output_dir.mkdir(parents=True, exist_ok=True)
    progress_file.parent.mkdir(parents=True, exist_ok=True)

    df = load_or_init_progress(progress_file)

    # Determine target OMEs
    if args.ome_list:
        ome_list = read_ome_list_file(Path(args.ome_list).resolve())
        for ome in ome_list:
            if not (df["OMEcode"] == ome).any():
                df = update_progress_row(df, ome, {"note": "Added for step4 selection"})
    else:
        # All ready per progress TSV
        ome_list = df.loc[df.apply(is_ready_for_step4, axis=1), "OMEcode"].astype(str).tolist()

    # Filter out already-done unless requested
    if not args.include_done:
        filtered: List[str] = []
        for ome in ome_list:
            row = df.loc[df["OMEcode"] == ome].iloc[0]
            if not step4_already_done(row):
                filtered.append(ome)
        ome_list = filtered

    if not ome_list:
        print("No genomes to process for Step 4 (after filtering).")
        print(f"Progress file: {progress_file}")
        return

    if not uniprot_path.exists():
        raise FileNotFoundError(f"UniProt protein evidence not found: {uniprot_path}")

    for ome in ome_list:
        row = df.loc[df["OMEcode"] == ome].iloc[0]
        genus = str(row.get("genus", "")).strip()

        if not is_ready_for_step4(row):
            print(f"Skipping {ome}: not ready for Step 4 (step2_done is empty).")
            continue

        masked_fasta = working_base / ome / "prep" / f"{ome}.masked.fasta"
        if not masked_fasta.exists():
            print(f"Skipping {ome}: missing masked fasta: {masked_fasta}")
            df = update_progress_row(df, ome, {"note": "Missing masked fasta for Step 4"})
            continue

        transcripts, proteins = gather_evidence_files(evidence_base, genus, uniprot_path)

        if not proteins:
            # You'll still have UniProt added later, but this warns you genus-specific protein evidence is absent
            print(f"[WARNING] No genus-specific protein evidence for genus '{genus}' ({ome}). UniProt will still be used.")

        out_dir = out_base / f"{ome}_output"
        final_out_dir = final_by_ome_dir / ome

        slurm_script_path = script_output_dir / f"{ome}_funannotate.sh"
        slurm_out = log_output_dir / f"{ome}_fun_predict_%j.o"
        slurm_err = log_output_dir / f"{ome}_fun_predict_%j.e"

        script_text = build_slurm_script(
            ome=ome,
            masked_fasta=masked_fasta,
            genus=genus,
            transcripts=transcripts,
            proteins=proteins,
            uniprot_path=uniprot_path,
            out_dir=out_dir,
            final_out_dir=final_out_dir,
            short_summary_dest=short_summary_dest,
            gff_dest=gff_dest,
            slurm_out=slurm_out,
            slurm_err=slurm_err,
            args=args,
        )
        slurm_script_path.write_text(script_text)

        print(f"[DONE] Script written: {slurm_script_path}")

        if args.submit:
            ok, job_or_err = sbatch_submit(slurm_script_path)
            if ok:
                df = update_progress_row(df, ome, {"step4_job": job_or_err, "step4_done": "", "note": ""})
                print(f"[✓] Submitted Step 4 for {ome}: Job {job_or_err}")
            else:
                df = update_progress_row(df, ome, {"step4_job": "FAILED", "step4_done": "FAILED", "note": f"step4 submission failed: {job_or_err}"})
                print(f"[✗] Failed to submit Step 4 for {ome}: {job_or_err}")

    df.to_csv(progress_file, sep="\t", index=False)
    print(f"\nProgress tracker updated: {progress_file}")


if __name__ == "__main__":
    main()
