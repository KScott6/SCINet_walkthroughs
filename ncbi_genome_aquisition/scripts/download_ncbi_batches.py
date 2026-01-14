#!/usr/bin/env python3

import argparse
import os
import re
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def check_exec_exists(exe: str) -> None:
    from shutil import which
    if which(exe) is None:
        die(f"Required executable not found in PATH: '{exe}'. Activate your ncbi_datasets env and try again.")


def run_cmd(cmd: List[str], cwd: Optional[Path] = None) -> Tuple[int, str, str]:
    p = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    return p.returncode, p.stdout, p.stderr


def read_metadata(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in [".tsv", ".tab"]:
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, dtype=str)


def ensure_dirs(base: Path) -> Dict[str, Path]:
    subdirs = {
        "fna": base / "fna",
        "gff": base / "gff",
        "gtf": base / "gtf",
        "faa": base / "faa",
        "cds": base / "cds",
        "rna": base / "rna",
        "logs": base / "logs",
        "dump_site": base / "download_dump_site",
        "chunks": base / "accession_chunks",
    }
    for p in subdirs.values():
        p.mkdir(parents=True, exist_ok=True)
    return subdirs


def chunk_list(items: List[str], chunk_size: int) -> List[List[str]]:
    return [items[i:i + chunk_size] for i in range(0, len(items), chunk_size)]


def write_lines(path: Path, lines: List[str]) -> None:
    path.write_text("\n".join(lines) + ("\n" if lines else ""))


def find_first_matching(directory: Path, patterns: List[str]) -> Optional[Path]:
    """
    Find first file in directory matching any regex pattern in `patterns`.
    Returns Path or None.
    """
    if not directory.exists():
        return None
    for p in directory.iterdir():
        if not p.is_file():
            continue
        name = p.name
        for pat in patterns:
            if re.fullmatch(pat, name):
                return p
    return None


def move_if_exists(src: Path, dest: Path) -> Optional[Path]:
    if src.exists() and src.is_file():
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dest))
        return dest
    return None


def extract_zip(zip_path: Path, dest_dir: Path) -> None:
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(dest_dir)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch-download NCBI genome assemblies from a metadata table (deduped NEW_ONLY.tsv), "
                    "optionally include annotations, and move files into tidy folders (fna/gff/gtf/faa/cds/rna)."
    )
    parser.add_argument("--metadata", required=True,
                        help="Metadata TSV/CSV that includes accessions (default col: assembly_acc). "
                             "Typically the *.NEW_ONLY.tsv output from fetch_ncbi_metadata_and_merge.py")
    parser.add_argument("--accession_col", default="assembly_acc",
                        help="Column containing accessions (default: assembly_acc)")
    parser.add_argument("--outdir", required=True,
                        help="Output directory (will create subfolders + logs + dump_site).")
    parser.add_argument("--chunk_size", type=int, default=10,
                        help="Number of accessions per datasets download batch (default: 10)")
    parser.add_argument("--sleep", type=float, default=5.0,
                        help="Seconds to sleep between downloads (default: 5.0)")
    parser.add_argument("--api_key", default=None,
                        help="NCBI API key for datasets downloads (optional). If not provided, uses env NCBI_API_KEY if set.")
    parser.add_argument("--with_annotation", action="store_true",
                        help="Request annotation-related payloads where available (gff3/gtf/protein/cds/rna).")
    parser.add_argument("--resume", action="store_true",
                        help="If set, skip accessions already listed in processed_accessions.log")
    args = parser.parse_args()

    check_exec_exists("datasets")

    meta_path = Path(args.metadata).resolve()
    if not meta_path.exists():
        die(f"Metadata file not found: {meta_path}")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    dirs = ensure_dirs(outdir)

    # Logs
    manifest_path = outdir / "download_manifest.tsv"
    failed_log = dirs["logs"] / "failed_accessions.log"
    processed_log = dirs["logs"] / "processed_accessions.log"
    stderr_log = dirs["logs"] / "datasets_stderr.log"
    updated_meta_path = outdir / "NEW_ONLY_with_paths.tsv"

    # Initialize logs
    if not manifest_path.exists():
        manifest_path.write_text(
            "accession\tannotation_available\tcds_location\tfna_location\tgff_location\tgtf_location\tfaa_location\trna_location\n"
        )
    if not failed_log.exists():
        failed_log.write_text("")
    if not processed_log.exists():
        processed_log.write_text("")
    if not stderr_log.exists():
        stderr_log.write_text("")

    # Load metadata
    meta = read_metadata(meta_path)
    if args.accession_col not in meta.columns:
        die(f"Accession column '{args.accession_col}' not found. Columns: {list(meta.columns)[:40]} ...")

    accs = (
        meta[args.accession_col]
        .dropna()
        .astype(str)
        .str.strip()
    )
    accs = accs[accs != ""].drop_duplicates().tolist()

    if not accs:
        die("No accessions found after cleaning input column.")

    # Resume logic
    processed_set = set()
    if args.resume and processed_log.exists():
        processed_set = set([l.strip() for l in processed_log.read_text().splitlines() if l.strip()])
        accs = [a for a in accs if a not in processed_set]

    if not accs:
        print("Nothing to do: all accessions already processed (resume mode).")
        return

    # Determine include string
    include_str = "genome"
    if args.with_annotation:
        include_str = "genome,rna,protein,cds,gff3,gtf"

    # API key
    api_key = args.api_key or os.environ.get("NCBI_API_KEY", "")
    api_key_arg = []
    if api_key:
        api_key_arg = ["--api-key", api_key]

    print(f"Accessions to process: {len(accs)}")
    print(f"chunk_size: {args.chunk_size}")
    print(f"with_annotation: {'yes' if args.with_annotation else 'no'}")
    sys.stdout.flush()

    # Chunk accessions
    chunks = chunk_list(accs, args.chunk_size)

    for i, chunk in enumerate(chunks, start=1):
        chunk_file = dirs["chunks"] / f"chunk_{i:04d}.txt"
        write_lines(chunk_file, chunk)

        # datasets download writes ncbi_dataset.zip in current working directory by default
        # We run in dump_site, and then unzip into dump_site.
        zip_path = dirs["dump_site"] / "ncbi_dataset.zip"

        # clean any prior zip
        if zip_path.exists():
            zip_path.unlink()

        cmd = [
            "datasets", "download", "genome", "accession",
            "--inputfile", str(chunk_file),
            "--include", include_str,
            "--filename", str(zip_path),
        ] + api_key_arg

        rc, out, err = run_cmd(cmd, cwd=dirs["dump_site"])

        # Always log stderr for debugging (even if success)
        with stderr_log.open("a") as fh:
            fh.write(f"\n===== CHUNK {i} =====\n")
            fh.write(f"CMD: {' '.join(cmd)}\n")
            fh.write(err.strip() + "\n")

        if rc != 0 or not zip_path.exists():
            # Log all accessions in chunk as failed
            with failed_log.open("a") as fh:
                for a in chunk:
                    fh.write(a + "\n")
            print(f"[WARN] Chunk {i} download failed; logged {len(chunk)} accessions.")
            time.sleep(args.sleep)
            continue

        # Unzip into dump_site
        # Remove any existing ncbi_dataset folder first to avoid confusion
        ncbi_dataset_dir = dirs["dump_site"] / "ncbi_dataset"
        if ncbi_dataset_dir.exists():
            shutil.rmtree(ncbi_dataset_dir)

        extract_zip(zip_path, dirs["dump_site"])
        zip_path.unlink()  # remove zip after extraction

        data_dir = dirs["dump_site"] / "ncbi_dataset" / "data"
        if not data_dir.exists():
            # catastrophic extraction failure; mark chunk failed
            with failed_log.open("a") as fh:
                for a in chunk:
                    fh.write(a + "\n")
            print(f"[WARN] Chunk {i} extracted but data dir missing; logged failures.")
            # cleanup
            if ncbi_dataset_dir.exists():
                shutil.rmtree(ncbi_dataset_dir)
            time.sleep(args.sleep)
            continue

        # Process each accession folder
        for acc in chunk:
            acc_dir = data_dir / acc
            if not acc_dir.exists():
                with failed_log.open("a") as fh:
                    fh.write(acc + "\n")
                continue

            cds_path = "NA"
            fna_path = "NA"
            gff_path = "NA"
            gtf_path = "NA"
            faa_path = "NA"
            rna_path = "NA"

            # cds_from_genomic.fna
            cds_src = acc_dir / "cds_from_genomic.fna"
            cds_dest = dirs["cds"] / f"{acc}_cds.fna"
            moved = move_if_exists(cds_src, cds_dest)
            if moved:
                cds_path = str(moved)

            # genome fasta: there can be multiple .fna; take the first that's not cds_from_genomic.fna
            # Typical name is genomic.fna; sometimes other names exist.
            fna_candidates = []
            for p in acc_dir.glob("*.fna"):
                if p.name == "cds_from_genomic.fna":
                    continue
                fna_candidates.append(p)

            if fna_candidates:
                # Prefer genomic.fna if present
                genomic = acc_dir / "genomic.fna"
                chosen = genomic if genomic in fna_candidates else fna_candidates[0]
                fna_dest = dirs["fna"] / f"{acc}.fna"
                moved = move_if_exists(chosen, fna_dest)
                if moved:
                    fna_path = str(moved)

            # gff / gff3
            gff_src = acc_dir / "genomic.gff"
            if gff_src.exists():
                gff_dest = dirs["gff"] / f"{acc}.gff"
                moved = move_if_exists(gff_src, gff_dest)
                if moved:
                    gff_path = str(moved)

            # gtf
            gtf_src = acc_dir / "genomic.gtf"
            if gtf_src.exists():
                gtf_dest = dirs["gtf"] / f"{acc}.gtf"
                moved = move_if_exists(gtf_src, gtf_dest)
                if moved:
                    gtf_path = str(moved)

            # proteins
            faa_src = acc_dir / "protein.faa"
            if faa_src.exists():
                faa_dest = dirs["faa"] / f"{acc}.faa"
                moved = move_if_exists(faa_src, faa_dest)
                if moved:
                    faa_path = str(moved)

            # rna
            rna_src = acc_dir / "rna.fna"
            if rna_src.exists():
                rna_dest = dirs["rna"] / f"{acc}_rna.fna"
                moved = move_if_exists(rna_src, rna_dest)
                if moved:
                    rna_path = str(moved)

            annotation_available = "yes" if gff_path != "NA" else "no"

            # Write manifest
            with manifest_path.open("a") as fh:
                fh.write(
                    f"{acc}\t{annotation_available}\t{cds_path}\t{fna_path}\t{gff_path}\t{gtf_path}\t{faa_path}\t{rna_path}\n"
                )

            # Mark processed
            with processed_log.open("a") as fh:
                fh.write(acc + "\n")

        # Cleanup extracted dataset folder
        if ncbi_dataset_dir.exists():
            shutil.rmtree(ncbi_dataset_dir)

        time.sleep(args.sleep)

    # Build updated metadata with paths + annotation_available
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    manifest = manifest.drop_duplicates(subset=["accession"], keep="last")

    meta_out = meta.copy()
    # Normalize accession col to string keys
    meta_out[args.accession_col] = meta_out[args.accession_col].astype(str).str.strip()

    merged = meta_out.merge(
        manifest,
        how="left",
        left_on=args.accession_col,
        right_on="accession"
    )

    # Keep a clean annotation_available column even if manifest missing
    if "annotation_available" not in merged.columns:
        merged["annotation_available"] = pd.NA

    # Optionally drop the extra 'accession' column from manifest merge
    merged.drop(columns=["accession"], inplace=True, errors="ignore")

    merged.to_csv(updated_meta_path, sep="\t", index=False)

    print("DOWNLOADS COMPLETE!")
    print(f"Manifest: {manifest_path}")
    print(f"Updated NEW_ONLY metadata with paths: {updated_meta_path}")
    print(f"Failed accessions log: {failed_log}")
    print(f"Processed accessions log: {processed_log}")
    print(f"Datasets stderr log: {stderr_log}")


if __name__ == "__main__":
    main()
