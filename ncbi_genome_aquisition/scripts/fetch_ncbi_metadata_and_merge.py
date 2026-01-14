#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Set, Tuple

import pandas as pd


# Default fields. I think I included all fields I care about. 
DEFAULT_FIELDS = ",".join([
    "current-accession",
    "assminfo-assembly-method",
    "assminfo-bioproject",
    "assminfo-biosample-accession",
    "assminfo-biosample-bioproject-accession",
    "assminfo-biosample-collected-by",
    "assminfo-biosample-collection-date",
    "assminfo-biosample-description-organism-name",
    "assminfo-biosample-description-organism-tax-id",
    "assminfo-biosample-description-title",
    "assminfo-biosample-geo-loc-name",
    "assminfo-biosample-host",
    "assminfo-biosample-ids-db",
    "assminfo-biosample-ids-value",
    "assminfo-biosample-last-updated",
    "assminfo-biosample-publication-date",
    "assminfo-biosample-strain",
    "assminfo-biosample-sub-species",
    "assminfo-biosample-submission-date",
    "assminfo-blast-url",
    "assminfo-level",
    "assminfo-name",
    "assminfo-refseq-category",
    "assminfo-release-date",
    "assminfo-sequencing-tech",
    "assminfo-status",
    "assminfo-type",
    "assmstats-contig-l50",
    "assmstats-contig-n50",
    "assmstats-gc-percent",
    "assmstats-genome-coverage",
    "assmstats-number-of-component-sequences",
    "assmstats-number-of-contigs",
    "assmstats-total-sequence-len",
    "organism-infraspecific-strain",
    "organism-tax-id",
    "type_material-display_text",
    "type_material-label",
])


def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def check_exec_exists(exe: str) -> None:
    from shutil import which
    if which(exe) is None:
        die(
            f"Required executable not found in PATH: '{exe}'. "
            f"Activate your conda env (datasets/dataformat) and try again."
        )


def read_taxa_file(path: Path) -> List[str]:
    """
    Accepts either:
      - CSV with a header (first column used)
      - headerless 1-column file (treated as data)
      - one-name-per-line text file
    """
    text_lines = [re.sub(r"\r$", "", l).strip() for l in path.read_text().splitlines()]
    text_lines = [l for l in text_lines if l]  # drop empty lines
    if not text_lines:
        return []

    # If it looks like a simple one-per-line list (no commas/tabs), treat as plain text
    if all(("," not in l and "\t" not in l) for l in text_lines):
        taxa = text_lines
    else:
        # Otherwise treat as delimited; assume header may or may not exist.
        # Read twice: once assuming header, once assuming no header, then pick the one with more rows.
        df_header = pd.read_csv(path, dtype=str)
        df_noheader = pd.read_csv(path, dtype=str, header=None)

        # First column values
        vals_header = df_header.iloc[:, 0].dropna().astype(str).tolist() if df_header.shape[1] >= 1 else []
        vals_noheader = df_noheader.iloc[:, 0].dropna().astype(str).tolist() if df_noheader.shape[1] >= 1 else []

        taxa = vals_noheader if len(vals_noheader) >= len(vals_header) else vals_header

    # trim + de-dupe preserving order
    seen = set()
    out = []
    for t in taxa:
        t = t.strip()
        if t and t not in seen:
            seen.add(t)
            out.append(t)
    return out


def read_accessions_file(path: Path) -> List[str]:
    accs: List[str] = []
    for line in path.read_text().splitlines():
        line = re.sub(r"\r$", "", line).strip()
        if line:
            accs.append(line)

    # de-dupe, preserve order
    seen = set()
    out: List[str] = []
    for a in accs:
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def sanitize_filename(s: str) -> str:
    s = s.strip()
    s = s.replace(" ", "_").replace("/", "_").replace(":", "_")
    s = re.sub(r"[^A-Za-z0-9_.()-]+", "", s)
    return s[:200] if len(s) > 200 else s


def run_cmd(cmd: List[str]) -> Tuple[int, str, str]:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return p.returncode, p.stdout, p.stderr


def fetch_one_taxon(taxon: str, fields: str) -> Optional[pd.DataFrame]:
    """
    Fetch metadata TSV via datasets->dataformat pipeline.
    Returns a DataFrame or None if nothing is returned.
    """
    cmd = [
        "bash", "-lc",
        f'datasets summary genome taxon "{taxon}" --mag exclude --as-json-lines | '
        f'dataformat tsv genome --fields "{fields}"'
    ]
    rc, out, err = run_cmd(cmd)
    if rc != 0:
        print(f"[WARN] Taxon fetch failed: {taxon}\n{err}", file=sys.stderr)

    out = out.strip()
    if not out:
        return None

    from io import StringIO
    df = pd.read_csv(StringIO(out), sep="\t", dtype=str)
    df["query_input"] = taxon
    df["query_type"] = "taxon"
    return df


def fetch_one_accession(acc: str, fields: str) -> Optional[pd.DataFrame]:
    cmd = [
        "bash", "-lc",
        f'datasets summary genome accession "{acc}" --as-json-lines | '
        f'dataformat tsv genome --fields "{fields}"'
    ]
    rc, out, err = run_cmd(cmd)
    if rc != 0:
        print(f"[WARN] Accession fetch failed: {acc}\n{err}", file=sys.stderr)

    out = out.strip()
    if not out:
        return None

    from io import StringIO
    df = pd.read_csv(StringIO(out), sep="\t", dtype=str)
    df["query_input"] = acc
    df["query_type"] = "accession"
    return df


def canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert display headers like:
      'Assembly BioProject Accession' -> 'assembly_bioproject_accession'
    while preserving uniqueness.
    """
    def slug(s: str) -> str:
        s = s.strip().lower()
        s = re.sub(r"[^a-z0-9]+", "_", s)
        s = re.sub(r"_+", "_", s).strip("_")
        return s

    cols = [slug(c) for c in df.columns]

    # De-duplicate if collisions occur
    seen = {}
    out_cols = []
    for c in cols:
        if c not in seen:
            seen[c] = 1
            out_cols.append(c)
        else:
            seen[c] += 1
            out_cols.append(f"{c}_{seen[c]}")

    df = df.copy()
    df.columns = out_cols
    return df


def normalize_accession_col(df: pd.DataFrame) -> pd.DataFrame:
    """
    dataformat headers differ across versions, e.g.:
      - current_accession
      - current-accession
      - Current Accession
      - current accession
    This function finds the column in a robust way and renames it to 'assembly_acc'.
    """
    def canon(s: str) -> str:
        return re.sub(r"[^a-z0-9]+", "", s.lower())

    canon_map = {canon(c): c for c in df.columns}

    # If canonicalized columns already exist (e.g. after canonicalize_columns),
    # this still works.
    candidates = [
        "currentaccession",
        "currentaccessionversion",
        # fallback options (rare)
        "assemblyaccession",
        "accession",
    ]

    found = None
    for key in candidates:
        if key in canon_map:
            found = canon_map[key]
            break

    if found is None:
        preview = list(df.columns)[:40]
        preview_canon = [canon(c) for c in preview]
        die(
            "Could not find current accession column in fetched metadata.\n"
            f"Columns seen: {preview} ...\n"
            f"Canonical preview: {preview_canon} ..."
        )

    df = df.copy()
    df.rename(columns={found: "assembly_acc"}, inplace=True)
    return df


def load_master_metadata(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in [".tsv", ".tab"]:
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, dtype=str)


def build_known_accessions(master: pd.DataFrame,
                          col_primary: str,
                          col_corresponding: str) -> Set[str]:
    if col_primary not in master.columns:
        die(f"Master metadata missing required column: {col_primary}")
    if col_corresponding not in master.columns:
        die(f"Master metadata missing required column: {col_corresponding}")

    s = pd.Series(
        pd.concat([master[col_primary], master[col_corresponding]], ignore_index=True),
        dtype="string"
    )
    s = s.dropna().astype(str).str.strip()
    s = s[s != ""]
    return set(s.unique().tolist())


def align_and_append(master: pd.DataFrame, new_rows: pd.DataFrame) -> pd.DataFrame:
    """
    Append new_rows to master by union of columns.
    Keeps master columns intact; adds any new columns from fetched metadata.
    """
    master = master.copy()
    new_rows = new_rows.copy()

    for c in master.columns:
        if c not in new_rows.columns:
            new_rows[c] = pd.NA
    for c in new_rows.columns:
        if c not in master.columns:
            master[c] = pd.NA

    # reorder new_rows to master col order
    new_rows = new_rows[master.columns.tolist()]

    combined = pd.concat([master, new_rows], ignore_index=True)

    # If assembly_acc exists, enforce uniqueness
    if "assembly_acc" in combined.columns:
        combined["assembly_acc"] = combined["assembly_acc"].astype("string").str.strip()
        combined = combined.drop_duplicates(subset=["assembly_acc"], keep="first")

    return combined


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Fetch NCBI genome metadata (datasets+dataformat) from taxa OR accessions, "
            "optionally deduplicate against an existing master metadata sheet, and write outputs."
        )
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--taxa_file", type=str, help="CSV or TXT with taxa names (first column used for CSV).")
    group.add_argument("--accessions_file", type=str, help="TXT with one NCBI assembly accession per line.")

    parser.add_argument("--master_metadata", type=str, default=None,
                        help="Existing metadata sheet to deduplicate against. If provided, writes NEW-only and UPDATED master.")
    parser.add_argument("--master_primary_col", type=str, default="assembly_acc",
                        help="Column in master metadata containing NCBI assembly accessions already in DB.")
    parser.add_argument("--master_corresponding_col", type=str, default="corresponding_ncbi_accession",
                        help="Column in master metadata for NCBI accessions corresponding to JGI genomes.")
    parser.add_argument("--outdir", type=str, required=True, help="Output directory.")
    parser.add_argument("--fields", type=str, default=DEFAULT_FIELDS, help="Comma-separated dataformat fields.")
    parser.add_argument("--prefix", type=str, default=None,
                        help="Prefix for output files (default derived from input filename).")
    parser.add_argument("--write_all_fetched", action="store_true",
                        help="Also write the full fetched metadata before deduplication (useful for auditing).")
    parser.add_argument("--canonicalize_headers", action="store_true",
                        help="Convert fetched headers to snake_case for stability across dataformat versions.")

    args = parser.parse_args()

    check_exec_exists("datasets")
    check_exec_exists("dataformat")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    if args.prefix:
        prefix = sanitize_filename(args.prefix)
    else:
        in_path = Path(args.taxa_file or args.accessions_file)
        prefix = sanitize_filename(in_path.stem)

    # Read inputs
    dfs: List[pd.DataFrame] = []
    if args.taxa_file:
        taxa = read_taxa_file(Path(args.taxa_file))
        if not taxa:
            die("No taxa found in taxa_file.")
        for t in taxa:
            df = fetch_one_taxon(t, args.fields)
            if df is not None and not df.empty:
                dfs.append(df)
        mode = "taxa"
    else:
        accs = read_accessions_file(Path(args.accessions_file))
        if not accs:
            die("No accessions found in accessions_file.")
        for a in accs:
            df = fetch_one_accession(a, args.fields)
            if df is not None and not df.empty:
                dfs.append(df)
        mode = "accessions"

    if not dfs:
        die("No metadata returned from NCBI for the provided inputs.")

    fetched = pd.concat(dfs, ignore_index=True)

    # Optionally canonicalize *all* headers so downstream merges are stable
    if args.canonicalize_headers:
        fetched = canonicalize_columns(fetched)

    # Normalize accession column to assembly_acc
    fetched = normalize_accession_col(fetched)

    # Add a standard source label
    if "source" not in fetched.columns:
        fetched["source"] = "ncbi"
    else:
        fetched["source"] = fetched["source"].fillna("ncbi")

    # Basic cleanup
    fetched["assembly_acc"] = fetched["assembly_acc"].astype("string").str.strip()
    fetched = fetched[fetched["assembly_acc"].notna() & (fetched["assembly_acc"] != "")]
    fetched = fetched.drop_duplicates(subset=["assembly_acc"], keep="first")

    # Optional audit output
    if args.write_all_fetched:
        all_fetched_path = outdir / f"{prefix}.{mode}.ALL_FETCHED.tsv"
        fetched.to_csv(all_fetched_path, sep="\t", index=False)

    # If no master metadata provided, just write fetched and exit
    if not args.master_metadata:
        fetched_path = outdir / f"{prefix}.{mode}.FETCHED.tsv"
        fetched.to_csv(fetched_path, sep="\t", index=False)
        print(f"Wrote fetched metadata (no dedup requested): {fetched_path}")
        return

    # Deduplicate against master
    master_path = Path(args.master_metadata).resolve()
    master = load_master_metadata(master_path)

    known = build_known_accessions(
        master,
        col_primary=args.master_primary_col,
        col_corresponding=args.master_corresponding_col
    )

    new_only = fetched[~fetched["assembly_acc"].isin(known)].copy()

    # Write new-only metadata
    new_only_path = outdir / f"{prefix}.{mode}.NEW_ONLY.tsv"
    new_only.to_csv(new_only_path, sep="\t", index=False)

    # Write a clean accession list for downstream downloads
    acc_list_path = outdir / f"{prefix}.{mode}.NEW_ONLY_accessions.txt"
    accs = new_only["assembly_acc"].dropna().astype(str).str.strip()
    accs = accs[accs != ""].drop_duplicates()
    acc_list_path.write_text("\n".join(accs) + "\n")
    print(f"Wrote NEW-only accession list: {acc_list_path}")

    # Append to master and write updated master
    updated = align_and_append(master, new_only)
    updated_path = outdir / f"{prefix}.{mode}.MASTER_UPDATED.csv"
    updated.to_csv(updated_path, index=False)

    print(f"Wrote NEW-only metadata: {new_only_path}")
    print(f"Wrote UPDATED master metadata: {updated_path}")
    print(f"New genomes to add: {new_only.shape[0]} (out of {fetched.shape[0]} fetched)")


if __name__ == "__main__":
    main()
