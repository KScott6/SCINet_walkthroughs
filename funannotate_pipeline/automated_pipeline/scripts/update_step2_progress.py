import os
import pandas as pd
from datetime import datetime

# === CONFIGURATION ===
project_dir = "/project/arsef/projects/bulk_genome_annotation"
working_dir = os.path.join(project_dir, "needs_annotation/working_files")
progress_file = os.path.join(project_dir, "progress/annotation_master_progress.tsv")

# Load progress file
df = pd.read_csv(progress_file, sep="\t", dtype=str).fillna("")

updated = 0

for i, row in df.iterrows():
    ome = row["OMEcode"]
    prep_dir = os.path.join(working_dir, ome, "prep")

    masked = os.path.join(prep_dir, f"{ome}.masked.fasta")
    fa = os.path.join(prep_dir, f"{ome}-families.fa")
    stk = os.path.join(prep_dir, f"{ome}-families.stk")

    if all(os.path.exists(path) for path in [masked, fa, stk]):
        mod_time = datetime.fromtimestamp(os.path.getmtime(masked)).strftime("%Y-%m-%d %H:%M:%S")
        if df.at[i, "step2_done"] != mod_time:
            df.at[i, "step2_done"] = mod_time
            updated += 1
    elif row["step2_job"] != "FAILED" and row["step2_done"] == "":
        print(f"[ ] Incomplete: {ome} (missing output files)")

print(f"[✓] Updated step2_done for {updated} genomes.")

# Save progress file
df.to_csv(progress_file, sep="\t", index=False)
print(f"[✓] Saved updated progress file: {progress_file}")
