import os
import pandas as pd
from datetime import datetime

# === CONFIGURATION ===
project_dir = "/project/arsef/projects/bulk_genome_annotation"
working_dir = os.path.join(project_dir, "needs_annotation/working_files")
progress_file = os.path.join(project_dir, "progress/annotation_master_progress.tsv")

# Load progress sheet
df = pd.read_csv(progress_file, sep="\t", dtype=str).fillna("")

updated = 0

for i, row in df.iterrows():
    ome = row["OMEcode"]
    sort_path = os.path.join(working_dir, ome, "prep", f"{ome}.sort.fasta")

    if os.path.exists(sort_path):
        # Get file modified time
        mod_time = datetime.fromtimestamp(os.path.getmtime(sort_path)).strftime("%Y-%m-%d %H:%M:%S")
        if df.at[i, "step1_done"] != mod_time:
            df.at[i, "step1_done"] = mod_time
            updated += 1

print(f"[✓] Updated step1_done for {updated} genomes.")

# Save progress file
df.to_csv(progress_file, sep="\t", index=False)
print(f"[✓] Saved updated progress file: {progress_file}")
