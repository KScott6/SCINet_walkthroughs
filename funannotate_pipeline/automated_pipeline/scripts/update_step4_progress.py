import os
import pandas as pd
from datetime import datetime

# === CONFIGURATION ===
project_dir = "/project/arsef/projects/bulk_genome_annotation"
progress_file = os.path.join(project_dir, "progress/annotation_master_progress.tsv")
gff_dir = os.path.join(project_dir, "needs_annotation/final_funannotate_results/gff3")
step_col = "step4_done"

# Load progress tracker
df = pd.read_csv(progress_file, sep="\t", dtype=str).fillna("")
updated = 0

# Iterate over all entries
for i, row in df.iterrows():
    ome = row["OMEcode"]
    gff_path = os.path.join(gff_dir, f"{ome}_run_1.gff3")

    if os.path.exists(gff_path):
        mod_time = datetime.fromtimestamp(os.path.getmtime(gff_path)).strftime("%Y-%m-%d %H:%M:%S")
        if df.at[i, step_col] != mod_time:
            df.at[i, step_col] = mod_time
            updated += 1
    elif row["step4_job"] != "FAILED" and row[step_col] == "":
        print(f"[ ] Incomplete: {ome} — No GFF3 output found")

# Save updates
print(f"[✓] Updated {updated} genomes in {step_col}.")
df.to_csv(progress_file, sep="\t", index=False)
print(f"[✓] Saved: {progress_file}")
