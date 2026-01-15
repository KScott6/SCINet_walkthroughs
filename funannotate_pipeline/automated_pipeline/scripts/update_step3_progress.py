import os
import pandas as pd
from datetime import datetime

# === CONFIGURATION ===
project_dir = "/project/arsef/projects/bulk_genome_annotation"
working_dir = os.path.join(project_dir, "needs_annotation/working_files")
progress_file = os.path.join(project_dir, "progress/annotation_master_progress.tsv")
step_col = "step3_done"

# Load progress tracker
df = pd.read_csv(progress_file, sep="\t", dtype=str).fillna("")
updated = 0

for i, row in df.iterrows():
    ome = row["OMEcode"]
    busco_dir = os.path.join(working_dir, ome, "busco", f"{ome}_prelim")

    if os.path.exists(busco_dir):
        mod_time = datetime.fromtimestamp(os.path.getmtime(busco_dir)).strftime("%Y-%m-%d %H:%M:%S")
        if df.at[i, step_col] != mod_time:
            df.at[i, step_col] = mod_time
            updated += 1
    elif row["step3_job"] != "FAILED" and row[step_col] == "":
        print(f"[ ] Incomplete: {ome} — BUSCO output not found")

print(f"[✓] Updated {updated} genomes in {step_col}.")

# Save updated progress file
df.to_csv(progress_file, sep="\t", index=False)
print(f"[✓] Saved: {progress_file}")
