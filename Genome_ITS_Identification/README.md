# BLAST-based ID of genomes

Here is a quick way of checking the identity or screening for contamination in your genomes using BLAST. This example uses ITS/"Internal Transcribed Spacer", but you can use whatever other region you want (TEF1, COX1, etc).

Your genomes of interest should all be in a single folder and named with a similar structure (same suffix, file extension).

The first script will create individual job scripts to BLAST a reference nucelotide sequence against your genomes:

```bash
python /project/arsef/scripts/genome_scripts/make_its_jobs.py \
  --ome-list /project/arsef/projects/ambrosiella/final_versions_ncbi/its_ID/ome_list.txt \
  --sample_column none \
  --asm-root /project/arsef/projects/ambrosiella/final_versions_ncbi/ncbi_assemblies \
  --asm-pattern "{ome}.fasta" \
  --its-out-root /project/arsef/projects/ambrosiella/final_versions_ncbi/its_ID/its_extraction \
  --cmd-root /project/arsef/projects/ambrosiella/final_versions_ncbi/its_ID/scripts \
  --log-root /project/arsef/projects/ambrosiella/final_versions_ncbi/its_ID/logs \
  --results-root /project/arsef/projects/ambrosiella/final_versions_ncbi/its_ID/final_results \
  --region-name its \
  --region-ref /project/arsef/projects/ambrosiella/final_versions_ncbi/its_ID/Ambrosiella_ITS_ref.fasta \
  --nt-db /project/arsef/databases/nt/nt \
  --genome-hit-count all \

  --submit
```

`--ome-list` single-column list or metadata sheet containing a column of all your ome codes. 

`--ome_column` specify the column name containing your ome codes. If no header, specify "--sample_column none"

`--asm-root` path where all your genome assemblies are stored

`--asm-pattern ` specify the structure of your genomes' file names. Examples:  --asm-pattern "{ome}_scaffolds.fasta" ; --asm-pattern "{ome}.fna"

`--its-out-root` where the extracted genome sequences should be placed

`--cmd-root` path where you want the job scripts to go

`--log-root` path where you want your slurm logs to go

`--results-root` path where you want the final results to go

`--region-name` name of the region you are investigating. Used to name job script, slurm output paths, reports, etc.

`--region-ref` path to a reference nucelotide sequence, used to BLAST against your genomes dtabases.

`--nt-db` path to your NCBI nt database

`--genome-hit-count` (integer or "all") How many genome-extracted sequences you want to blast against the NCBI nt database. Recommended that you choose "all", especially if you're screening for fungal contamination. 

`--submit` include this flag to automatically submit jobs

<br>

After all your ITS jobs are complete, compile all your results (per genome) using the generated command:

For me, it was: 

```bash
bash /project/arsef/projects/collab/costa_rica/its_id/scripts/combine_its_nt_results.sh
```