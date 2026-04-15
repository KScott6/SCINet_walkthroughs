
# Genome assembly from PacBio Hifi data

This pipeline assumes you are a SCINet user and are in the arsef project. This pipeline also assumes you have demultiplexed PacBio reads, one .bam/.bam.bpi file set per sample.

This pipeline takes raw HiFi reads from one or more samples and performs:

1. BAM → FASTQ conversion
2. Read QC/filtering (fastplong)
3. Genome assembly (hifiasm / flye)
4. (Optional downstream steps: QC, annotation)

All steps are driven by a single metadata file (TSV/CSV). 

<br>

**Required columns in the metadata file**

sample_name
   * Used for all output naming
   * Must not contain spaces or special characters (use _ or - only)

barcode 
   * Used to locate BAM files in bam2fastq step

genome_size_Mbp
   * Used in assembly step

ncbi_tax_ID 
   * Used in NCBI FCS genome cleaning step
   * get the TaxIDs from NCBI Taxonomy database
   * Be as specific as you're able to be (to genus level preferred)

The species column in the example sheet is not used by any of my scripts and only exists for your benefit. 


## Step 1 : BAM → FASTQ conversion (pbtk / bam2fastq v3.5.0)

Convert PacBio BAM files into FASTQ format and renames using the provided sample_name. The output fastq reads will be stored in a subfolder labeled after the sample_name.

Script : make_bam2fastq_jobs.py

Inputs :
* metadata file
* directory of BAM files (barcode.bam)
  
Output :
* renamed FASTQ files


```bash
python make_bam2fastq_jobs.py \
  --metadata sample_metadata.tsv \
  --bam_dir [path to folder containing bam files] \
  --script_dir [path to desired job scripts folder] \
  --log_dir [path to desired slurm log folder] \
  --fastq_dir [path to desired parent fastq read output folder] \
  --threads 24 \
  --submit
```

Example useage:

```bash
python /project/arsef/scripts/genome_assembly_scripts/pacbio/make_bam2fastq.py \
  --sample_sheet /project/arsef/projects/genome_assembly/April_PacBio/barcodes_isolates.txt \
  --bam_dir /project/arsef/projects/genome_assembly/April_PacBio/raw_reads \
  --script_dir /project/arsef/projects/genome_assembly/April_PacBio/commands/bam2fastq \
  --log_dir /project/arsef/projects/genome_assembly/April_PacBio/logs/bam2fastq \
  --fastq_dir /project/arsef/projects/genome_assembly/April_PacBio/raw_fastq \
  --threads 24 \
  --time 24:00:00 \
  --account arsef \
  --conda_env /project/arsef/environments/pacbio \
  --command_name "bam2fastq" \
  --submit
```

Explanation of key options:

`--metadata`	Input metadata table with sample names and barcodes

`--bam_dir`	Location of input BAM files (requires all .bam and .bam.pbi files to be in the same folder)

`--fastq_dir`	Desired output parent FASTQ directory

`--threads`	Threads for bam2fastq

`--submit`	Submit jobs immediately

`--overwrite`	Re-run even if output exists


Notes:

* BAM files must match barcode exactly
* Output naming is based on sample_name
* Script will:
* will skip missing BAM files
* will skip samples if FASTQ already exists (unless --overwrite)
* Log directory is created automatically before submission

<br>

---

## Step 2 : Read QC and filtering (fastplong v0.4.1)

PacBio data is usually great quality, but it's still worth running it through fastplong to filter out any missed adapter sequences. Also creates a nice report of your read data.

Script : make_fastplong_jobs.py

Inputs
* FASTQ files from Step 1
* desired output folder

Outputs:
* filtered reads (quality and adapter content)
* one html/json report per sample

```bash
python make_fastplong_jobs.py \
  --metadata sample_metadata.tsv \
  --input_dir [path to folder containing fastq files]  \
  --script_dir [path to desired job scripts folder] \
  --log_dir [path to desired slurm log folder] \
  --output_dir [path to desired parent filtered read and report output folder]  \
  --threads 16 \
  --submit
```

Example usage:

```bash
python /project/arsef/scripts/genome_assembly_scripts/pacbio/make_fastplong.py \
  --metadata /project/arsef/projects/genome_assembly/April_PacBio/sample_data.txt \
  --sample_column sample_name \
  --input_dir /project/arsef/projects/genome_assembly/April_PacBio/raw_fastq/ \
  --script_dir /project/arsef/projects/genome_assembly/April_PacBio/commands/fastplong/ \
  --log_dir /project/arsef/projects/genome_assembly/April_PacBio/logs/fastplong/ \
  --output_dir /project/arsef/projects/genome_assembly/April_PacBio/fastplong/ \
  --threads 24 \
  --time 24:00:00 \
  --account arsef \
  --conda_env /project/arsef/environments/pacbio \
  --submit
```

Explanation of key options:

`--metadata`	Metadata file

`--sample_column`	Column names containing sample naming info (default: sample_name)

`--input_dir`	FASTQ input location

`--output_dir`	QC output location

`--report_suffix`	Suffix for HTML/JSON reports; prefix will be sample_name

`--threads`	Threads (max 16; fastplong will not use more than this)

`--overwrite`	Force re-run

<br>

---

## Genome assembly

Two assembly software options are available in this pipeline:  Flye and hifiasm. Run both if you are unsure of which to use, then compare. I have found that some fungal taxa work better with one over the other. 

Both of these software require that you report an estimated genome size (in Mbp) in the metadata file. The provided size does not need to be perfect, just a ballpark.

Script: make_assembly_jobs.py

Inputs:
* metadata file (tsv/csv) with at least sample_name and genome_size_Mbp info
* input read directory (from either fastplong or raw_fastq)
* optional: prefix file with subset of prefixes to run

Outputs:
* genome assemblies and all associated files
* renamed final assemblies to specificed output folder

This script can detect input reads in both a flat layout (input_dir/OME.filt.fastq.gz) or nested (input_dir/OME/OME.filt.fastq.gz).

If you followed this pipeline exactly, you want to specify "--input_suffix .filt.fastq.gz" so the script knows to use the high-quality fitered reads (.filt.fastq.gz), and not the filtered out reads (.failed.fastq.gz). 

<br>

Explanation of key options:

`--metadata`	Metadata table (required)

`--assembly_software`	flye or hifiasm

`--flye_scaffold` Turn on scaffolding option (Flye only)

`--assembly_suffix_append` [string]  Include a character string to add to your genome assemblies.

`--input_dir`	Location of FASTQ reads

`--input_suffix`	Recommended: .filt.fastq.gz (default highest priority)

`--script_dir`	Output location for job scripts

`--log_dir`	SLURM logs (auto-created)

`--assembly_dir`	Per-sample working directories

`--final_assembly_dir`	Centralized final FASTA output

`--threads`	CPU threads

`--mem_per_cpu`	Memory per CPU

`--submit`	Submit jobs immediately

`--overwrite	`Re-run even if output exists

<br>

### Option A : Flye (v2.9.2-b1786)

Flye has the option to automatically scaffold your genome -- to turn on this option, include tthe paramter --flye_scaffold in your command. I also recommend you add a custom suffix to your genome, so that non-scaffolded flye genome assemblies do not get overwritten/confused (--assembly_suffix_append).

```bash
python make_assembly_jobs.py \
  --metadata sample_metadata.tsv \
  --assembly_software flye \
  --input_dir fastplong \
  --input_suffix .filt.fastq.gz \
  --script_dir commands/flye \
  --log_dir logs/flye \
  --assembly_dir flye \
  --final_assembly_dir final_assemblies \
  --threads 24 \
  --submit
```

Example usage:

```bash
python /project/arsef/scripts/genome_assembly_scripts/pacbio/make_assembly_jobs.py \
  --metadata /project/arsef/projects/genome_assembly/April_PacBio/sample_data.txt \
  --assembly_software flye \
  --input_dir /project/arsef/projects/genome_assembly/April_PacBio/fastplong \
  --input_suffix .filt.fastq.gz \
  --script_dir /project/arsef/projects/genome_assembly/April_PacBio/commands/flye \
  --log_dir /project/arsef/projects/genome_assembly/April_PacBio/logs/flye \
  --assembly_dir /project/arsef/projects/genome_assembly/April_PacBio/flye \
  --final_assembly_dir /project/arsef/projects/genome_assembly/April_PacBio/final_assemblies \
  --input_suffix .filt.fastq.gz \
  --threads 24 \
  --time 72:00:00 \
  --mem_per_cpu 8000 \
  --account arsef \
  --conda_env /project/arsef/environments/pacbio \
  --submit
```

Example with scaffolding option:

```bash
python /project/arsef/scripts/genome_assembly_scripts/pacbio/make_assembly_jobs.py \
  --metadata /project/arsef/projects/genome_assembly/April_PacBio/sample_data.txt \
  --assembly_software flye \
  --input_dir /project/arsef/projects/genome_assembly/April_PacBio/fastplong \
  --input_suffix .filt.fastq.gz \
  --script_dir /project/arsef/projects/genome_assembly/April_PacBio/commands/flye2 \
  --log_dir /project/arsef/projects/genome_assembly/April_PacBio/logs/flye2 \
  --assembly_dir /project/arsef/projects/genome_assembly/April_PacBio/flye2 \
  --final_assembly_dir /project/arsef/projects/genome_assembly/April_PacBio/final_assemblies \
  --threads 24 \
  --mem_per_cpu 8000 \
  --flye_scaffold \
  --assembly_suffix_append _scaffold \
  --submit
```

<br>

### Option B : Hifiasm (v0.21.0-r686)

```bash
python make_assembly_jobs.py \
  --metadata sample_metadata.tsv \
  --assembly_software hifiasm \
  --input_dir fastplong \
  --input_suffix .filt.fastq.gz \
  --script_dir commands/hifiasm \
  --log_dir logs/hifiasm \
  --assembly_dir hifiasm \
  --final_assembly_dir final_assemblies \
  --threads 24 \
  --submit
```

Example usage:

```bash
python /project/arsef/scripts/genome_assembly_scripts/pacbio/make_assembly_jobs.py \
  --metadata /project/arsef/projects/genome_assembly/April_PacBio/sample_data.txt \
  --assembly_software hifiasm \
  --input_dir /project/arsef/projects/genome_assembly/April_PacBio/fastplong \
  --input_suffix .filt.fastq.gz \
  --script_dir /project/arsef/projects/genome_assembly/April_PacBio/commands/hifiasm \
  --log_dir /project/arsef/projects/genome_assembly/April_PacBio/logs/hifiasm \
  --assembly_dir /project/arsef/projects/genome_assembly/April_PacBio/hifiasm \
  --final_assembly_dir /project/arsef/projects/genome_assembly/April_PacBio/final_assemblies \
  --input_suffix .filt.fastq.gz \
  --threads 24 \
  --time 72:00:00 \
  --mem_per_cpu 8000 \
  --account arsef \
  --conda_env /project/arsef/environments/pacbio \
  --submit
```

<br>

---

## Assessing basic genome stats with QUAST (v5.2.0)

After assemblies are complete, QUAST can be run on all final genome assemblies together to generate comparative assembly statistics.

This step creates **one QUAST job** that includes all assembly FASTA files found in the specified parent directory.

Script: make_quast_job.py

Input:

* parent directory containing the final assembly fasta files
* metadata file (optional)

Output:

* QUAST output files

```bash
python make_quast_job.py \
  --assembly_parent_dir final_assemblies \
  --script_dir commands/quast \
  --log_dir logs/quast \
  --quast_base_dir quast \
  --threads 16 \
  --submit
```

Example usage:

```bash
python /project/arsef/scripts/genome_assembly_scripts/pacbio/make_quast_job.py \
  --assembly_parent_dir /project/arsef/projects/genome_assembly/April_PacBio/final_assemblies \
  --script_dir /project/arsef/projects/genome_assembly/April_PacBio/commands/quast \
  --log_dir /project/arsef/projects/genome_assembly/April_PacBio/logs/quast \
  --quast_base_dir /project/arsef/projects/genome_assembly/April_PacBio/quast \
  --run_name run1 \
  --threads 16 \
  --submit
```

Explanation of key options:

`--assembly_parent_dir`	Folder containing final assemblies

`--script_dir`	Output location for QUAST job script

`--log_dir`	SLURM log directory

`--quast_base_dir`	Base directory for QUAST output folders

`--run_name`	Optional custom name for this QUAST run (if no name provided, will be labeled with a timestamp)

`--threads`	Number of CPU threads

`--metadata` If provided, the script will only include the genomes with sample names present in the metadata sheet.

`--submit`	Submit job immediately

`--overwrite`	Allow reuse of an existing output folder

<br>

---

## Assessment of genome completeness with BUSCO (v5.8.3)

After assemblies are generated, BUSCO can be run on each genome to assess completeness using a lineage-specific BUSCO database.

This step creates **one BUSCO job per assembly** and copies each BUSCO short summary into a customizable centralized summary folder for easy comparison.

Script: make_busco_jobs.py

Input:

* parent directory containing the final assembly fasta files
* metadata file (optional)

Output:

* BUSCO output files (including short_summary file)

```bash
python make_busco_jobs.py \
  --assembly_parent_dir final_assemblies \
  --script_dir commands/busco \
  --log_dir logs/busco \
  --busco_out_dir busco \
  --short_summary_dir busco/__busco_short_summaries \
  --busco_lineage /project/arsef/databases/funannotate_databases/fungi_odb10 \
  --threads 16 \
  --submit
```

Example useage:

```bash
python /project/arsef/scripts/genome_assembly_scripts/pacbio/make_busco_jobs.py \
  --assembly_parent_dir /project/arsef/projects/genome_assembly/April_PacBio/final_assemblies \
  --metadata /project/arsef/projects/genome_assembly/April_PacBio/sample_data.txt \
  --sample_column sample_name \
  --script_dir /project/arsef/projects/genome_assembly/April_PacBio/commands/busco \
  --log_dir /project/arsef/projects/genome_assembly/April_PacBio/logs/busco \
  --busco_out_dir /project/arsef/projects/genome_assembly/April_PacBio/busco \
  --short_summary_dir /project/arsef/projects/genome_assembly/April_PacBio/busco/__busco_short_summaries \
  --busco_lineage /project/arsef/databases/funannotate_databases/fungi_odb10 \
  --threads 16 \
  --time 24:00:00 \
  --account arsef \
  --conda_env /project/arsef/environments/funannotate \
  --force_busco \
  --submit
```

Explanation of key options:

`--assembly_parent_dir`	Folder containing final assemblies

`--busco_lineage`	BUSCO lineage/database path

`--script_dir`	Output location for BUSCO job scripts

`--log_dir`	SLURM log directory

`--busco_out_dir`	Parent BUSCO output directory

`--short_summary_dir`	Folder for copied/renamed BUSCO summary files

`--threads`	CPU threads

`--submit`	Submit jobs immediately

`--overwrite`	Re-create script for sample, even if BUSCO summary already exists

`--force_busco`	Add -f (force overwrite) to BUSCO command

<br>
<br>

---


## Screening and cleaning genomes for contaminants with NCBI FCS 

After assembly and quality assessment, genomes can be screened for contamination using the NCBI Foreign Contamination Screen (FCS-GX). This software is automatically run on any genome you try to submit to NCBI, so doing this yourself will save a bit of time and effort.

This workflow consists of **two steps**:

1. **Screening** — identifies potentially contaminant contigs/regions  
2. **Cleaning** — removes or splits flagged contaminant regions based on the screening report  

FCS-GX is run in two phases because screening is computationally intensive (>=600GB) and benefits from batching multiple genomes per job, while cleaning is lightweight and best handled with one genome per job.

<br>

### Genome screening

Script: make_fcs_screen_jobs.py

This step depends on the proper NCBI taxID metadata provided in the sample metadata sheet. The more specific you can be with the TaxID, the better.

```bash
python make_fcs_screen_jobs.py \
  --metadata sample_data.txt \
  --sample-col sample_name \
  --taxid-col ncbi_tax_ID \
  --assembly-parent-dir final_assemblies \
  --assembly-software hifiasm \
  --gxdb /project/arsef/databases/gxdb \
  --assemblies-per-job 8 \
  --threads 32 \
  --mem 600G \
  --out-root fcs_gx/screen \
  --cmd-root commands/fcs_screen \
  --log-root logs/fcs_screen \
  --submit
```

Example command:

```bash
python /project/arsef/scripts/genome_assembly_scripts/pacbio/make_fcs_screen_jobs.py \
  --metadata /project/arsef/projects/genome_assembly/April_PacBio/sample_data.txt \
  --sample-col sample_name \
  --taxid-col ncbi_tax_ID \
  --assembly-parent-dir /project/arsef/projects/genome_assembly/April_PacBio/final_assemblies \
  --gxdb /project/arsef/databases/gxdb \
  --assemblies-per-job 8 \
  --threads 32 \
  --time 48:00:00 \
  --mem 600G \
  --out-root /project/arsef/projects/genome_assembly/April_PacBio/fcs_gx/screen \
  --cmd-root /project/arsef/projects/genome_assembly/April_PacBio/commands/fcs_screen \
  --log-root /project/arsef/projects/genome_assembly/April_PacBio/logs/fcs_screen \
  --submit
```

Explanation of key options:

`--assembly-parent-dir`	Folder containing assemblies to screen

`--assembly-software`	Optional filter for flye or hifiasm assemblies (if you don't specify either, all will be run)

`--gxdb`	Path to FCS GX database

`--assemblies-per-job`	Number of genomes screened sequentially per job

`--taxid-col`	Metadata column name containing NCBI taxids

`--submit`	Submit generated jobs immediately

<br>

### Genome cleaning

This step applies the FCS action report to produce cleaned assemblies.

Two cleaning modes are supported: 

`exclude`	Remove contaminant contigs/regions entirely

`split`	Split flagged contigs instead of removing them

(I recommend using `exclude`, espeically if you're trying to prepare for NCBI submission)

```bash
python make_fcs_clean_jobs.py \
  --assembly-parent-dir final_assemblies \
  --screen-out-root fcs_gx/screen \
  --final-clean-dir final_assemblies/fcs_cleaned_exclude \
  --metadata sample_data.txt \
  --sample-col sample_name \
  --assembly-software hifiasm \
  --mode exclude \
  --out-root fcs_gx/clean \
  --cmd-root commands/fcs_clean \
  --log-root logs/fcs_clean \
  --submit
```


Example command:

```bash
python /project/arsef/scripts/genome_assembly_scripts/pacbio/make_fcs_clean_jobs.py \
  --assembly-parent-dir /project/arsef/projects/genome_assembly/April_PacBio/final_assemblies \
  --screen-out-root /project/arsef/projects/genome_assembly/April_PacBio/fcs_gx/screen \
  --final-clean-dir /project/arsef/projects/genome_assembly/April_PacBio/final_assemblies/fcs_cleaned_exclude \
  --metadata /project/arsef/projects/genome_assembly/April_PacBio/sample_data.txt \
  --sample-col sample_name \
  --assembly-software hifiasm \
  --mode exclude \
  --out-root /project/arsef/projects/genome_assembly/April_PacBio/fcs_gx/clean \
  --cmd-root /project/arsef/projects/genome_assembly/April_PacBio/commands/fcs_clean \
  --log-root /project/arsef/projects/genome_assembly/April_PacBio/logs/fcs_clean \
  --threads 1 \
  --mem 8G \
  --submit
```

Explanation of key options:

`--screen-out-root`	Folder containing FCS screening outputs

`--final-clean-dir`	Folder for final cleaned FASTA files

`--mode`	exclude or split

`--overwrite`	Re-run cleaning even if cleaned output exists


Please consider:

* FCS can occasionally over-call contamination in highly repetitive or unusual genomes
  
* Always review reports before accepting cleaned assemblies blindly
  
* Retain original uncleaned assemblies for comparison
  
* For publication-quality assemblies, manual inspection of flagged regions is recommended

* This process is entirely dependent on the taxonomy data **you** provide