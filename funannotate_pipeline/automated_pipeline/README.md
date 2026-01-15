# Automated Funannotate Pipeline (SCINet / ARSEF)

This walkthrough describes how to run the Funannotate pipeline for one or more fungal genome annotations on SCINet. The input of this pipeline should be *unrestricted* publicly available genomes from NCBI or JGI, or your own genome assemblies which have *passed the NCBI foreign contamination screening*. The output of this pipeline will be the structural annotation (.gff) for each of your input genomes. 

This walkthrough will get you one step closer to incorporating genomes assemblies/annotations into our lab's shared MycoTools database.

In short, this pipeline will automatically:

1) standardize input assemblies

2) softmask repeats de novo (RepeatModeler + funannotate mask)

3) (completely optional) run BUSCO training for AUGUSTUS

4) run "funannotate predict" with genus-specific evidence + UniProt

5) track progress and allow safe restart/resubmission using a master progress sheet


The scripts associated with this pipeline have been uploaded to GitHub, but they are already freely available to SCINet users in the arsef project. 

<br>

Internal notes (planned updates):

* upload the automated driver and monitor script for large (>100) bulk genome anntotations.
* incorporate streamlined commands which handle the entire pipeline automatically
* incorporate automated addition of genus names associated with accessions at step 1 (modify step1 script, update)
* fix genome update count when using update scripts (step2 esp)

---

## Read before starting

### SCINet details

You will submit all python commands from the SCINet login node. They will automatically generate job scripts that will be automatically submitted to the SCINet compute nodes. 

This pipeline assumes that you have the miniconda module loaded before running each of the commands:

```bash
module load miniconda
```

This pipeline also assumes you have access to the shared Funannotate conda environment (/project/arsef/environments/funannotate_working/) and all Funannotate accessory software.

### Master progress tracker 

This file should already exist (/project/arsef/projects/bulk_genome_annotation/progress/annotation_master_progress.tsv). But if the scripts cannot find this file, a new one will be generated. 

This file is central to the pipeline - it tracks where each genome is in the pipeline and each script uses it as a reference when determing which genomes are ready for the next step. SCINet job IDs are stored in this progress tracker, as well as the time/dates when each step was verified as being successfully completed.

Each genome is referred to by its "ome code" (OMEcode). The ome code is a SIMPLE (no spaces or special characters besides . or _) prefix, which is used to label each of the intermediate and final output files generated from this pipeline. I recommend using either the NCBI accession number or the unique isolate name as each genome's ome code.

There is a column to record the genus of each genome (genus). This column is referenced during the "funannotate predict" step to determine if we have any extra protein/transcript evidence available for that particular genus. 

---

## Step 0 - Prepare input genomes

1) Input genomes needs to be in one folder. 

If you are running this pipeline immediately after the [NCBI genome acquisition pipeline](https://github.com/KScott6/SCINet_walkthroughs/tree/main/ncbi_genome_aquisition), you can simply provide the fna folder path. For example:

`/project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/ncbi_downloads/fna`

2) Rename genome input files.
   
The file names need to be in this exact format:

`<OMECode>.fna`

If your assemblies have suffixes like _cleanedassembly or the wrong extension, rename them first. This is how I rename all the .fasta files to .fna in a particular folder:

```
cd /project/arsef/projects/bulk_genome_annotation/needs_annotation/fusariums

for f in *.fasta; do
  mv "$f" "${f%.fasta}.fna"
done
```

3) Make sure input genomes are in their final form.

If you are adding genomes of your own creation (rather than downloaded from NCBI), I highly recommend that they are in their FINAL form. Submit them to NCBI so they confirm the genomes do not have contamination, or process them with [NCBI FCS](https://github.com/ncbi/fcs) yourself. 

4) Create your input list.

Create a single-column text file listing all the ome codes of the genomes you want to submit to the pipeline.

This is easy to make if you already have nicely formatted genome assembly names in a folder:

```
cd [input genome folder]
ls -1 | sed 's/\.[^.]*$//' > input_omecodes.txt
```

<br>

## Step 1 - Funannotate sort

This step runs "funannotate sort" for each genome, which sorts the contigs in the genome assembly by length (bp) and removes contigs below a specified length. This step also renames each of the contigs in an assembly and standardizes the fasta format.

Provide the path to the input genome folder, as well as the list of ome codes you made. 

Here is an example:

```bash
module load miniconda

python /project/arsef/projects/bulk_genome_annotation/commands/generate_step1_sort_scripts.py \
--fna_input_dir /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/ncbi_downloads/fna \
--ome_list /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/input_omes.txt
```

Useful options:

`--no_submit` Generate scripts without submitting.

`--include_done` Include genomes even if step1_done is already set. In case you need to re-run something.


You will need to monitor these jobs with:

```bash
squeue --me
```

Once a the sort job is complete for a genome or genomes, you can update the progress tracker sheet:

```bash
module load miniconda # if you haven't loaded it already

python /project/arsef/projects/bulk_genome_annotation/commands/update_step1_progress.py
```

This will verify outputs exist and populate the step1_done column for genomes that finished successfully.

You should also fill in the genus information for each genome (required for Step 4 evidence). After Step 1, open: 

/project/arsef/projects/bulk_genome_annotation/progress/annotation_master_progress.tsv

Fill the genus column for each genome you plan to annotate.

This is critical because Step 4 uses genus to locate evidence files.

<br>

## Step 2 - RepeatModeler + softmask

This step is dependent on step1 being successfully completed and marked as complete on the progress tracker sheet.

This step builds a RepeatModeler database from the sorted genome assembly, runs RepeatModeler to create a custom de novo repeat library, and softmasks the genome assembly with "funannotate mask".

The important files generated from this step are copied safely into: needs_annotation/working_files/<OME>/prep/

```bash
python /project/arsef/projects/bulk_genome_annotation/commands/generate_step2_mask_scripts.py \
  --ome_list /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/input_omes.txt
```

Useful options (although you'll probably never change the defaults):

`--no_submit` Will generate the slurm job scripts, but not submit them. 

`--include_done` or `--force` If either of these options are specified, will submit a job even if the progress sheet indicates that a genome has already successfully completed this step.

`--nodes` Specify how many nodes to ask for per job in the slurm job script. Default is 1 node. 

`--ntasks_per_node` Specify how many threads to ask for per job in the slurm job script. Default is 48 threads. This 

`--time` Specify how much time to ask for per job in the slurm job script. Default is 168 hours (7 days).

`--mem_per_cpu` Specify how much memory per CPU (MB) to ask for per job in the slurm job script. Default is 5000 MB. 

After these jobs are completed, you need to run the proper updating script to mark which genomes have successfully completed this step.

```bash
python /project/arsef/projects/bulk_genome_annotation/commands/update_step2_progress.py
```

Just as before, the script will check for the presence of the necessary output files, then indicate that time/date when the job was verified as completed. 

<br>

## Step 3 - BUSCO training (optional, I prefer to not do this step)

Step 3 is not required for Step 4 and can be safely skipped for genome annotation.

Initially, I designed this pipeline to use the output of the BUSCO training as part of the input for the actual annotation command. I found that it complicated things without much benefit to the final annotation, so I removed this step from the pipeline. 

You can run this command if you really want quick BUSCO results, but I recommend waiting until you've incorporated the final genome assembly/annotation into MycoTools before running BUSCO. That way, you have uniformly formatted input files and you can add the BUSCO results to the lab genome metadata database. 

<br>

## Step 4 - Genome Structural Annotation

For this step to run, there are several dependencies:

* step2 is complete and update_step2_progress.py has been run

* the masked.fasta assembly exists

* the genus column is filled with the genus name of each genome


This step runs "funannotate predict" on the masked genome assembly. 

It also provides additional protein/transcript evidence specific to the genome's specified genus, if we have any available. This protein/transcript evidence was downloaded from JGI, for a maximum of 5 sets of protein/transcript files maximum (if present). The UniProt database will always be provided as protein evidence for each genome, regardless if there are extra protein/transcript data files available for that particular genus.

By default, generate_step4_funannotate_scripts.py generates jobs scripts but does not submit them unless you include --submit in your command.

```bash
python /project/arsef/projects/bulk_genome_annotation/commands/generate_step4_funannotate_scripts.py \
  --ome_list /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/input_omes.txt \
  --submit
```

Once the annotation jobs are complete, you can update the progress sheet with this command:

```bash
python /project/arsef/projects/bulk_genome_annotation/commands/update_step4_progress.py
```

Step 4 uses an intermediate output folder on /90daydata (scratch storage) and then copies the final results to:

* needs_annotation/final_funannotate_results/by_ome/<OME>/

* needs_annotation/final_funannotate_results/busco_short_summaries/

* needs_annotation/final_funannotate_results/gff3/

<br> 

And the bulk genome annotation pipeline is now complete!

You probably want to add these assemblies and annotations to the lab's MycoTools database now. 



