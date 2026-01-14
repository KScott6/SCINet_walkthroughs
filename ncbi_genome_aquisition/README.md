# NCBI genome acquisition

This walkthrough describes how to download genome assemblies, optional annotation files, and structured metadata from NCBI in bulk. This is the first step to incorporating public genome data into our lab's shared MycoTools database. 

This walkthrough will assume you are on the arsef SCINet project and are trying to download to our shared genome database, but you can make this work for your own purposes with minimal changes. 

Ideally, we would use MycoTools to directly download NCBI accessions of interest. However, MycoTools does not automatically store as much metadata as I'd like, and it only incorporates entries that have both an assembly and a corresponding annotation file. So we need to perform this workaround instead.

Please be aware of any limitations/restrictions on any of the NCBI accessions you want to incorporate. 

<br> 

You will provide:

1) A list of taxon names OR a list of NCBI assembly accessions

2) the current version of the MycoTools metadata database
   

And the scripts in this walkthrough will:

1) Pull requested genomes from NCBI (while verifying they are not duplicates of a genome already in the database) 

2) Retrieve structured metadata for each accession and append the metadata from the new genomes into a shared metadata file. 

<br>

---

<br>

## Step 1 - get a list of new genome accessions to download

You will need access to the ncbi_datasets software. If you are on SCINet you can use the shared ncbi_datasets conda environment. If you are doing your own independent genome pull, you need to make your own ncbi_datasets conda environment.

With the fetch_ncbi_metadata_and_merge.py script, you will provide either a list of taxa or a list of accessions (txt or csv). 

If you provide the path to the lab's genome metadata sheet (make sure you're using the most recent version!), the script will automatically deduplicate any incoming genome data. If you do NOT provide a link to a metadata sheet, the script will fetch metadata for any genome matching your desired taxa/accessions, without performing deduplication.

For this example, since the incoming genomes will most likely need to be annotated before being incorporated into MycoTools, I will run this script in the project folder that is associated with our bulk genome annotation pipeline. 

```bash
mkdir /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26
cd /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26

source activate /project/arsef/environments/ncbi_datasets
```

Use the fetch_ncbi_metadata_and_merge.py to take incoming search parameters (taxa names or accessions) and prepare for the download. You can give a custom prefix to your output files with --prefix.

You can provide a csv or txt with a list of taxa names (e.g. Ambrosiella, Fusarium virguliforme, etc.), like so:

```bash
python /project/arsef/scripts/fetch_ncbi_metadata_and_merge.py \
  --taxa_file /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/desired_taxa.txt \
  --master_metadata /project/arsef/databases/mycotools/MTDB_metadata_COMPLETE_07.08.25.csv \
  --outdir /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/ncbi_metadata_by_taxa_py \
  --prefix new_genomes \
  --write_all_fetched
```

Or you can provide a csv or txt of NCBI accessions:

```bash
python /project/arsef/scripts/fetch_ncbi_metadata_and_merge.py \
  --accessions_file /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/desired_accessions.txt \
  --master_metadata /project/arsef/databases/mycotools/MTDB_metadata_COMPLETE_07.08.25.csv \
  --outdir /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/ncbi_metadata_by_acc \
  --prefix new_genomes \
  --write_all_fetched
```

Either option will result in several output files:

1) <prefix>.<accessions/taxa>.ALL_FETCHED.tsv : all of the accessions in NCBI that match your search parameters.
   
2) <prefix>.<accessions/taxa>.MASTER_UPDATED.csv : an unpolished but updated metadata file that includes all the old genomes info, as well as any new unique genome accessions and their metadata. 
   
3) <prefix>.<accessions/taxa>.NEW_ONLY.tsv : a list of the accessions that passed the deduplication step.
   
<br>

## Step 2 - download the genome assemblies (optionally, their annotations as well)

Now you can provide the <prefix>.<accessions/taxa>.NEW_ONLY.tsv file to the download_ncbi_batches.py script. This script will download these new accessions will automatically. You will probably want to run this step as a job if you are downloading a large number of accessions.

Note that the --with_annotation flag controls whether annotation files are requested. All genome assemblies are downloaded regardless.

You will also need to provide a path to the output folder you want to create and download to - this output folder is created automatically.

Here is an example command where I provide the NEW_ONLY.tsv output from my first command, and point the output to a folder I want to make in the bulk genome annotation folder. Make sure to provide your NCBI API key. 

```bash
source activate /project/arsef/environments/ncbi_datasets # activate ncbi_datasets environment if it's not already activated

python /project/arsef/scripts/download_ncbi_batches.py \
  --metadata /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/ncbi_metadata_by_taxa_py/new_genomes.taxa.NEW_ONLY.tsv \
  --outdir /project/arsef/projects/bulk_genome_annotation/needs_annotation/1.14.26/ncbi_downloads \
  --with_annotation \
  --api_key "" # include your API key
```

Congrats! This should result in your genomes and annotations being downloaded in one place. The genome assemblies will be in the "fna" folder and if there were any annotations, they will be in "gff". 

You probably want to annotate your assemblies now - check out the [SCINet Funannotate walkthrough](https://github.com/KScott6/SCINet_walkthroughs/tree/main/funannotate_pipeline/basic_pipeline) to get one step closer to incorporating these genomes into the lab MycoTools database. 

<br>

Notes:

If you are downloading many accessions, you should probably submit this as a job. If you run into any issues with this step, there are many helpful log files that are generated in your specified output folder. 

Although there were no explicitly stated NCBI genome download limits, I started to get intermittent download denials around my 1000th accession and I got completely cut off around my 25000th accession. Between the cutoffs and the occasional normal download failures, it took three attempts to download all the genomes from my list of 5438 accessions.

Of the 5438 genomes I downloaded, only 1447 had associated annotation files (~26%). All data downloaded for these accessions (faa, gff, fna, etc) was only 327GB.

There are a few discrepancies (<30 accessions) between the initial list of accessions I thought I would have and my final list of accessions. I suspect this may be due to different versions of genomes, or version reporting errors. There are so many opportunities where a genome can fall through the cracks in this pipeline - if this ever got to be a more polished workflow this would need to have many extra checks.