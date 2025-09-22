# SCINet Funnannotate pipeline walkthrough (for one/few genomes)

This walkthrough will guide you through annotating a **single genome assembly**. If you have many genome assemblies, consider using my automated funannotate scripts instead. 

Funannotate is a collection of software used to structurally and functionally annotate (fungal) genomes. I have downloaded and set up a funannotate pipeline so that any ARSEF/SCINet user can easily annotate a genome assembly. 

Please see the bottom of this walkthrough for software version information as well as citation information. 

---

## First-time users

Edit your .condarc file to include the following:

```
pkgs_dirs:  

    - /project/arsef/environments/pkgs
```

Obtain the GeneMark key. This key is only good for 400 days and will need to be replaced **February 12, 2026**.

```bash
cp /project/arsef/environments/funannotate/__external_software/gm_key.gz ~
gunzip ~/gm_key.gz
mv ~/gm_key ~/.gm_key
```

Add the software dependencies to your PATH. 

```bash
echo 'export PATH=/project/arsef/environments/funannotate/__external_software/:$PATH' >>~/.bash_profile
echo 'export PATH=/project/arsef/environments/funannotate/__external_software/eggnog-mapper/:$PATH' >>~/.bash_profile
echo 'export PATH=/project/arsef/environments/funannotate/__external_software/signalp-5.0b/bin/:$PATH' >>~/.bash_profile
echo 'export PATH=/project/arsef/environments/funannotate/__external_software/gmes_linux_64_4/:$PATH' >>~/.bash_profile
echo 'export FUNANNOTATE_DB=/project/arsef/databases/funannotate_databases/' >>~/.bash_profile
echo 'export GENEMARK_PATH=/project/arsef/environments/funannotate/__external_software/gmes_linux_64_4/' >>~/.bash_profile
source ~/.bash_profile
```

Test to see if the conda environment and funannotate pipeline works:

```bash
module load miniconda/24.7.1-2
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# to activate environment
source activate /project/arsef/environments/funannotate

# testing
funannotate check --show-versions

```

If you see the funannotate help menu, it means you are able to use the funannotate pipeline. 

You can try a funannotate test run with provided test data, but I've found that you will probably run into issues due to SCINet security (possible anti-scraping measures).

```bash
funannotate test -t all --cpus 16
```

Once you have everything set up, now all you have to do to access the funannotate environment on a fresh session is run:

```
module load miniconda/24.7.1-2

source activate /project/arsef/environments/funannotate
```

You need to run "**source** activate [env]", not "conda activate [env]". See this SCINet page for more info: [Do not run conda init](https://scinet.usda.gov/guides/software/conda)

To get out of the funannotate environment:

```
conda deactivate
```

---

## Step 0:  File acquisition and folder setup

First, obtain your genome assemblies (.fa, .fasta).

I recommend renaming your genome assembly so that it has a short, unique, and informative name. This name will be used to refer to the assembly and all resulting output files. 

For example, if I was annotating a *Fusarium oxysporum* genome, I would refer to it as **fusoxy1**, and the genome assembly would be **fusoxy1_genomic.fa**.

Then, you need to make a project folder for your genome annotation. For this set of examples, I assume you are working out of the "funannotate_test" directory. Make sure to change the name of "ome" to your chosen genome name.

The following can be run on the login node:

```bash
ome="fusoxy1"
fun_dir="/project/arsef/projects/funannotate_test"
project_dir=${fun_dir}/${ome}

mkdir -p $project_dir
cd $project_dir
mkdir -p data data/evidence/protein data/evidence/EST.transcript commands logs prep repmod funannotate

```

Copy your genome into the /data folder.

```bash
cp fusoxy1_genomic.fasta $project_dir/data
```

---

## Step 0.5: Prepare genome assembly (clean); **OPTIONAL, for haploid genomes only**

The "funannotate clean" option is recommend to only be used for haploid genomes. This step uses Mummer to identify and remove duplicated contigs. Many of our fungal genomes are NOT haploid, or have some duplication inherent to the genome that we want to study, so I would recommend avoiding this step. 

If you do run the "funannotate clean" step, you should submit this command as a job as it can take a while to run.

```bash
ome="fusoxy1"
fun_dir="/project/arsef/projects/funannotate_test"
project_dir="${fun_dir}/${ome}

funannotate clean -i "${project_dir}/data/${ome}_genomic.fasta" -o "${project_dir}/prep/${ome}.clean.fasta" --minlen 1000
```

If Step 0.5 is skipped, the input to "funannotate sort" will still be ${ome}_genomic.fasta. If it is run, the input will be ${ome}.clean.fasta

---

## Step 1:  Prepare genome assembly (sort)

This step sorts the contigs in the genome assembly by size, and renames the assembly headers.

Unless you have a large (>100Mbp) genome, this step is pretty quick and can be run on the login node. 

I have provided an example job script, see [sort.sh](commands/sort.sh) -- make sure to edit the ome name in the SLURM options (#SBATCH). Instead of "fusoxy1", it should be your specified ome name instead.

This is the command you will run:

```bash
### CHANGE OPTIONS HERE ##
ome="fusoxy1"
fun_dir="/project/arsef/projects/funannotate_test"
project_dir="${fun_dir}/${ome}

module load miniconda/24.7.1-2
source activate /project/arsef/environments/funannotate

cd ${project_dir}
funannotate sort -i "${project_dir}/data/${ome}_genomic.fasta" -o "${project_dir}/prep/${OMEcode}.clean.fasta" --minlen 1000
```

### Explanation of Parameters

- `-i "${project_dir}/data/${ome}_genomic.fasta"`  
  Input genome assembly in FASTA format. This should be the raw or cleaned genome you copied into the `data/` folder.

- `-o "${project_dir}/prep/${ome}.sort.fasta"`  
  Output file for the sorted assembly. Contigs will be renamed in a standardized format and ordered by length. Saved in the `prep/` folder.

- `--minlen 1000`  
  Minimum contig length to retain. Any contigs shorter than 1000 bp will be excluded from the output. Edit value to your desired minimum contig length. 
  
---

## Step 2 : Softmask your prepared genome

Next, you need to softmask the prepared genome assembly. We will create a *de novo* repeat library and use that to softmask your genome using RepeatModeler. Since this step can be computationally demanding, you will need to submit these commands as a job. See my example job [mask.sh](commands/mask.sh) -- again, make sure to edit the ome name in the SLURM options (#SBATCH) before you run this yourself.

Since these commands produce a lot of intermediate files, we will output them to the temporary 90daydata folder, then copy over only the most important files to the permanent folder. 

```bash
### CHANGE OPTIONS HERE ##
ome="fusoxy1"
fun_dir="/project/arsef/projects/funannotate_test"
project_dir="${fun_dir}/${ome}"
threads="12"

module load miniconda/24.7.1-2
source activate /project/arsef/environments/funannotate


# Build database from sorted genome assembly
mkdir "/90daydata/arsef/${ome}"
cd "/90daydata/arsef/${ome}"
BuildDatabase -name "${ome}" "${project_dir}/prep/${ome}.sort.fasta"


# Create repeat library
RepeatModeler -threads ${threads} -database "/90daydata/arsef/${ome}/${ome}" -engine ncbi
# copy over output files
cp "/90daydata/arsef/${ome}/${ome}-families.fa" "/90daydata/arsef/${ome}/${ome}-families.stk" "/90daydata/arsef/${ome}/${ome}-rmod.log" "${project_dir}/prep/"


# Softmask genome assembly
# you need to be in the same folder as your sorted assembly and repeat library in order for this all to work
cd "${project_dir}/prep/"
funannotate mask -i "${project_dir}/prep/${ome}.sort.fasta" -m repeatmodeler -l "${project_dir}/prep/${ome}-families.fa" --cpus ${threads} -o "${ome}.masked.fasta" 

```

NOTE: if you are getting strange errors, such as "No such file or directory:  [your masked output file name]", it's probably because you need to install/upgrade the h5py module. Run "pip install h5py" and try again.

### Explanation of Parameters

- `-i "${project_dir}/prep/${ome}.sort.fasta"`  
  Input genome assembly in FASTA format. This is the sorted assembly produced in the previous step.

- `-m repeatmodeler`  
  Specifies the repeat masking method. Here we use the custom repeat library generated by RepeatModeler.

- `-l "${project_dir}/prep/${ome}-families.fa"`  
  Path to the repeat family library generated by RepeatModeler. This file contains the repeat sequences used for masking.

- `--cpus ${threads}`  
  Number of CPU threads to use. Increasing this can speed up masking for larger genomes.

- `-o "${ome}.masked.fasta"`  
  Output FASTA file with repeats softmasked. This masked assembly is used as input for downstream annotation steps.

---

## Step 3 : Funannotate structural annotation 

In this step, Funannotate integrates multiple sources of evidence to predict gene models on the masked genome assembly. This step produces the primary set of gene models that will be carried forward for downstream analyses, so it is typically the most computationally intensive and time-consuming part of the pipeline. You will need to submit this as a job - see my example [predict.sh](/commands/predict.sh) -- again, make sure to edit the ome name in the SLURM options.

There are three things you need to specify for each funannotate run:

1. Transcript and protein evidence

2. The BUSCO linage that best fits your genome 

3. A number to specify the iteration of funannotate you're running (this is in case you want to compare one funannotate run to another, for the same genome)

The process combines *ab initio predictions* (e.g., from Augustus) along with *extrinsic evidence* such as transcript alignments and protein homology from files you provide. I have already made the UniProt database available for you on SCINet, and you should include this as protein evidence. I also encourage you to find other sources of transcript and protein evidence on [JGI](https://mycocosm.jgi.doe.gov/mycocosm/home). 

Try to get transcript evidence (1-5 files) from organisms from the **same genus** as your target genome. Look for JGI files labeled "refined transcripts" or "espressed sequence tags". Avoid files that say "allTranscripts".

You can provide many protein evidence files, and it's less important that they come from a closely-related organism.

Store all your evidence files in the appropriate /evidence subfolder. As you can see in my example, provide the paths to each evidence file separated by a space.


```{bash, eval = FALSE}
### CHANGE OPTIONS HERE ##
ome="fusoxy1"
threads=24
busco_lineage="ascomycota_odb10"
run_number="1"
fun_dir="/project/arsef/projects/funannotate_test"
project_dir="${fun_dir}/${ome}"

module load miniconda/24.7.1-2
source activate /project/arsef/environments/funannotate

cd "${project_dir}"

funannotate predict \
  -i "./prep/${ome}.masked.fasta" \
  -s "${ome}_run_${run_number}" \
  --transcript_evidence ./data/evidence/EST.transcript/FoxII5_GeneCatalog_transcripts_20201108.nt.fasta \
  --protein_evidence ./data/evidence/protein/FoxFo5176_GeneCatalog_proteins_20201111.aa.fasta \
  ./data/evidence/protein/FoxII5_GeneCatalog_proteins_20201108.aa.fasta \
  /project/arsef/databases/funannotate_databases/uniprot_sprot.fasta \
  --cpus ${threads} \
  --optimize_augustus \
  --busco_db "${busco_lineage}" \
  -o ./funannotate/
```

### Explanation of Parameters

- `-i "./prep/${ome}.masked.fasta"`  
  Input genome assembly in FASTA format. This should be the repeat-masked assembly produced in the previous step.

- `-s "${ome}_run_${run_number}"`  
  Strain/project name used to label the output files. Including the run number makes it easy to track multiple annotation attempts for the same genome.

- `--transcript_evidence ./data/evidence/EST.transcript/...fasta`  
  Transcript FASTA file(s) providing RNA-seq–based evidence to support gene model prediction. Multiple files can be included.

- `--protein_evidence ./data/evidence/protein/...fasta ... uniprot_sprot.fasta`  
  Protein FASTA file(s) providing homology-based evidence to support gene model prediction. In this example, species-specific protein sets and the universal UniProt/SwissProt database are included.

- `--cpus ${threads}`  
  Number of CPU threads to use. Annotation is computationally intensive, so allocate plenty of resources.

- `--optimize_augustus`  
  Instructs Funannotate to optimize Augustus parameters for gene prediction using the available training data.

- `--busco_db "${busco_lineage}"`  
  BUSCO lineage dataset to use for training and quality assessment (e.g., `ascomycota_odb10`). This should match the lineage most relevant to your organism.

- `-o ./funannotate/`  
  Output directory for annotation results. Funannotate will generate subdirectories with predicted gene models, annotation files (GFF3, GBK), protein and transcript FASTAs, and summary reports.

---

## Funannotate predict output explained

After `funannotate predict` completes, you should find the following important files in the output directory (`./funannotate/`):

- `predict_results/*.gff3`  
  Genome annotation in GFF3 format (coordinates of predicted genes, transcripts, and features).

- `predict_results/*.gbk`  
  Genome annotation in GenBank format (compatible with NCBI submission pipelines).

- `predict_results/*.faa`  
  Predicted protein sequences in FASTA format.

- `predict_results/*.fna`  
  Predicted transcript (nucleotide) sequences in FASTA format.

- `predict_results/*.tbl`  
  Feature table format used for NCBI genome submissions.

- `predict_misc/`  
  Additional intermediate and support files used during prediction.

- `predict_results/*.txt`  
  Summary report including statistics on gene counts, transcript counts, and evidence support.

These files are the primary outputs you will use for downstream analyses (functional annotation, comparative genomics, submission, etc.).


This pipeline is fine and easy if you are only annotating one or a handful of genomes. However, when you are trying to annotate many genomes all at once it becomes a hassle to generate, organize, and submit these different job scripts, and keep track of where each genome is in the pipeline. 

I created an automated genome annotation pipeline manager script, available on SCINet: see [Funannotate pipeline manager](). **NOTE: I will upload these scripts soon**

This pipeline manager keeps track of all genomes in the pipeline, and works well to generate the annotations necessary for the genome assembly/annotation manager software [MycoTools](https://github.com/xonq/mycotools).

---

## Citations and versions for software/databases used

Funannotate is a wrapper software (incorporates many other software), so you need to cite funannotate itself as well as all the software dependencies. 

 - `Funannotate v1.8.17`
	Palmer, J.M. and Stajich, J., 2020. Funannotate v1. 8.1: Eukaryotic genome annotation. Zenodo, 2020(4054262). https://doi.org/10.5281/zenodo.1134477

- `Augustus v3.5.0`
	Stanke, M., Schöffmann, O., Morgenstern, B. and Waack, S., 2006. Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources. BMC bioinformatics, 7(1), p.62.

- `BUSCO v2`
	Seppey, M., Manni, M. and Zdobnov, E.M., 2019. BUSCO: assessing genome assembly and annotation completeness. In Gene prediction: methods and protocols (pp. 227-245). New York, NY: Springer New York.

Note: Funannotate uses a "modified BUSCO2 script to identify conserved orthologs" (funannotate manual). You have access to BUSCO5 in the funannotate_environment and this is the version that automatically runs when you call "busco", but BUSCO2 is the version used by funannotate.

- `EVidence Modeler`
	Haas, B.J., Salzberg, S.L., Zhu, W., Pertea, M., Allen, J.E., Orvis, J., White, O., Buell, C.R. and Wortman, J.R., 2008. Automated eukaryotic gene structure annotation using EVidenceModeler and the Program to Assemble Spliced Alignments. Genome biology, 9(1), p.R7.
	
- `Exonerate v2.4.0`
	Slater, G.S.C. and Birney, E., 2005. Automated generation of heuristics for biological sequence comparison. BMC bioinformatics, 6(1), p.31.
		
- `GeneMark-ES Suite v4.71_lic`
	Lomsadze, A., Ter-Hovhannisyan, V., Chernoff, Y.O. and Borodovsky, M., 2005. Gene identification in novel eukaryotic genomes by self-training algorithm. Nucleic acids research, 33(20), pp.6494-6506.

- `GlimmerHMM v3.0.4`
	Majoros, W.H., Pertea, M. and Salzberg, S.L., 2004. TigrScan and GlimmerHMM: two open source ab initio eukaryotic gene-finders. Bioinformatics, 20(16), pp.2878-2879.

- `minimap2 v2.28-r1209`
	Li, H., 2018. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), pp.3094-3100.
	
- `Snap v2006-07-28`
	Korf, I., 2004. Gene finding in novel genomes. BMC bioinformatics, 5(1), p.59.

- `tbl2asn v25.8`
	(https://www.ncbi.nlm.nih.gov/genbank/table2asn/)
	
- `tRNAScan-SE v2.0.12`
	Chan, P.P., Lin, B.Y., Mak, A.J. and Lowe, T.M. (2021) "tRNAscan-SE 2.0: 
improved detection and functional classification of transfer RNA genes",
Nucleic Acids Res. 49:9077–9096.
https://doi.org/10.1093/nar/gkab688

- `RepeatModeler v2.0.6`
	Flynn, J.M., Hubley, R., Goubert, C., Rosen, J., Clark, A.G., Feschotte, C. and Smit, A.F., 2020. RepeatModeler2 for automated genomic discovery of transposable element families. Proceedings of the National Academy of Sciences, 117(17), pp.9451-9457.




