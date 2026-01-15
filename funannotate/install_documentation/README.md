# Documentation of SCINet Funannotate pipeline install

Funannotate is a collection of software used to structurally and functionally annotate (fungal) genomes. Unfortunately, the copy of Funannotate on SCINet will not allow you to perform annotations with custom parameters. So, I set up a funannotate pipeline able to be used on SCINet, with customizable options. 

This page contains the documentation of how I set everything up- software dependencies, etc. I set up a shared conda environment for funannotate; packages and dependencies are saved to the `/project/arsef` folder rather than my home folder.

---

## Creating shared conda environment

All funannotate packages are available in a shared environment folder. If you want to export the funannotate environment to a new HPC or make a copy, try cloning it from the environment yaml:

[Creating an environment from a yaml](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

[Sharing a conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#sharing-an-environment)


In order to download/use shared libraries, you need to first edit your .condarc to include the following:


pkgs_dirs:  

    - /project/arsef/environments/pkgs
    
    
To create the environment and install the base funannotate software:

```bash
module load miniconda/24.7.1-2
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create "python>=3.6,<3.9" --prefix /project/arsef/environments/funannotate_env

source activate /project/arsef/environments/funannotate_env

conda install funannotate=1.8.17 #need to specify this EXACT version or it will not work
```

---

## Installing Funannotate databases

```bash
echo "export FUNANNOTATE_DB=/project/arsef/databases/funannotate_databases/" > /project/arsef/environments/funannotate_working/etc/conda/activate.d/funannotate.sh   

echo "unset FUNANNOTATE_DB" > /project/arsef/environments/funannotate_working/etc/conda/deactivate.d/funannotate.sh 

pip install distro

funannotate setup -i all -d /project/arsef/databases/funannotate_databases/

```

I had some trouble installing the BUSCO databases via Funannotate, so I manually downloaded many of our commonly used fungal / insect BUSCO databases and made them available. 

```bash
busco --download_path /project/arsef/databases/busco_datasets/ --download eukaryota

cp -R agaricales_odb10/ agaricomycetes_odb10/ arachnida_odb10/ arthropoda_odb10/ ascomycota_odb10/ basidiomycota_odb10/ boletales_odb10/ diptera_odb10/ dothideomycetes_odb10/ eukaryota_odb10/ eurotiomycetes_odb10/ eurotiales_odb10/ fungi_odb10/ glomerellales_odb10/ helotiales_odb10/ hemiptera_odb10/ hymenoptera_odb10/ hypocreales_odb10/ insecta_odb10/ leotiomycetes_odb10/ lepidoptera_odb10/ nematoda_odb10/ pleosporales_odb10/ polyporales_odb10/ saccharomycetes_odb10/ sordariomycetes_odb10/ stramenopiles_odb10/ tremellomycetes_odb10/ /project/arsef/databases/funannotate_databases/

```

I have most fungal and insect/arthropod BUSCO odb10 datasets stored in /project/arsef/databases/funannotate_databases/ and they should be freely available. If you find yourself wanting a dataset that I haven't already downloaded, you can do one of two things:

1) Check the separate BUSCO dataset folder (/project/arsef/databases/busco_datasets/lineages/)

I downloaded all obd_10 eukaryora BUSCO datasets 1/30/25, as time goes on I'll probably add more datasets. If your desired dataset is stored in the busco_datasets/lineage/ folder, just copy the folder into the funannotate database folder **and change permissions on that folder**. Like so:

```bash
cp -R /project/arsef/databases/busco_datasets/lineages/aves_odb10 /project/arsef/databases/funannotate_databases/

chmod -R 770 /project/arsef/databases/funannotate_databases/aves_odb10
```

Or you can download a new dataset via BUSCO. It's just one extra step. However, it is **still important to change permissions** once you have downloaded the file. 

```bash
# downloading aves_odb10
busco --download_path /project/arsef/databases/busco_datasets/ --download aves_odb10

cp -R /project/arsef/databases/busco_datasets/lineages/aves_odb10 /project/arsef/databases/funannotate_databases/

chmod -R 770 /project/arsef/databases/funannotate_databases/aves_odb10
```

---

## Installing software dependencies

I am storing all manually downloaded software in: /project/arsef/environments/funannotate_env/__external_software/

### GeneMark

Note:  the GeneMark key will expire once every 400 days. You will need to get a new key from the [GeneMark website)[https://exon.gatech.edu/] once your old one expires (remove old .gm_key file, go to the website, fill out their form, download the key, upload the key to your home folder, and rename the new file to .gm_key). 

**The .gm_key I provided will expire February 12, 2026.**


I downloaded genemark GeneMark-ES/ET/EP+ ver 4.72_lic (LINUX 64 kernel 3.10 - 5) via their [website](https://genemark.bme.gatech.edu/license_download.cgi)

```bash
cd /project/arsef/environments/funannotate/__external_software/
tar xvzf /project/arsef/environments/funannotate/__external_software/gmes_linux_64_4.tar.gz
cp /project/arsef/environments/funannotate/__external_software/gm_key.gz ~
gunzip ~/gm_key.gz
mv ~/gm_key ~/.gm_key # key will expire 02/12/2026
```

GeneMark often uses the wrong perl installation. This is such a common problem that there is actually an official GeneMark script to change all the GeneMark scripts to use the proper perl. This is how you fix it:

```bash
cd /project/arsef/environments/funannotate/__external_software/gmes_linux_64_4/
which perl # make note of path
# /project/arsef/environments/funannotate/bin/perl
perl ./change_path_in_perl_scripts.pl /project/arsef/environments/funannotate/bin/perl
```

### eggnog-mapper (emapper.py)

```bash
git clone https://github.com/eggnogdb/eggnog-mapper.git

# also downloaded eggnog_proteins.dmnd database, taxonomy database, and diamond database via the following:
eggnog-mapper/download_eggnog_data.py
# saved to: /eggnog-mapper/data/

# for emapper.py
pip install psutil
pip install Bio
```

### SignalP

I downloaded signalp (v5, Linux) via their [website](https://services.healthtech.dtu.dk/services/SignalP-6.0/)

```bash
cd /project/arsef/environments/funannotate/__external_software/
tar xvzf /project/arsef/environments/funannotate/__external_software/signalp-5.0b.Linux.tar.gz
```

## EvidenceModeler

Hopefully this fixes the EvidenceModeler error that pops up about 1.5hr into a funannotate predict run

```bash
cd /project/arsef/environments/funannotate/opt/evidencemodeler-2.1.0
find . -name "evidence_modeler.pl" # take note of path
chmod +x ./EvmUtils/evidence_modeler.pl
ln -s ./EvmUtils/evidence_modeler.pl ./evidence_modeler.pl
```

## InterProScan

```bash
cd /project/arsef/environments/funannotate/__external_software/
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.73-104.0/interproscan-5.73-104.0-64-bit.tar.gz
tar xvzf /project/arsef/environments/funannotate/__external_software/interproscan-5.73-104.0-64-bit.tar.gz
cd /project/arsef/environments/funannotate/__external_software/interproscan-5.73-104.0/
python3 setup.py -f interproscan.properties
```

## Assorted software via conda

```bash
# for repeatmasker/repeatmodeler 
conda install bioconda::repeatmodeler

# for BUSCO
conda install bioconda::busco

# older version of h5py distributed in one of the packages, make sure to upgrade or repeatmodeler won't work
pip install --upgrade h5py

```

## Putting software in PATH

```bash
echo 'export PATH=/project/arsef/environments/funannotate/__external_software/:$PATH' >>~/.bash_profile
echo 'export PATH=/project/arsef/environments/funannotate/__external_software/eggnog-mapper/:$PATH' >>~/.bash_profile
echo 'export PATH=/project/arsef/environments/funannotate/__external_software/signalp-5.0b/bin/:$PATH' >>~/.bash_profile
echo 'export PATH=/project/arsef/environments/funannotate/__external_software/gmes_linux_64_4/:$PATH' >>~/.bash_profile
echo 'export PATH=/project/arsef/environments/funannotate/__external_software/interproscan-5.73-104.0/:$PATH' >>~/.bash_profile
echo 'export FUNANNOTATE_DB=/project/arsef/databases/funannotate_databases/' >>~/.bash_profile
echo 'export GENEMARK_PATH=/project/arsef/environments/funannotate/__external_software/gmes_linux_64_4/' >>~/.bash_profile
source ~/.bash_profile
```

## antiSMASH environment

I tried to install antiSMASH into the funannotate environment, but had conflicting libraries. Hopefully these conflicts will be fixed with Funannotate2, but in the meantime I have also created a shared conda environment for antiSMASH, if this is needed for your analyses.

```bash
cd /project/arsef/environments/
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create --prefix /project/arsef/environments/antismash biopython==1.76
conda activate /project/arsef/environments/antismash
conda install prodigal hmmer hmmer2 blast diamond fasttree meme==4.11.2 glimmerhmm
download-antismash-databases --database-dir /project/arsef/databases/antismash_databases/
antismash --check-prereqs --database /project/arsef/databases/antismash_databases/
```









