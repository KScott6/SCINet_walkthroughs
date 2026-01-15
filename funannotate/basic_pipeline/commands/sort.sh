#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --partition=ceres
#SBATCH --job-name=fun_sort
#SBATCH --account=arsef
#SBATCH --output=/project/arsef/projects/funannotate_test/fusoxy1/logs/sort_%j.o
#SBATCH --error=/project/arsef/projects/funannotate_test/fusoxy1/logs/sort_%j.o


### CHANGE OPTIONS HERE ##
ome="fusoxy1"
fun_dir="/project/arsef/projects/funannotate_test"
project_dir="${fun_dir}/${ome}

module load miniconda/24.7.1-2
source activate /project/arsef/environments/funannotate

cd ${project_dir}
funannotate sort -i "${project_dir}/data/${ome}_genomic.fasta" -o "${project_dir}/prep/${OMEcode}.clean.fasta" --minlen 1000
