#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --partition=ceres
#SBATCH --job-name=fun_mask
#SBATCH --account=arsef
#SBATCH --output=/project/arsef/projects/funannotate_test/fusoxy1/logs/mask_%j.o
#SBATCH --error=/project/arsef/projects/funannotate_test/fusoxy1/logs/mask_%j.e

### CHANGE OPTIONS HERE ##
ome="fusoxy1"
fun_dir="/project/arsef/projects/funannotate_test"
project_dir="${fun_dir}/${ome}"
threads=12

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