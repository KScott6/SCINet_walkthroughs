#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=ceres
#SBATCH --job-name=fun_predict
#SBATCH --account=arsef
#SBATCH --output=/project/arsef/projects/funannotate_test/fusoxy1/logs/predict_%j.o
#SBATCH --error=/project/arsef/projects/funannotate_test/fusoxy1/logs/predict_%j.e

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
