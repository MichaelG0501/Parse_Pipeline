#!/bin/bash
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=12:00:00
#PBS -N Auto_parse_qc
#PBS -koed

echo "$(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline
ZIP=/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline/8da03a32-a9ad-4938-9338-cbd1731f445d_filtered_matrices.zip

cd "${WD}"
mkdir -p parse_outs parse_outs/input parse_outs/logs parse_outs/plots parse_outs/summary parse_outs/by_samples

if [ ! -d parse_outs/input/output_combined ]; then
  unzip -oq "${ZIP}" -d parse_outs/input
fi

Rscript Auto_parse_QC_Pipeline.R
Rscript analysis/summary/Auto_parse_summary.R
Rscript analysis/plotting/Auto_parse_qc_heatmap.R

echo "$(date +%T)"
