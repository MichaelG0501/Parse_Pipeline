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

WD=/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/novogene/Auto_parse_qc_pipeline
ZIP=/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/novogene/8da03a32-a9ad-4938-9338-cbd1731f445d_filtered_matrices.zip

cd "${WD}"
mkdir -p parse_qc_outs parse_qc_outs/input parse_qc_outs/logs parse_qc_outs/plots parse_qc_outs/summary parse_qc_outs/by_samples

if [ ! -d parse_qc_outs/input/output_combined ]; then
  unzip -oq "${ZIP}" -d parse_qc_outs/input
fi

Rscript QC_Pipeline.R
Rscript summary.R
Rscript ../Visualization_scripts/qc_heatmap.R

echo "$(date +%T)"
