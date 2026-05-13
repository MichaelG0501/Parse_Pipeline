#!/bin/bash
set -euo pipefail

WD=/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline
cd $WD

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate velocity
python analysis/trajectory/Auto_prepare_velocity_inputs.py

jid_sort_a=$(qsub analysis/trajectory/Auto_filter_sort_A.pbs)
jid_sort_b=$(qsub analysis/trajectory/Auto_filter_sort_B.pbs)
jid_vel_a=$(qsub -W depend=afterok:${jid_sort_a} analysis/trajectory/Auto_run_velocyto_A.pbs)
jid_vel_b=$(qsub -W depend=afterok:${jid_sort_b} analysis/trajectory/Auto_run_velocyto_B.pbs)
jid_vis=$(qsub -W depend=afterok:${jid_vel_a}:${jid_vel_b} analysis/trajectory/Auto_run_scvelo_visualisation.pbs)

echo "Submitted filter/sort A: ${jid_sort_a}"
echo "Submitted filter/sort B: ${jid_sort_b}"
echo "Submitted velocyto A: ${jid_vel_a}"
echo "Submitted velocyto B: ${jid_vel_b}"
echo "Submitted dependent scVelo visualisation: ${jid_vis}"
