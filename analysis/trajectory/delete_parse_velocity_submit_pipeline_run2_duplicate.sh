#!/bin/bash
set -euo pipefail

# DELETE CANDIDATE: duplicate of parse_velocity_submit_pipeline.sh retained only
# so the user can remove it manually.

WD=/rds/general/ephemeral/project/spatialtranscriptomics/ephemeral/Parse_Pipeline
cd $WD

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate velocity
python analysis/trajectory/parse_velocity_prepare_inputs.py

jid_sort_a=$(qsub analysis/trajectory/parse_velocity_filter_sort_NACT1090_A.pbs)
jid_sort_b=$(qsub analysis/trajectory/parse_velocity_filter_sort_NACT1090_B.pbs)
jid_vel_a=$(qsub -W depend=afterok:${jid_sort_a} analysis/trajectory/parse_velocity_run_velocyto_NACT1090_A.pbs)
jid_vel_b=$(qsub -W depend=afterok:${jid_sort_b} analysis/trajectory/parse_velocity_run_velocyto_NACT1090_B.pbs)
jid_vis=$(qsub -W depend=afterok:${jid_vel_a}:${jid_vel_b} analysis/trajectory/parse_velocity_run_scvelo_visualisation.pbs)

echo "Submitted filter/sort A: ${jid_sort_a}"
echo "Submitted filter/sort B: ${jid_sort_b}"
echo "Submitted velocyto A: ${jid_vel_a}"
echo "Submitted velocyto B: ${jid_vel_b}"
echo "Submitted dependent scVelo visualisation: ${jid_vis}"
