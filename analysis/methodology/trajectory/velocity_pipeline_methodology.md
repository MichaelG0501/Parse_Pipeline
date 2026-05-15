# Parse RNA Velocity Pipeline Methodology

Scripts:

- `analysis/trajectory/parse_velocity_prepare_inputs.py`
- `analysis/trajectory/parse_velocity_filter_sort_NACT1090_A.pbs`
- `analysis/trajectory/parse_velocity_filter_sort_NACT1090_B.pbs`
- `analysis/trajectory/parse_velocity_run_velocyto_parse_bam.py`
- `analysis/trajectory/parse_velocity_run_velocyto_NACT1090_A.pbs`
- `analysis/trajectory/parse_velocity_run_velocyto_NACT1090_B.pbs`
- `analysis/trajectory/parse_velocity_scvelo_visualise.py`
- `analysis/trajectory/parse_velocity_run_scvelo_visualisation.pbs`
- `analysis/trajectory/parse_velocity_submit_pipeline.sh`

Status: active RNA velocity workflow.

## Purpose

This workflow prepares Parse/Trailmaker BAM inputs, runs velocyto on NACT1090 A/B non-PDO cells, and generates scVelo visualizations for the six Parse treatment-response samples.

## Inputs

Trailmaker inputs:

- `data/trailmaker/output_NACT1090_A_MKDL260004725-1A_23JFJWLT3_L6_REP_CLEAN/`
- `data/trailmaker/output_NACT1090_B_MKDL260004726-1A_23JFJWLT3_L6_REP_CLEAN/`

BAM inputs:

- `.../process/barcode_headAligned_anno.bam`

Reference input:

- `/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz`

Download if absent:

- `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz`

## Processing

1. `parse_velocity_prepare_inputs.py`
   - Extract non-PDO barcode lists and metadata.
   - Write combined metadata tables.
   - Re-map GTF chromosome names to Parse/Trailmaker BAM convention.
   - Download and convert RepeatMasker hg38 annotations if absent.
   - Create BAM symlinks under the velocity output folder.
2. `parse_velocity_filter_sort_NACT1090_[A|B].pbs`
   - Filter each BAM to non-PDO barcodes.
   - Sort and index coordinate BAMs.
3. `parse_velocity_run_velocyto_NACT1090_[A|B].pbs`
   - Run the Parse-aware velocyto wrapper on each filtered BAM.
4. `parse_velocity_scvelo_visualise.py`
   - Load A/B loom files.
   - Merge with non-PDO metadata.
   - Run scVelo preprocessing, velocity estimation, and embedding visualization.
   - Save velocity figures, cell metadata, and H5AD object.
5. `parse_velocity_submit_pipeline.sh`
   - Submits the above PBS jobs with dependencies.

## Outputs

Output root:

`parse_outs/Auto_velocity/`

Output tiers:

- `intermediate`: `coord/*.bam`, `looms/*/*.loom`, `Auto_scvelo_nonPDO.h5ad`
- `tables`: barcode lists and velocity metadata CSV files
- `figures`: `figures/Auto_velocity_*.png`, UMAP quality PNG files
- `logs`: PBS stdout/stderr files under `logs/`
- `reports`: no combined report currently

## Downstream Use

Terminal velocity workflow. No active downstream R script consumes the velocity outputs.
