# SUR1090 Count Matrix Preparation Methodology

Script: `analysis/input_preparation/parse_prepare_sur1090_count_matrices.R`

Status: active input-preparation utility.

## Purpose

This script converts SUR1090 treated and untreated PDO count matrices from CSV format into the `DGE_filtered` layout expected by the Parse QC pipeline.

## Inputs

External input directory:

`/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/00_counts_matrix_all`

Files:

- `SUR1090_Treated_PDO.csv`
- `SUR1090_Untreated_PDO.csv`

The CSV files are expected to have genes in the first column and cell barcodes in the remaining column names.

## Processing

For each SUR1090 sample:

1. Read the CSV with `data.table::fread()`.
2. Store the first column as gene names.
3. Convert the remaining numeric matrix from genes x cells to cells x genes.
4. Convert the transposed matrix to `dgCMatrix`.
5. Write the sparse count matrix as Matrix Market format.
6. Write gene metadata with `gene_name` and empty `gene_id`.
7. Compute per-cell `nCount_RNA` and `nFeature_RNA` from the original genes x cells matrix.
8. Write cell metadata with `bc_wells`, `nCount_RNA`, `nFeature_RNA`, and placeholder `percent.mt`.

## Outputs

Output root:

`parse_outs/input/output_combined/<sample>/DGE_filtered/`

Output tiers:

- `intermediate`: `count_matrix.mtx.gz`, `all_genes.csv.gz`, `cell_metadata.csv.gz`
- `logs`: `parse_outs/logs/run_summaries/parse_prepare_sur1090_count_matrices_*.txt`

## Downstream Use

The QC pipeline consumes these `DGE_filtered` folders before downstream `analysis/` scripts run.
