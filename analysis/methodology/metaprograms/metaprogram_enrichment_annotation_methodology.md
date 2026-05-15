# Parse Metaprogram Enrichment Annotation Methodology

Script: `analysis/metaprograms/parse_metaprogram_enrichment_annotation.R`

Status: active terminal annotation workflow.

## Purpose

This script annotates selected Parse metaprograms against Hallmark, GO biological process, 3CA pan-cancer metaprograms, and developmental-stage reference signatures.

## Inputs

- `parse_outs/Auto_parse_metaprograms/Auto_parse_MP_outs_default.rds`
- 3CA metaprograms: `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`
- Developmental enrichment RDS files: `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/*.rds`
- MSigDB Hallmark via `msigdbr`
- GO:BP via `clusterProfiler` and `org.Hs.eg.db`

## Processing

1. Load reference gene sets.
2. Load selected Parse metaprogram object.
3. Filter negative-silhouette MPs.
4. Filter low sample-coverage MPs.
5. Run enrichment for each retained MP:
   - GO biological process.
   - Hallmark.
   - 3CA MPs.
   - Each developmental custom reference.
6. Save the full enrichment object.
7. Generate multi-page PDF and individual PNG heatmaps.

## Outputs

Output tiers:

- `intermediate`: `Auto_parse_cluster_enrich.rds`
- `reports`: `Auto_parse_enrichment_annotation.pdf`
- `figures`: `Auto_parse_enrich_<reference>.png`
- `logs`: `parse_outs/logs/run_summaries/parse_metaprogram_enrichment_annotation_*.txt`

## Downstream Use

Terminal figure/report workflow. No active downstream script consumes the enrichment plots.
