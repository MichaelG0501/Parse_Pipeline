# High-Resolution MP TCGA Survival Volcano Methodology

Script: `analysis/metaprograms/parse_highres_mp_tcga_survival_volcano.R`

Status: active terminal metaprogram clinical-association figure workflow.

## Purpose

This workflow tests whether selected high-resolution Parse metaprograms are associated with overall survival in the TCGA ESCA cohort. It reproduces the reference TCGA mode from `scRef_Pipeline/analysis/clinical/tcga_mp_state_survival_reg_noreg.R`: `(b) GSVA on whole TCGA profile (legacy/reference)`.

## Inputs

- Strict high-resolution MP genes: `parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_selected_mp_genes_nMP117.rds`
- Strict trend summary: `parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_trend_summary_nMP117.csv`
- Strict high-resolution MP preliminary labels: `parse_outs/Auto_parse_highres_metaprogram_trends/Auto_parse_highres_top_3CA_noncellcycle_nMP117.csv`
- Legacy T2/T4-high MP genes: `parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_selected_mp_genes_nMP117.rds`
- Legacy T2/T4-high MP preliminary labels: `parse_outs/Auto_parse_highres_metaprogram_trends/Auto_T2T4_gt_T0eR4_filter/Auto_parse_highres_T2T4_top_3CA_noncellcycle_nMP117.csv`
- TCGA metadata compatibility copy: `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/tcga_esca_meta.rds`
- TCGA whole-profile TPM compatibility copy: `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/cibersortx/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt`
- Reconstructed TCGA fallback paths under `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/TCGA/esca_gdc_reconstruction/`

## Methods

The script builds three MP groups:

- `strict_increase`: retained strict mean/median MPs with `treatment_direction == "increase"`.
- `strict_decrease`: retained strict mean/median MPs with `treatment_direction == "decrease"`.
- `legacy_t2t4_high`: MPs retained by the legacy T2/T4-high comparison workflow.

Each MP gene list is intersected with the TCGA whole-profile TPM matrix and retained for GSVA if at least five genes overlap. GSVA is run on the whole TCGA expression profile with Gaussian kernel settings, matching the reference whole-TCGA mode. The primary requested plots use the continuous Cox proportional hazards model in EAC primary tumour samples:

`Surv(OS_time, OS_event) ~ GSVA_score`

For parity with the reference clinical script, the output table and multi-page PDF also include the median split and upper-versus-lower-quartile split models. The volcano x-axis is `log2(HR)` and the y-axis is `-log10(p)`, matching the reference clinical script. MP point labels use the preliminary MP name/description from the first non-cell-cycle 3CA hit. Only nominal Cox p-values are reported; no multiple-testing correction column is written.

## Outputs

Output root:

`parse_outs/highres_mp_tcga_survival/`

Output tiers:

- `intermediate/`: cached whole-TCGA GSVA score matrix and filtered gene sets.
- `tables/`: Cox model results and gene-set overlap summaries.
- `figures/`: one multi-page volcano PDF plus individual PDF/PNG volcano plots for increasing, decreasing, and T2T4-high MPs. Unsuffixed group plots are the primary continuous-Cox plots; `_median` and `_q1q4` suffixed plots are reference split-model companions.
- `reports/`: short text summary of inputs, sample counts, and significance counts.
- `logs/`: standard run summary under `parse_outs/logs/run_summaries/`.

## Interpretation Notes

These plots are TCGA bulk-expression survival associations, not Parse single-cell differential tests. The whole-profile GSVA mode is intentionally retained as a legacy/reference comparator because it includes the complete TCGA tumour profile rather than malignant-only deconvolution.
