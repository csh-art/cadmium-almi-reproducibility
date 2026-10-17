# Cross-cohort validation (Figure 4)

This module archives the cohort-specific result tables and plotting provenance used for Figure 4 cross-cohort evaluation.

Cohorts: GSE111016 (discovery), GSE111010 (main validation), GSE111006 (same-SuperSeries sensitivity), and GSE167186 (independent phenotype-heterogeneity sensitivity-validation). Raw expression matrices are not pooled.

## Final panel-label note
The final manuscript labels the seven-pathway concordance panel as **Figure 4D** and the four-cohort hub-gene forest as **Figure 4E**. Historical analysis/output filenames retain the earlier internal naming (`Figure4E_pathway...` for the final panel D source and `Figure4D_multicohort_forest...` for the final panel E source). These historical names are preserved for provenance.

## Final panel D — pathway concordance
`source_data/Figure4E_pathway_replication_source_data.csv` is the locked 7-pathway x 4-cohort (28-row) table. The original standalone plotting code is archived under `archived_original_code/`. A repository-relative portable plotting script is provided under `scripts/01_plot_Figure4D_pathway_concordance_portable.R`.

## Final panel E — multi-cohort hub-gene forest
The archived 44-row source is `source_data/Figure4D_multicohort_forest_source_data.csv`. The original historical report documents an incomplete capture of the old plotting script. However, all 41 available estimates retain both the unshrunk log2 fold-change and Wald statistic. The missing sensitivity-cohort standard errors are reconstructed from the DESeq2 Wald identity and checked against the 22 rows with explicitly stored standard errors. See `audit/Figure4E_reconstruction_audit.md`.

`source_data/Figure_04E_multicohort_forest_reconstructed_plot_data.csv` contains the reconstruction audit table, and `scripts/02_plot_Figure4E_multicohort_forest_reconstructed.R` is a transparent portable reconstruction of the final panel. It is not represented as the lost original code.

## Composite
The final historical aspect-ratio-preserving A-E assembly code is archived under `archived_original_code/Figure_04_final_plain_PRESERVE_ASPECT_code_original.txt`. It confirms the final layout: A/B in row 1, C/D in row 2, and E full width.

## Archived source data
- `07_GSE111010_predefined_11_hub_gene_validation.csv`
- `12_GSE111006_Figure4_replication_input.csv`
- `06_GSE167186_Figure4_replication_input.csv`
- `04_Figure4D_multicohort_hub_replication_input.csv`
- `Figure4E_pathway_replication_source_data.csv`
- `Figure4D_multicohort_forest_source_data.csv`
- `Figure_04E_multicohort_forest_reconstructed_plot_data.csv`

## Status
The Figure 4 analytical/source-data evidence chain is now closed for repository purposes. The only caveat is that the original historical plotting script for the four-cohort forest was not recoverable; a clearly labeled reconstruction is supplied instead. Final public-release cleanup should convert any remaining archived absolute-path scripts into optional provenance-only material and locally test the portable scripts in the recorded R environment.
