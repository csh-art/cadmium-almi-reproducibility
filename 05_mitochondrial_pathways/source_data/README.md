# Human-readable Figure 5 source data

These CSV files were exported on 2026-10-17 from the frozen M05 Figure 5 objects by the repository audit script. They are provided so the key numerical inputs to Figure 5 can be inspected without loading RData/RDS objects.

Contents:
- `01_group_df.csv`: 40-sample GSE111016 group table (20 control, 20 sarcopenia-related).
- `02_Figure5B_ssGSEA_scores.csv`: four-pathway ssGSEA scores for 40 samples.
- `03_Figure5B_statistics.csv`: group comparison statistics for Figure 5B.
- `04_Figure5C_correlations.csv`: 44 hub-gene × pathway Spearman correlations.
- `05_Figure5D_target_GSEA_statistics.csv`: four target mitochondrial GSEA results.
- `06_Figure5D_running_ES.csv`: running enrichment-score traces.
- `07_Figure5D_gene_set_hits.csv`: gene-set hit positions for the enrichment curves.
- `08_Figure5D_rank_vector.csv`: ranked gene-level metric used for Figure 5D.
- `09_hub_genes.csv`: 11 consensus hub genes and hub classes.
- `10_figure5a_display.csv`: 18 displayed Figure 5A GO-BP terms.
- `10_pathway_size_table.csv`: pathway-size audit for the four mitochondrial pathways.
- `11_Figure5C_pathways_after_hub_removal.csv`: pathway membership after removal of all 11 consensus hub genes.
- `12_Figure5C_ssGSEA_scores_after_hub_removal.csv`: hub-removed pathway scores used for Figure 5C correlations.
- `13_Figure5A_display_terms_from_core.csv`: duplicate audit export of the same 18 Figure 5A display rows as `10_figure5a_display.csv` (byte-identical in this audit).
- `14_Figure5A_full_GO_BP_GSEA_result.csv`: full stored GO-BP GSEA result (5,880 rows).

Audit checks confirm that no consensus hub gene remains in `11_Figure5C_pathways_after_hub_removal.csv`. The original frozen R objects remain archived under `02_analysis_files/`.
