============================================================
Figure 3A original candidate-gene screening record
Figure 3A 原始候选基因筛选记录
============================================================

Created: 2026-07-16 19:01:44

AUTHORITATIVE SCREENING LOGIC / 权威筛选逻辑

1. CTD cadmium-related genes:
   interaction_count > 10
   Result: 138 genes.

2. GSE111016 differential-expression analysis:
   DESeq2 design: ~ sarcopenia
   Comparison: sarcopenia yes versus no
   Reference group: no
   All 65,217 count-matrix rows entered DESeq2.
   No low-expression prefiltering was performed.

3. Nominal DEG pool displayed in Figure 3A:
   pvalue < 0.05
   Result: 3,342 Ensembl-level DESeq2 results.
   This is the value displayed in the right circle of Figure 3A.

4. Original candidate-prioritization rule:
   pvalue < 0.05 AND abs(log2FoldChange) > 0.1

5. Intersection:
   CTD 138 genes were intersected with prioritized
   GSE111016 genes after mapping to official gene symbols.
   Result: 23 candidate genes.

IMPORTANT INTERPRETATION / 重要解释

The 23 genes are not the direct intersection obtained using
pvalue < 0.05 alone. The direct nominal-P intersection contains
24 genes, including SLC11A2.

SLC11A2 satisfies nominal pvalue < 0.05 and is present in the
CTD 138-gene set, but its absolute log2 fold change is below
0.1. It is therefore excluded automatically by the original
minimal effect-size criterion; it was not manually deleted.

Figure 3A should therefore be interpreted as an integration and
candidate-prioritization schematic rather than a strict
unfiltered two-set Venn diagram.

中文说明：

Figure 3A右侧的3342表示DESeq2分析中nominal P < 0.05的
全部结果。最终23个候选基因在此基础上进一步满足
|log2FC| > 0.1，并与CTD的138个基因在Gene Symbol层面
取交集。SLC11A2由于|log2FC|未超过0.1而自动排除，
并非后期人工删除。

FINAL 23 CANDIDATE GENES / 最终23个候选基因

ABCG2, ATF4, CASP3, CASP9, CCL2, CYCS, CYP1A1, DDIT3, FOXO1, FOXO3, HIF1A, HMGA2, MAP1LC3B, MFN2, MMP9, MTF1, PHGDH, PKM, PPARGC1A, SLC16A1, SP1, STAT3, VEGFA

FILES

00_original_screening_parameters.csv
    Complete thresholds and methodological settings.

01_input_file_manifest.csv
    Source paths, file sizes, modification times and MD5 values.

02_CTD_cadmium_interaction_count_gt10_138_genes.csv/txt
    Reproduced CTD 138-gene set.

03_GSE111016_aligned_phenotype_20_no_20_yes.csv
    Verified sample alignment and group assignments.

04_GSE111016_DESeq2_all_results_65217_Ensembl.csv
    Complete unfiltered DESeq2 result table.

05_GSE111016_nominal_P005_DEGs_3342_Ensembl.csv
    Nominal DEG pool displayed as n = 3342 in Figure 3A.

06_GSE111016_original_priority_P005_absLog2FC_gt0.1_Ensembl.csv
    Original prioritized Ensembl-level DEG set.

07_GSE111016_original_priority_P005_absLog2FC_gt0.1_unique_Symbol.csv
    Prioritized DEG set after symbol mapping and deduplication.

08_CTD138_nominal_P005_direct_overlap_24_genes_AUDIT.csv
    Audit-only direct overlap before the |log2FC| criterion.

09_CTD138_original_priority_overlap_all_Ensembl.csv
    CTD-prioritized overlap before final symbol deduplication.

10_FINAL_original_23_candidate_genes.csv/txt
    Authoritative final candidate-gene set.

11_SLC11A2_threshold_boundary_AUDIT.csv
    Explicit record explaining the SLC11A2 boundary case.

12_original_screening_summary.csv
    Key counts and reproducibility checks.

13_original_Figure3A_screening_record.rds
    Consolidated R object containing all important tables.

14_GSE111016_DESeq2_dds_unfiltered_65217.rds
14_GSE111016_DESeq2_results_object.rds
    Saved DESeq2 objects for future reuse.

15_sessionInfo_and_package_versions.txt
    R environment and package versions.

Do not replace this original screening rule with a human-only
recount of the current CTD interaction table, and do not apply
low-expression filtering before DESeq2 when reproducing Figure 3A.

============================================================
