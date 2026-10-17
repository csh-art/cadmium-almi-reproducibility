# Figure 4E reconstruction audit

The historical `Figure4D_multicohort_forest_report.txt` was retained for provenance but is **not** treated as the final execution authority. Its own audit reports that the original plotting-code capture was incomplete and that only 22 available estimates had an explicit stored `lfcSE`.

The archived 44-row `Figure4D_multicohort_forest_source_data.csv` nevertheless contains `log2FC_unshrunk` and `Wald_stat` for all 41 available gene-cohort estimates. For rows with an explicit `lfcSE`, `abs(log2FC_unshrunk / Wald_stat)` reproduces the recorded `lfcSE` to floating-point precision (maximum absolute difference = 1.221e-15, n = 22). Therefore, the sensitivity-cohort standard errors and 95% Wald confidence intervals can be reconstructed without changing any effect estimate or p value.

Checks performed in v6:
- source rows: 44
- unique genes: 11
- cohorts: 4
- available estimates: 41
- unavailable GSE167186 estimates: 3
- final SE available after reconstruction: 41/41 available estimates
- duplicate uploaded source-data copies had identical MD5: `ffcc4361266ccc4eb6dddf419cc64e5c`; one copy was retained

The portable script `scripts/02_plot_Figure4E_multicohort_forest_reconstructed.R` is explicitly labeled as a **reconstruction**, not as the lost original plotting script.
