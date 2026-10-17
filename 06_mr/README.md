# Gene-level MR (M11_MR / Figure 6)

## Current status

The **frozen source/provenance evidence chain is closed** for M11/Figure 6. Final result tables, outcome metadata, GTEx target-gene exposure summaries, exact/position-recovered outcome matching, harmonization, EUR-reference overlap, primary and sensitivity clumping, sex-stratified extraction, heterogeneity, Steiger directionality, and targeted QC tables are archived.

This close-out does **not** mean that large public GTEx V11, ALM GWAS, or 1000 Genomes files are redistributed. It also does not relabel reconstructed validation code as recovered original upstream code.

## Final source tables

- `source_data/Table_MR_instrument_characteristics_and_results.tsv`
- `source_data/Table_MR_Wald_ratio_combined_ALM_final.tsv`
- `source_data/Table_MR_Wald_ratio_sex_stratified_ALM.tsv`
- `source_data/Table_MR_Steiger_directionality.tsv`
- `source_data/MR_10_final_independent_instruments_combined_ALM.tsv`
- `source_data/Table_MR_female_male_effect_difference_final.tsv`
- `source_data/Table_MR_LD_clumping_sensitivity.tsv`

## Frozen upstream chain

- `MR_05`: 342 strict cis-eQTL candidate rows.
- `MR_06`: 327 exact combined-ALM matches; 15 variants initially missing.
- `MR_06c`: 6 position-recovered variants.
- `MR_07`: 333 harmonized candidates (327 exact + 6 recovered).
- EUR-reference overlap: 309 present, 24 absent.
- `MR_10`: four independent retained instruments.

The four retained instruments are ATF4/rs882565, CASP9/rs933704, CASP3/rs4647601, and MAP1LC3B/rs7865.

Primary PLINK clumping used `r2 < 0.001` within 10 Mb. Relaxing the threshold to `r2 < 0.01` retained the same four lead variants.

## Result anchors

- F statistics: 311.65, 106.35, 72.20, and 43.46.
- ATF4 combined Wald ratio: beta=-0.016503, 95% CI -0.028419 to -0.004586, nominal p=0.006642, Bonferroni p=0.026568.
- Female-male heterogeneity for ATF4: p=0.291184.
- Steiger supports expression -> ALM for all four retained instruments.
- rs7865 is a palindromic C/G variant, but the frozen audit records exact same alignment and allele frequencies far from 0.5.

## Code status

Historical final Figure 6A-C plotting code and the Supplementary Table S5 generator are preserved under `archived_final_code/`.

`02_validate_frozen_M11_chain.R` is a **reconstructed repository validator**, not recovered original upstream code. It validates the archived small-file chain and independently recomputes the combined Wald ratios/F statistics.

## Pleiotropy limitation

The archived GWAS Catalog table is a contextual lookup screen. Because each eligible gene retained only one independent instrument, the analysis cannot support multi-instrument heterogeneity, MR-Egger, leave-one-out, or comparable formal horizontal-pleiotropy diagnostics.

## Public-data policy

Do not redistribute the full GTEx V11 archive/parquet, full ALM GWAS summary statistics, or the full 1000 Genomes reference panel. Preserve accession/source metadata, hashes where available, frozen small intermediates, and exact analysis outputs instead.
