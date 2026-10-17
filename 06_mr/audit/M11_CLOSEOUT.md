# M11 / Figure 6 close-out

The small frozen MR_05–MR_11 intermediates requested during the repository audit are now archived.

The preserved chain is:

`MR_05 strict cis-eQTL candidates (342)`
→ `MR_06 exact combined-ALM matches (327) + 15 initially missing`
→ `MR_06c position recovery (6)`
→ `MR_07 harmonized candidates (333)`
→ `MR_09 EUR-reference clump mapping (309; 24 absent)`
→ `MR_10 four independent instruments`
→ combined/sex-stratified Wald-ratio tables, F statistics, and Steiger directionality.

Primary PLINK clumping (`r2 < 0.001`, 10 Mb) and the `r2 < 0.01` sensitivity retain the same four instruments.

The source/provenance evidence chain is therefore closed for the public archive. Large public GTEx V11, ALM GWAS, and 1000 Genomes files are intentionally not redistributed. The original analysis was built interactively; the repository includes archived final plotting/table code and a **reconstructed** base-R validator for the frozen small-file chain. The reconstructed validator is not presented as recovered original upstream code.
