# Cadmium–ALMI reproducibility repository

Reproducibility materials for the manuscript integrating blood cadmium exposure, appendicular lean mass index (ALMI), skeletal-muscle transcriptomics, pathway analyses, and exploratory gene-level Mendelian randomization.

## Release status

**GitHub-ready public archive — `v1.0.0-submission` source package.**

The repository contains the frozen source/provenance materials needed to audit the manuscript's NHANES, candidate-gene, PPI/hub, cross-cohort, mitochondrial-pathway, and gene-level MR results. Large third-party public resources are not redistributed when a stable source/acquisition route and local provenance are more appropriate.

This repository distinguishes three kinds of code:

1. **Portable repository scripts** — repository-relative and intended for reuse/validation.
2. **Reconstructed scripts** — clearly labeled reconstructions created from archived analysis plans/results when original standalone code was not retained.
3. **Archived historical code** — preserved verbatim for provenance; these files may contain old Windows paths and are not presented as portable entry points.

## Repository layout

```text
00_documentation/                 study overview, parameters, status tables, manifests
01_nhanes/                        reconstructed NHANES pipeline, audits, aggregate source data
02_candidate_gene_integration/    Figure 3 candidate-gene source data and provenance
03_ppi_hub_prioritization/        STRING/PPI and cytoHubba hub-prioritization materials
04_cross_cohort_validation/       Figure 4 cross-cohort source data and portable scripts
05_mitochondrial_pathways/        Figure 5 frozen objects, source data, audit and provenance
06_mr/                            Figure 6 / Supplementary Table S5 MR source/provenance chain
07_figure_table_source_data/      compact manuscript figure/table source data
08_supplementary/                 supplementary-analysis artifacts
 environment/                     software and package-version records
 release/                         GitHub/Zenodo release guidance
```

## Key reproducibility anchors

- NHANES: 2011–2018 adults aged 20–59; final analytic N = 7,628; primary Model 3 N = 6,575.
- Candidate screening: CTD 138 genes; 3,342 nominal GSE111016 DEGs; 23 final candidates.
- Hub prioritization: 11 consensus hubs; strict hubs ATF4 and CASP9.
- Figure 5: 40 GSE111016 samples; four mitochondrial pathways; 44 hub-depleted gene–pathway correlations; 5,880 GO-BP GSEA rows archived.
- MR: 342 strict cis-eQTL candidate rows → 333 harmonized candidates → 309 EUR-reference-present variants → four independent instruments. Primary `r² < 0.001` and sensitivity `r² < 0.01` clumping retained the same four instruments.

## Quick validation

Run these from the repository root in R:

```r
# Figure 3 frozen-source anchors (base R)
source("02_candidate_gene_integration/scripts/01_validate_frozen_figure3_sources.R", encoding="UTF-8")

# Figure 5 frozen-source anchors (base R)
source("05_mitochondrial_pathways/scripts/01_validate_frozen_figure5_sources.R", encoding="UTF-8")

# Figure 6 / M11 frozen MR chain (base R)
source("06_mr/scripts/02_validate_frozen_M11_chain.R", encoding="UTF-8")
```

The Figure 4 portable plotting scripts are under `04_cross_cohort_validation/scripts/`. The reconstructed NHANES raw-XPT pipeline is under `01_nhanes/scripts/`; it requires locally acquired NHANES XPT files.

## Raw-data policy

Raw NHANES XPT files are **not redistributed**. Obtain them from the original NHANES source and either set `NHANES_RAW_XPT_DIR` or place them under `data_raw/01_raw_xpt/` with cycle subfolders. The exact 52-file SHA-256 manifest is archived at `01_nhanes/metadata/raw_xpt_manifest_sha256.tsv`.

Likewise, the full GTEx V11 eQTL archive, full ALM GWAS summary files, and the full 1000 Genomes LD reference are not redistributed in the MR module. The repository instead preserves source metadata, hashes where available, small frozen intermediates, PLINK clumping outputs, and final statistical tables.

## Known historical provenance limitations

The exact historical STRING release/access date and complete web-query settings used for Figure 3C were not retained. The repository reports this explicitly and does not infer settings from output files. See `02_candidate_gene_integration/PROVENANCE_LIMITATIONS.md`.

Some original/recovered historical scripts retain the absolute Windows paths from the original workstation. They are stored only under directories whose names begin with `archived_`; portable scripts do not depend on those paths.

## Interpretation boundaries

- NHANES is cross-sectional.
- The GEO skeletal-muscle cohorts did not measure cadmium exposure.
- Cross-layer convergence is candidate/pathway prioritization, not proof of a single cadmium → gene → muscle causal chain.
- Gene-level MR uses one independent cis-eQTL per eligible gene and is not a mediation analysis of the blood cadmium–ALMI association.

## License and third-party data

Repository-authored code and documentation are released under the MIT License. Third-party datasets, database exports, and derived materials remain subject to the terms of their original sources; inclusion here does not relicense upstream data.

## Citation

Citation metadata are provided in `CITATION.cff`. After the GitHub repository is created and Zenodo archives the tagged release, add the final repository URL/DOI in GitHub/Zenodo metadata; the source package itself contains no placeholder URL or DOI.

## Release workflow

Upload this directory to GitHub, verify the repository file list, then create a tagged release named `v1.0.0-submission`. Zenodo archiving should be performed only after that GitHub release exists. See `release/GITHUB_UPLOAD_INSTRUCTIONS.md` and `release/GITHUB_ZENODO_RELEASE_CHECKLIST.md`.
