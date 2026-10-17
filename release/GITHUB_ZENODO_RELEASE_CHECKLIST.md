# GitHub + Zenodo release checklist

## GitHub upload readiness — completed in v13

- [x] Six analytical modules populated with source/provenance materials.
- [x] Release placeholders removed from `CITATION.cff` and repository documentation.
- [x] MIT code/documentation license added.
- [x] Public executable R scripts use repository-relative paths or documented environment variables.
- [x] Historical absolute-path code separated into clearly named `archived_*` provenance directories.
- [x] No participant-level NHANES dataset is included.
- [x] No full GTEx V11, ALM GWAS, or 1000 Genomes reference dataset is included.
- [x] GitHub 100 MB/file limit checked; no repository file exceeds it.
- [x] `.gitignore` exceptions added so the canonical Figure 5 `.RData` and scientific PLINK `.log` files are not accidentally omitted.
- [x] Local user-home strings sanitized from environment records.
- [x] Known STRING metadata limitation documented without guessing.
- [x] Public-release static audit recorded in `GITHUB_READY_AUDIT_v13.md`.

## After creating the GitHub repository

- [ ] Push the `v13_GitHub_READY` contents to the repository root.
- [ ] Confirm the canonical Figure 5 `.RData` and three MR `.log` files are visible on GitHub.
- [ ] Optionally add the final GitHub repository URL to `CITATION.cff` as `repository-code`.
- [ ] Create Git tag/release `v1.0.0-submission`.
- [ ] Connect/enable Zenodo for the repository and archive that release.
- [ ] Add the Zenodo DOI to repository metadata and the manuscript Code Availability statement.

## Maintenance recommendation (not a blocker to GitHub upload)

A fresh raw-to-result rerun in the recorded software environment is recommended before future methodological extensions. The present public release distinguishes frozen-source validation, reconstructed code, and archived historical code rather than claiming that every historical interactive step was rerun during repository packaging.
