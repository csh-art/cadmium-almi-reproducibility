# GitHub upload instructions

This package is designed so that the **contents of this directory** become the GitHub repository root. Do not upload the enclosing ZIP as the only repository file.

## GitHub Desktop / web workflow

1. Create a new empty GitHub repository, for example `cadmium-almi-reproducibility`.
2. Keep it **Private** for the first upload if you want one last visual check; switch to Public after inspection.
3. Add all files from this directory to the repository root.
4. Confirm that `README.md`, `LICENSE`, `CITATION.cff`, `.gitignore`, and folders `01_nhanes` through `06_mr` are visible.
5. Specifically confirm that these ignored-extension exceptions were uploaded:
   - `05_mitochondrial_pathways/02_analysis_files/canonical_core/M05_Figure5_CANONICAL_CORE_FINAL.RData`
   - `06_mr/provenance/clumping/MR_08_EUR_reference_overlap.log`
   - `06_mr/provenance/clumping/MR_09_PLINK_clump_r2_0.001.log`
   - `06_mr/provenance/clumping/MR_12_PLINK_clump_r2_0.01.log`
6. When satisfied, create a GitHub Release/tag named `v1.0.0-submission`.
7. Archive that tagged release with Zenodo and record the DOI.

## Command-line workflow

From inside the extracted directory:

```bash
git init
git add .
git status
git commit -m "Reproducibility repository v1.0.0-submission"
git branch -M main
git remote add origin <YOUR_GITHUB_REPOSITORY_URL>
git push -u origin main
```

The angle-bracket value above is an instruction for the local command and is intentionally not embedded as repository metadata.
