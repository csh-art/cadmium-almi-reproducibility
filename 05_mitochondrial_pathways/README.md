# M05 — Mitochondrial pathway analysis (Figure 5)

## Public-release status

The frozen Figure 5 source/provenance evidence chain is archived and audited. The repository contains the canonical input/core objects, final plotting snapshots, recovered historical upstream functions, audit reports, and human-readable Figure 5A-D source-data exports.

The 2026-10-17 audit confirmed 25 required core objects, 40 samples (20 control / 20 sarcopenia-related), four locked mitochondrial pathways, 11 consensus hubs, two strict hubs (`ATF4`, `CASP9`), 44 Figure 5C correlations, complete hub removal from Figure 5C pathway memberships, and 5,880 GO-BP GSEA rows.

## Code organization

- `scripts/01_validate_frozen_figure5_sources.R`: portable validator for the human-readable frozen source data.
- `archived_recovered_functions/`: seven functions recovered from the historical R workspace. These preserve historical absolute Windows paths and are retained verbatim as provenance, not as portable entry points.
- `archived_original_code/`: final Figure 5A-D plotting-code snapshots.
- `archived_audit_code/`: the local audit/export script used during repository reconstruction.

The original 2-byte failed reconstruction-code capture is retained under `historical_incomplete_capture/` to preserve provenance.

## Locked anchors

- 40 GSE111016 samples: 20 control and 20 sarcopenia-related.
- 11 consensus hubs; strict hubs: `ATF4`, `CASP9`.
- Detailed pathways: oxidative phosphorylation, respiratory electron transport chain, mitochondrial ATP synthesis, mitochondrial translation.
- Figure 5C: 44 gene–pathway correlations after simultaneously removing all 11 consensus hub genes from the relevant pathway gene sets.

The frozen evidence chain can be audited from the archived source data. A fresh end-to-end rerun in a new environment is recommended for future maintenance but is not represented as having been performed inside this release build.
