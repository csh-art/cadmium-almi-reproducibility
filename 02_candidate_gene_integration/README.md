# 02 — Cadmium–sarcopenia candidate-gene integration (Figure 3)

## Public-release status

The frozen Figure 3 evidence chain is archived and suitable for public audit. The repository contains the exact historical GSE111016 count/phenotype inputs used for the reconstructed screening, complete DESeq2 results, threshold audits, the final 23-gene list, enrichment source tables, compact expression matrices, and STRING source exports.

The 156 MB historical DESeq2 `dds` object is intentionally not included because it exceeds GitHub's normal per-file limit and is regenerable from the archived counts and phenotype inputs.

## Locked screening logic

1. Historical CTD cadmium interaction count: `interaction_count > 10` → **138 genes**.
2. GSE111016: **65,217 rows**, **40 samples** (20 control / 20 sarcopenia-related).
3. DESeq2 design: `~ sarcopenia`, `yes` versus `no`; no low-expression prefilter.
4. Nominal `p < 0.05` → **3,342** results (1,771 up / 1,571 down).
5. Candidate prioritization additionally required `|log2FC| > 0.1`.
6. Mapping to official symbols and intersecting with CTD 138 → **23 candidate genes**.
7. `SLC11A2` was excluded automatically because `|log2FC| = 0.0998187`.

## Code organization

- `scripts/01_validate_frozen_figure3_sources.R`: portable repository validator.
- `archived_original_code/`: original Figure 3A and panel/composite code snapshots. These retain historical Windows paths for provenance and are not advertised as portable entry points.

## STRING limitation

Exact historical STRING release/access-date and full query settings were not retained. This is documented explicitly in `PROVENANCE_LIMITATIONS.md`; no undocumented setting is inferred from the archived output.
