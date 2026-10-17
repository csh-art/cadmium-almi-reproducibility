# Recovered Figure 5 functions

These R functions were recovered on 2026-10-17 from the private historical
`project_workspace.RData` by the repository audit script. They are preserved as
historical executable provenance, not rewritten substitutes.

Most important upstream function:
- `run_figure5_abcd_final.R`: reconstructs the Figure 5 analysis from the
  canonical GSE111016 expression/group input and the stored GO-BP GSEA object,
  generates Figure 5A-D analysis objects, and freezes
  `M05_Figure5_CANONICAL_CORE_FINAL.RData`.
- `run_ssgsea.R`: the recovered ssGSEA helper.
- `rebuild_figure5_bcd.R`: regenerates final Figure 5B-D panels from the frozen
  canonical core.

The recovered code contains original Windows absolute paths. It is retained
verbatim for provenance. Repository-relative portability cleanup should be done
in a separate script rather than silently editing these historical files.
