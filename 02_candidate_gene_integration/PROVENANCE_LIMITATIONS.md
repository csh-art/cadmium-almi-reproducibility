# Figure 3 provenance notes and known historical limitations

The Figure 3 source-data chain is archived at the level needed to audit the manuscript results: CTD screening, complete GSE111016 DESeq2 results, the final 23-gene candidate list, GO/KEGG tables, compact expression matrices, and the STRING interaction/node-degree exports are present.

## Historical STRING metadata not recoverable

The exact STRING database release/access date and the complete historical web-query settings were not retained. The archived interaction table has a minimum `combined_score` of 0.400, but the repository does **not** infer the original query threshold solely from that observed minimum. This limitation is reported rather than reconstructed from guesswork.

The retained source tables directly support the manuscript-safe statement that 22 of the 23 candidate genes were present in the archived final network and that `VEGFA` is absent from the retained node-degree table. They do not establish why `VEGFA` was absent.

## Code provenance

Original Figure 3A and panel/composite plotting-code snapshots are stored under `archived_original_code/`. They preserve historical local Windows paths and are intentionally **not** presented as portable entry points.

For repository-level validation without those local paths, use:

```r
source("02_candidate_gene_integration/scripts/01_validate_frozen_figure3_sources.R", encoding="UTF-8")
```

The validator checks the frozen source-data anchors used in the manuscript; it does not claim to reproduce unrecoverable historical STRING web-session metadata.
