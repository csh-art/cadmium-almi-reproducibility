# 03 — PPI and hub-gene prioritization

## Figure 3C PPI

`figure3_ppi/` contains the archived Cytoscape session, vector SVG, original 23-gene STRING query list, full and de-duplicated interaction exports, and the retained 22-node degree table.

The manuscript-safe interpretation is: **22 of the 23 candidate genes were retained in the archived final STRING-derived network; `VEGFA` is absent from the retained node-degree table.** The repository does not infer why `VEGFA` was absent.

The exact historical STRING release/access date and complete web-query settings were not retained. This provenance limitation is disclosed rather than guessed.

## Figure 4 hub prioritization

`figure4_hub_prioritization/` contains the consensus source table, locked 11-hub discovery results, Figure 4A-C source data, and final plotting-code snapshots. Legacy Degree/MCC/EPC top-10 exports are retained under `audit/legacy_cytoHubba_top10/`.

The 11 consensus hubs are defined by support from at least 2 of 3 cytoHubba algorithms; ATF4 and CASP9 are the two 3/3 strict hubs.
