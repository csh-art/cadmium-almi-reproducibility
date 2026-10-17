# GitHub-ready static release audit — v13

Generated: 2026-10-17

| Check | Status |
|---|---|
| CITATION.cff parses and has 5 authors + MIT | PASS |
| Release placeholders absent | PASS |
| No email addresses embedded in public text/code | PASS |
| No local C:/Users/... strings remain | PASS |
| No absolute Windows paths in public executable R scripts | PASS |
| No auto-install calls in public executable R scripts | PASS |
| No file exceeds GitHub 100 MB limit | PASS |
| Figure 5 canonical RData explicitly unignored | PASS |
| Scientific PLINK logs explicitly unignored | PASS |
| Figure 3 frozen source count anchors pass | PASS |
| Figure 5 frozen source count anchors pass | PASS |
| M11 final validation report present | PASS |

- Largest file: `05_mitochondrial_pathways/source_data/06_Figure5D_running_ES.csv` (8.08 MB).
- Files >100 MB: 0.
- Public executable R scripts with absolute Windows paths: 0.
- Public executable R scripts with package auto-install calls: 0.

## Scope note
Historical archived/provenance code may intentionally retain original workstation paths or historical package-install suggestions. Those files are not public execution entry points and are labeled by directory/documentation.

## Known non-blocking provenance limitation
Exact historical STRING release/access date and complete Figure 3C web-query settings were not retained. This is explicitly documented rather than reconstructed by assumption.

## Verdict
**READY FOR GITHUB UPLOAD**
