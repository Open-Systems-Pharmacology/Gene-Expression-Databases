# Software and Package Versions

This repository is built with R scripts and SQLite outputs. The versions below were captured from the active Linux build environment used to validate the Olive baboon integration for Bgee release 15.2.

## Core Software

| Component | Version |
| --- | --- |
| OS | Linux 5.14.0-570.21.1.el9_6.x86_64 |
| bash | 5.1.8(1)-release |
| git | 2.47.1 |
| R | 4.4.1 |
| R platform | x86_64-pc-linux-gnu |
| RSQLite | 2.4.3 |

## R Packages

| Package | Version |
| --- | --- |
| BgeeDB | 2.32.0 |
| biomaRt | 2.62.1 |
| DBI | 1.2.3 |
| RSQLite | 2.4.3 |
| dplyr | 1.2.1 |
| tidyr | 1.3.1 |
| readr | 2.1.5 |
| stringr | 1.5.2 |
| ggplot2 | 4.0.0 |
| ggrepel | 0.9.6 |
| scales | 1.4.0 |
| foreach | 1.5.2 |
| doParallel | 1.0.17 |
| here | 1.0.1 |
| xlsx | 0.6.5 |
| rlang | 1.2.0 |
| tibble | 3.3.0 |
| tidyselect | 1.2.1 |

## Validated Olive Baboon Build

The following species-specific artifacts were validated in this environment:

- BioMart dataset: `panubis_gene_ensembl`
- Bgee dataset folder: `BgeeDBs/Papio_anubis_Bgee_15_2`
- BioMart tables: `Baboon_olive_Annotations`, `Baboon_olive_ADME`
- Generated PK-Sim DB: `PK-Sim DBs/Baboon_olive/GENEDB_baboon_olive_ADME_ONLY_BgeeRelease_15_2.expressionDB`
- Generated DB size: 46 MB
- Generated DB content summary:
  - 32 SQLite tables
  - 1,058 genes
  - 1,058 gene variants
  - 553,763 expression values

## Notes

- This document captures the environment used for a validated build, not an enforced lockfile.
- If package versions change, update this file after rerunning the relevant build and validation commands.
