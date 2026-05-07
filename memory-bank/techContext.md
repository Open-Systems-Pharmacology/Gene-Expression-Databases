# Tech Context — Gene-Expression-Databases

## Language & Runtime
- **R** version ≥ 4.4.1 (2024-06-14)
- **RScript** for batch execution
- **Bash** for utility shell scripts

## R Package Dependencies (pinned in DESCRIPTION + renv.lock)

| Package | Version | Purpose |
|---------|---------|---------|
| here | 1.0.1 | Reproducible path resolution |
| dplyr | 1.2.1 | Data manipulation (lazy eval via tbl/collect) |
| readr | 2.1.5 | CSV reading |
| ggplot2 | 4.0.0 | Visualization (violin, XY scatter) |
| ggrepel | 0.9.6 | Non-overlapping gene labels on plots |
| scales | 1.4.0 | Color/scale utilities |
| tidyr | 1.3.1 | Data reshaping |
| DBI | 1.2.3 | Database interface (parameterized queries) |
| RSQLite | 2.4.3 | SQLite driver for `.expressionDB` files |
| BgeeDB | 2.32.0 | Bgee release 15.2 data access (Bioconductor) |
| biomaRt | 2.62.0 | Ensembl BioMart annotation download (Bioconductor) |
| stringr | 1.5.1 | String utilities (str_escape for regex safety) |
| tidyselect | 1.2.1 | Column selection helpers |
| rlang | 1.1.4 | R language utilities |
| foreach | 1.5.2 | Parallel iteration |
| doParallel | 1.0.17 | PSOCK parallel backend |

**Repositories**: CRAN + Bioconductor  
**Environment files**: `Code/04-qualification/DESCRIPTION` + `Code/04-qualification/renv.lock`

## External Data Services
- **Bgee** (bgee.org): RNA-Seq expression data, release 15.2 (2024-05-21)
  - Accessed via BgeeDB R package
  - Pre-cached locally as `BgeeDBs/{Species}_bgee.db` SQLite files
  - Download timeout: 60×60×60 seconds (large human dataset ~65 GB)
- **Ensembl BioMart** (ensembl.org): Gene annotations
  - Accessed via biomaRt R package
  - Cached locally in `BioMarts/All_Species_BioMarts.DB`

## Output Formats
- **`.expressionDB`**: SQLite database, PK-Sim compatible, 32 tables per species
- **`.expressionDB.tar.gz`**: Compressed archive for distribution (~50% size reduction)
- **PNG (300 DPI)**: Publication-ready qualification plots
- **CSV**: Intermediate data and validation summary exports

## Development Setup
```bash
# Clone repo
git clone https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases.git
cd Gene-Expression-Databases

# (Optional) Restore R environment
cd Code/04-qualification
Rscript -e "renv::restore()"

# Regenerate all databases
Rscript Code/00-pipeline/MakeAllDBs.R
```

## CI/CD & GitHub Actions
- Shared ARC runners: `runs-on: ['atmos-aws-arc-runner-set']`
- No automated pipeline currently; manual execution model
- PR validation via CodeRabbit automated review bot

## Git Conventions
- Branch naming: `{issue-number}-{description}-{version}`
- Commit style: conventional commits (`fix:`, `feat:`, `refactor:`, `chore:`, `enhance:`)
- PR target: master

## Key File Paths

| Path | Description |
|------|-------------|
| `Code/00-pipeline/MakeAllDBs.R` | Main orchestration entrypoint |
| `Code/01-biomart-prep/PrepareBioMarts.R` | BioMart annotation builder |
| `Code/02-db-generation/GeneratePKsimDB.R` | Core DB generation function |
| `Code/03-helpers/helper_Species.R` | Species name/category definitions |
| `Code/03-helpers/helper_SQL_Commands.R` | SQLite DDL + view definitions |
| `Code/03-helpers/helper_SQL_Queries.R` | Parameterized query helpers |
| `Code/03-helpers/helper_Relative_Expression.R` | Ensembl ID lookup utility |
| `Code/03-helpers/helper_All_Bgee_organs.R` | Organ/age lookup table builder |
| `Code/04-qualification/01_config/` | Shared config TXT files |
| `Code/04-qualification/level1-technical-validation/` | Level 1 QC script |
| `Code/04-qualification/level2-human-old-vs-new/` | Level 2 QC script |
| `Code/04-qualification/level3-cross-species/` | Level 3 QC script |
| `Code/05-utilities/helper_compress_DBs.sh` | Compress ADME DBs to .tar.gz |
| `Code/05-utilities/helper_upload_release_asset.sh` | Upload large files to GitHub Releases |
| `BgeeDBs/BgeeAges.txt` | Life stage lookup table (age → Bgee stage ID) |
| `BgeeDBs/BgeeOrgans.txt` | Organ/tissue lookup table |

## Technical Constraints
- **COMPUTE_IN_RAM=TRUE** loads full Human dataset (~65 GB) — requires 64-bit R + large memory
- **Mouse ADME archive** (~113 MB) exceeds GitHub's 100 MB per-file hard limit — must be uploaded as release asset only
- **BgeeDB species format**: Must use underscores (`Canis_lupus_familiaris`), not spaces
- **BioMart dataset IDs**: Short prefix format (e.g., `mnemestrina_gene_ensembl` for Pig-tailed macaque)
- **dplyr**: Requires bare column names in `arrange()`, not strings
- **PSOCK workers**: Must explicitly `.export` variables and `.packages` libraries
