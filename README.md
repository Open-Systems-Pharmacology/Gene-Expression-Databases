# Gene-Expression-Databases

Repository hosting the latest Gene Expression Databases to be used with PK-Sim.
Databases are built from bulk RNA-Seq data provided by [Bgee](https://www.bgee.org/) (release 15.2) and cover 18 species including humans, preclinical pharmacological species, and animal health species.

Find the latest pre-built databases in the [releases section](https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/releases).

## Supported Species

| Category | Species |
|---|---|
| **Human** | Human (*Homo sapiens*) |
| **PharmaSpecies** | Mouse, Rat, Rabbit, Guinea pig, Dog, Minipig, Monkey (*M. mulatta*), Monkey (*M. fascicularis*), Monkey (Pig-tailed) |
| **AnimalHealthSpecies** | Cattle, Horse, Cat, Chicken, Goat, Sheep, Turkey, Zebrafish |

## Data and Process Flow

The repository follows a structured pipeline to build PK-Sim expression databases. Below is the end-to-end data flow:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     EXTERNAL DATA SOURCES                                   │
│  ┌──────────────────┐                           ┌──────────────────────┐   │
│  │ Ensembl BioMart  │                           │    Bgee RNA-Seq      │   │
│  │   (Annotations)  │                           │     (Expression)     │   │
│  └──────────────────┘                           └──────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘
              │                                             │
              ▼                                             ▼
┌─────────────────────────────────┐      ┌──────────────────────────────────┐
│   CODE/01-biomart-prep/         │      │   CODE/02-db-generation/         │
│   PrepareBioMarts.R             │      │   GeneratePKsimDB.R              │
│                                 │      │                                  │
│ • Download species gene annot.  │      │ • Download TPM data from Bgee    │
│ • Map NCBI / Ensembl / Symbol   │      │ • Join with BioMart annotations  │
│ • Identify ADME genes           │      │ • Normalize to PK-Sim schema     │
│ • Map human orthologs           │      │ • Create SQLite expression DB    │
└─────────────────────────────────┘      └──────────────────────────────────┘
              │                                             │
              │ Output: BioMarts/All_Species_BioMarts.DB   │
              │ (36 tables: Annotations + ADME per species) │
              │                                             │
              └─────────────────────────┬───────────────────┘
                                        │
                    ┌───────────────────▼───────────────────┐
                    │   PK-Sim DBs/{Species}/               │
                    │   GENEDB_{species}_*.expressionDB     │
                    │                                       │
                    │ Two variants per species:             │
                    │ • ADME_ONLY (ADME genes only)        │
                    │ • Full (all genes with data)         │
                    └───────────────────┬───────────────────┘
                                        │
              ┌─────────────────────────┴─────────────────────────┐
              │                                                   │
              ▼                                                   ▼
      ┌──────────────────┐                               ┌──────────────────┐
      │   CODE/05-       │                               │   CODE/04-       │
      │   utilities/     │                               │   qualification/ │
      │                  │                               │                  │
      │ • Compress to    │                               │ • Technical QC   │
      │   .tar.gz        │                               │   (Bgee vs DB)   │
      │ • Upload to      │                               │ • Biological QC  │
      │   GitHub         │                               │   (ADME profiles)│
      │   Release        │                               │                  │
      └──────────────────┘                               └──────────────────┘
              │                                                   │
              ▼                                                   ▼
      [Distribution]                                    [Validation Reports]
```

### Process Stages

1. **Setup & Configuration** (`CODE/03-helpers/`)
   - Species definitions (PharmaSpecies, AnimalHealthSpecies)
   - SQL schema and view definitions
   - Tissue-to-container mappings

2. **BioMart Annotation Preparation** (`CODE/01-biomart-prep/PrepareBioMarts.R`)
   - Downloads species-specific gene catalogs from Ensembl
   - Maps gene identifiers (Ensembl ID, NCBI ID, gene symbol, synonyms)
   - Identifies ADME-relevant genes (CYPs, UGTs, SLCs, ABCs, etc.)
   - Maps human orthologs for cross-species homology
   - Output: `BioMarts/All_Species_BioMarts.DB`

3. **Expression Database Generation** (`CODE/02-db-generation/GeneratePKsimDB.R`)
   - Downloads RNA-Seq TPM data from Bgee for each species
   - Joins TPM data with BioMart annotations
   - Normalizes to PK-Sim-compatible SQLite schema
   - Generates both ADME-only and full-genome databases
   - Output: `PK-Sim DBs/{Species}/GENEDB_{species}_*.expressionDB`

4. **Quality Control & Validation** (`CODE/04-qualification/`)
   - **Technical Qualification**: Confirms TPM values match source Bgee data
   - **Biological Qualification**: Compares ADME expression profiles across releases

5. **Distribution** (`CODE/05-utilities/`)
   - Compresses ADME databases to `.tar.gz` for GitHub releases
   - Uploads artifacts to GitHub Release pages

### Parallel Execution

Steps 1–3 run in parallel for non-human species using `foreach + doParallel`:
- All PharmaSpecies process concurrently
- All AnimalHealthSpecies process concurrently
- Human processes separately (as basis for ortholog mapping)

Each species generates 2 databases (ADME-only + full), so 18 species × 2 = 36 databases total.

---

## How to Run

### Prerequisites

- **R** ≥ 4.2 (64-bit strongly recommended — annotation processing is memory intensive)
- **Internet access** to download data from Bgee and Ensembl BioMart (large files; up to 65 GB for human)
- Sufficient disk space (~10–20 GB working space per species for intermediate files)

Install required R packages:

```r
install.packages(c(
  "here", "DBI", "RSQLite", "dplyr", "tidyr", "readr", "stringr",
  "ggplot2", "ggrepel", "scales", "foreach", "doParallel", "xlsx"
))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("BgeeDB", "biomaRt"))
```

### Execution Order

All scripts should be run from the repository root. The easiest way is to source `Code/00-pipeline/MakeAllDBs.R` which orchestrates the full pipeline:

```r
# Set the working directory to the repository root, then:
source("Code/00-pipeline/MakeAllDBs.R")
```

The pipeline runs in this order:

1. **`Code/01-biomart-prep/PrepareBioMarts.R`** — Downloads gene annotations from Ensembl BioMart for each species (Ensembl IDs, NCBI IDs, gene symbols, ortholog mappings). Human must be run first as it is the basis for homology mapping. Results are stored in `BioMarts/` as a SQLite annotation DB.

2. **`Code/02-db-generation/GeneratePKsimDB.R`** — Main function. Downloads RNA-Seq TPM data from Bgee for a given species, maps it against the BioMart annotations, and writes a PK-Sim compatible SQLite expression database to `PK-Sim DBs/{Species}/`.

3. **`Code/05-utilities/helper_compress_DBs.sh`** — Compresses each ADME-only database to `.tar.gz` for distribution.

4. **`Code/04-qualification/Qualification_BgeeDB_2_PKSimDB.R`** — Technical validation: confirms that TPM values in the generated PK-Sim DB match the raw Bgee source data for a reference experiment (GSE30611).

5. **`Code/04-qualification/Qualification_PKSimDB.R`** — Biological validation: compares relative expression profiles of key ADME proteins (CYPs, ABCs, SLCs, UGTs) between the new DB and the previous release.

### Proxy / Timeout

For large downloads (human dataset is ~65 GB), the timeout is set automatically to 3600 seconds. If you are behind a proxy, uncomment and set these lines at the top of `Code/00-pipeline/MakeAllDBs.R`:

```r
Sys.setenv("http_proxy" = "http://PROXY:PORT")
Sys.setenv("ftp_proxy" = "http://PROXY:PORT")
```

### Building a Single Species

To generate the database for one species only:

```r
library(here)
setwd(here::here())
PATH <- getwd()

source("Code/03-helpers/helper_Species.R")
source("Code/01-biomart-prep/PrepareBioMarts.R")
source("Code/02-db-generation/GeneratePKsimDB.R")

# Download gene annotations
PrepareBioMarts(SPECIE = "Rat")

# Build ADME-only database
GeneratePKsimDB(SPECIE = "Rat", ADME_ONLY = TRUE, RELEASE = "15_2", COMPUTE_IN_RAM = TRUE)

# Build full database (all genes)
GeneratePKsimDB(SPECIE = "Rat", ADME_ONLY = FALSE, RELEASE = "15_2", COMPUTE_IN_RAM = TRUE)
```

Valid values for `SPECIE`: `Human`, `Mouse`, `Rat`, `Rabbit`, `Guineapig`, `Dog`, `Minipig`, `Monkey_mulatta`, `Monkey_fascicularis`, `Monkey_PigTailed`, `Cattle`, `Horse`, `Cat`, `Chicken`, `Goat`, `Sheep`, `Turkey`, `Zebrafish`

## Database Contents

Each species produces two SQLite files (distributed as `.tar.gz`):

| File | Contents |
|---|---|
| `GENEDB_{species}_ADME_ONLY_BgeeRelease_15_2.expressionDB` | ADME-relevant genes only (CYPs, UGTs, SLCs, ABCs, etc.) |
| `GENEDB_{species}_BgeeRelease_15_2.expressionDB` | All genes with RNA-Seq data |

### Large File Strategy (GitHub Release Assets)

Some compressed databases can exceed GitHub's per-file push limit for git history. In those cases, publish the `.tar.gz` as a GitHub Release asset instead of committing it to the branch.

Current case:

- `PK-Sim DBs/Mouse/GENEDB_mouse_ADME_ONLY_BgeeRelease_15_2.expressionDB.tar.gz` is release-asset only.

Use the helper script to create/update a release and upload the asset:

```bash
export GH_TOKEN="<github-token-with-repo-scope>"
Code/05-utilities/helper_upload_release_asset.sh \
  v3.0.2 \
  "PK-Sim DBs/Mouse/GENEDB_mouse_ADME_ONLY_BgeeRelease_15_2.expressionDB.tar.gz" \
  "OSP Expression DB v3.0.2"
```

The script creates a draft release if the tag does not already exist, then uploads the file as an asset.

### Database Schema

The SQLite databases use the following key tables:

| Table | Description |
|---|---|
| `tab_genes` | One row per gene; primary key `gene_id` |
| `tab_gene_variants` | Gene isoforms/variants; links `gene_id` → `variant_id` |
| `tab_gene_names` | All name types per gene (columns: `gene_id`, `gene_name`, `name_type`). `name_type` values include `SYMBOL`, `ENSEMBL`, `GENE_ID` (NCBI), `OFFICIAL_FULL_NAME`, `SYNONYM`, `HOMOLOG_SYMBOL` |
| `tab_expression_data_values` | TPM expression values; links to `variant_id` and `data_source_id` |
| `tab_expression_data_records` | One record per RNA-Seq library; stores tissue, gender, health state, age, sample source |
| `tab_expression_data_properties` | Key–value metadata per record (e.g. `data_base`, experiment ID) |
| `tab_expression_data_ages` | Age ranges associated with records |
| `tab_container_tissue` | Maps Bgee tissue ontology terms to PK-Sim organ/container names |
| `tab_dts_properties` | Column metadata used by PK-Sim UI (labels, tooltips, dimensions) |

### Querying the Database

Helper functions are provided in `Code/03-helpers/helper_SQL_Queries.R`:

```r
library(DBI); library(RSQLite); library(dplyr)

conn <- DBI::dbConnect(RSQLite::SQLite(), "PK-Sim DBs/Rat/GENEDB_rat_ADME_ONLY_BgeeRelease_15_2.expressionDB")

source("Code/03-helpers/helper_SQL_Queries.R")

# Find a gene by symbol, synonym, or other identifier
genes <- get_proteins_by_name(name = c("Cyp3a2", "Cyp1a1"), conn = conn)

# Retrieve expression values for a gene
expr <- get_expression_data_by_gene_id(
  P_ID   = genes |> dplyr::filter(has_data == 1) |> dplyr::pull(gene_id),
  conn   = conn
)

DBI::dbDisconnect(conn)
```

## Code Organization

The `Code/` folder is organized by processing stage for clarity and maintainability:

```
Code/
├── 00-pipeline/              # Main orchestration
│   └── MakeAllDBs.R          # Pipeline entry point (runs all species)
│
├── 01-biomart-prep/          # BioMart annotation preparation
│   └── PrepareBioMarts.R     # Download and prepare gene annotations
│
├── 02-db-generation/         # Core expression database generation
│   └── GeneratePKsimDB.R     # Main DB generation function
│
├── 03-helpers/               # Core utilities and configuration
│   ├── helper_Species.R      # Defines species lists and constants
│   ├── helper_SQL_Commands.R # Table DDL, indices, view definitions
│   ├── helper_SQL_Queries.R  # Query helpers for the databases
│   ├── helper_All_Bgee_organs.R  # Organ/age extraction helpers
│   ├── helper_Relative_Expression.R  # Expression profile utilities
│   ├── tab_container_tissue.txt      # Tissue-to-container mapping
│   └── tab_dts_properties.txt        # PK-Sim UI column metadata
│
├── 04-qualification/         # Quality control and validation
│   ├── Qualification_BgeeDB_2_PKSimDB.R  # Technical validation
│   └── Qualification_PKSimDB.R           # Biological validation
│
├── 05-utilities/             # Utilities and maintenance
│   ├── helper_compress_DBs.sh    # Compress DBs to .tar.gz
│   ├── helper_upload_release_asset.sh  # Upload to GitHub releases
│   ├── helper_update_container_mapping.R   # Update tissue mappings
│   ├── helper_update_indizes.R.r           # Index maintenance
│   ├── helper_SQL_Queries.cs              # C# SQL helpers (legacy)
│   ├── helper_compress_DBs_7zip.bat       # Windows compression alt
│   └── Qualification/              # Additional QC data
```

## Code of conduct
Everyone interacting in the Open Systems Pharmacology community (codebases, issue trackers, chat rooms, mailing lists etc...) is expected to follow the Open Systems Pharmacology [code of conduct](https://github.com/Open-Systems-Pharmacology/Suite/blob/master/CODE_OF_CONDUCT.md).

## Contribution
We encourage contribution to the Open Systems Pharmacology community. Before getting started please read the [contribution guidelines](https://github.com/Open-Systems-Pharmacology/Suite/blob/master/CONTRIBUTING.md). If you are contributing code, please be familiar with the [coding standards](https://github.com/Open-Systems-Pharmacology/Suite/blob/master/CODING_STANDARDS.md).

## License
Gene-Expression-Databases is released under the [GPLv2 License](LICENSE).

All trademarks within this document belong to their legitimate owners.
