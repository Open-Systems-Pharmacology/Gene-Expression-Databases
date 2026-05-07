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

All scripts should be run from the repository root. The easiest way is to source `Code/MakeAllDBs.R` which orchestrates the full pipeline:

```r
# Set the working directory to the repository root, then:
source("Code/MakeAllDBs.R")
```

The pipeline runs in this order:

1. **`Code/PrepareBioMarts.R`** — Downloads gene annotations from Ensembl BioMart for each species (Ensembl IDs, NCBI IDs, gene symbols, ortholog mappings). Human must be run first as it is the basis for homology mapping. Results are stored in `BioMarts/` as a SQLite annotation DB.

2. **`Code/GeneratePKsimDB.R`** — Main function. Downloads RNA-Seq TPM data from Bgee for a given species, maps it against the BioMart annotations, and writes a PK-Sim compatible SQLite expression database to `PK-Sim DBs/{Species}/`.

3. **`Code/helper_compress_DBs.sh`** — Compresses each ADME-only database to `.tar.gz` for distribution.

4. **`Code/Qualification_BgeeDB_2_PKSimDB.R`** — Technical validation: confirms that TPM values in the generated PK-Sim DB match the raw Bgee source data for a reference experiment (GSE30611).

5. **`Code/Qualification_PKSimDB.R`** — Biological validation: compares relative expression profiles of key ADME proteins (CYPs, ABCs, SLCs, UGTs) between the new DB and the previous release.

### Proxy / Timeout

For large downloads (human dataset is ~65 GB), the timeout is set automatically to 3600 seconds. If you are behind a proxy, uncomment and set these lines at the top of `MakeAllDBs.R`:

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

source("Code/PrepareBioMarts.R")
source("Code/GeneratePKsimDB.R")
source("Code/helper_Species.R")

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
Code/helper_upload_release_asset.sh \
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

Helper functions are provided in `Code/helper_SQL_Queries.R`:

```r
library(DBI); library(RSQLite); library(dplyr)

conn <- DBI::dbConnect(RSQLite::SQLite(), "PK-Sim DBs/Rat/GENEDB_rat_ADME_ONLY_BgeeRelease_15_2.expressionDB")

source("Code/helper_SQL_Queries.R")

# Find a gene by symbol, synonym, or other identifier
genes <- get_proteins_by_name(name = c("Cyp3a2", "Cyp1a1"), conn = conn)

# Retrieve expression values for a gene
expr <- get_expression_data_by_gene_id(
  P_ID   = genes |> dplyr::filter(has_data == 1) |> dplyr::pull(gene_id),
  conn   = conn
)

DBI::dbDisconnect(conn)
```

## Code Structure

| File | Purpose |
|---|---|
| `Code/MakeAllDBs.R` | Pipeline entry point — runs all species in parallel |
| `Code/GeneratePKsimDB.R` | Core DB generation function |
| `Code/PrepareBioMarts.R` | Gene annotation download and preparation |
| `Code/helper_Species.R` | Defines species lists (`PharmaSpecies`, `AnimalHealthSpecies`, `ALL_SPECIE`) |
| `Code/helper_SQL_Commands.R` | Table DDL, index definitions, SQL views |
| `Code/helper_SQL_Queries.R` | Query helper functions for the output databases |
| `Code/helper_All_Bgee_organs.R` | Organ/age export helpers |
| `Code/helper_Relative_Expression.R` | Utility to extract relative expression by Ensembl ID |
| `Code/helper_update_container_mapping.R` | Updates tissue-to-container mapping table |
| `Code/helper_compress_DBs.sh` | Compresses ADME-only DBs to `.tar.gz` |
| `Code/Qualification_BgeeDB_2_PKSimDB.R` | Validates Bgee source data == PK-Sim DB data |
| `Code/Qualification_PKSimDB.R` | Compares ADME expression profiles old vs new DB |
| `Code/tab_container_tissue.txt` | Static tissue-to-container mapping lookup |
| `Code/tab_dts_properties.txt` | Static column metadata for PK-Sim UI |

## Code of conduct
Everyone interacting in the Open Systems Pharmacology community (codebases, issue trackers, chat rooms, mailing lists etc...) is expected to follow the Open Systems Pharmacology [code of conduct](https://github.com/Open-Systems-Pharmacology/Suite/blob/master/CODE_OF_CONDUCT.md).

## Contribution
We encourage contribution to the Open Systems Pharmacology community. Before getting started please read the [contribution guidelines](https://github.com/Open-Systems-Pharmacology/Suite/blob/master/CONTRIBUTING.md). If you are contributing code, please be familiar with the [coding standards](https://github.com/Open-Systems-Pharmacology/Suite/blob/master/CODING_STANDARDS.md).

## License
Gene-Expression-Databases is released under the [GPLv2 License](LICENSE).

All trademarks within this document belong to their legitimate owners.
