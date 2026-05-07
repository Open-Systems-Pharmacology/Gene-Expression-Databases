# System Patterns — Gene-Expression-Databases

## Architecture Overview

```
External Sources (Ensembl BioMart, Bgee RNA-Seq)
    ↓
Stage 01: PrepareBioMarts → BioMarts/All_Species_BioMarts.DB (36 annotation tables)
    ↓
Stage 02: GeneratePKsimDB → PK-Sim DBs/{Species}/*.expressionDB (32 tables per DB)
    ├─ Parallel: PharmaSpecies (9 species via PSOCK cluster)
    ├─ Parallel: AnimalHealthSpecies (8 species via PSOCK cluster)
    └─ Sequential: Human (first, basis for ortholog mapping)
    ↓
Stage 04: Qualify → Validation plots & CSV exports (3-level QC framework)
    ↓
Stage 05: Compress & Upload → GitHub Release assets (*.tar.gz)
```

## Code Structure (`/Code/` — 6-stage pipeline)

```
Code/
├── 00-pipeline/           # MakeAllDBs.R — main orchestration
├── 01-biomart-prep/       # PrepareBioMarts.R — BioMart annotation download
├── 02-db-generation/      # GeneratePKsimDB.R — core DB generation
├── 03-helpers/            # Utilities: helper_*.R, SQL commands, static lookup tables
├── 04-qualification/      # QC: 3-level validation scripts + plots + config
└── 05-utilities/          # Maintenance: compress, upload, index updates
```

## Key Design Patterns

### Resource Safety — on.exit() Guards
All DB connections and working-directory changes use `on.exit()` for cleanup:
```r
old_wd <- setwd(new_path)
on.exit(setwd(old_wd), add = TRUE)

conn <- DBI::dbConnect(...)
on.exit(DBI::dbDisconnect(conn), add = TRUE)

old_timeout <- getOption("timeout")
options(timeout = 60 * 60 * 60)
on.exit(options(timeout = old_timeout), add = TRUE)
```

### Parallel Processing — PSOCK Clusters
- Species processed in parallel via `foreach` + `doParallel`
- Workers require explicit `.export` (free variables) and `.packages` (libraries)
- `seq_along()` always used (not `1:length()`) to handle empty vectors safely

### SQL Safety — Parameterized Queries
- **Never** use `sprintf()` for SQL string interpolation
- Always use `DBI::dbGetQuery(conn, sql, params = list(...))` with `?` placeholders
- IN clauses: `paste(rep("?", n), collapse = ",")` for variable-length lists

### Helper Wrapping — No Side Effects on Source
Files intended to be `source()`d must wrap executable code in functions:
```r
# helper_All_Bgee_organs.R
build_bgee_lookup_tables <- function(PATH) { ... }
# MakeAllDBs.R calls it explicitly:
build_bgee_lookup_tables(PATH)
```

### Species Name Format
- BgeeDB requires underscore-delimited format: `Canis_lupus_familiaris`
- Never use spaces (causes silent BgeeDB lookup failure)
- BioMart dataset IDs use short prefix format: `mnemestrina_gene_ensembl` (NOT `nnemestrina`)

### Error Signaling
- `stop()` for errors (not `simpleError()`)
- `message()` for informational output (not `simpleMessage()`)
- `switch()` blocks always include `default = stop(paste("Unknown SPECIE:", SPECIE))`

### Typed NAs in dplyr
- Use `NA_real_` (not `NA`) for numeric columns to avoid type coercion in `bind_rows()`
- Use `NA_character_` for character columns

### Regex Safety
- `stringr::str_escape(name)` before building regex patterns from gene names
- Gene names can contain regex metacharacters (`.`, `*`, `(`, `)`)

## Qualification Framework (3-Level)

```
Code/04-qualification/
├── 01_config/              # Shared configs (proteins-validation.txt, container-mapping.txt, cross-species-genes.txt)
├── level1-technical-validation/   # Bgee source == OSP DB integrity checks
│   └── Qualification_BgeeDB_2_PKSimDB.R
├── level2-human-old-vs-new/       # Old RT-PCR vs New Bgee ADME profiles
│   └── Qualification_PKSimDB.R
├── level3-cross-species/          # 40-gene ADME panel across species
│   └── Qualification_CrossSpecies.R
├── level1/ level2/ level3/        # Output dirs (02_data/, 03_plots/)
├── DESCRIPTION                    # R package spec for reproducibility
└── renv.lock                      # Pinned package versions
```

### Level 1: Technical Validation
- Compares Bgee GSE30611 ERX011211 sample data against OSP expressionDB
- Generates 7 X-Y scatter plots by ADME family with gene labels (ggrepel)
- Two scenarios: OSP DB available (TPM vs sample count) or unavailable (TPM vs detection flag)

### Level 2: Biological Validation
- Compares Old DB (RT-PCR) vs New DB (Bgee TPM) for ADME genes
- Generates 4 composite violin plots by family (CYP, ABC, UGT, SLC)
- Gene symbol mapping: `variant_name <- if_else(!is.na(gene_name) & nzchar(gene_name), gene_name, variant_name)`
- Exports `validation_summary.csv` with correlation statistics

### Level 3: Cross-Species Validation
- 40 ADME genes compared across key preclinical species
- One plot per gene, showing relative expression by tissue

## SQLite Database Schema
- **32 tables per DB**: tab_genes, tab_gene_variants, tab_gene_names, tab_expression_data_values, tab_expression_data_records, tab_expression_data_properties, tab_container_tissue, tab_dts_properties, + views + indices
- **Two DB variants**: ADME_ONLY (faster, smaller) and Full (all genes)
- **Field types**: `total_count` = `"real"` (DOUBLE) to preserve TPM decimal precision

## File Size Management
- GitHub per-file hard limit: 100 MB
- Mouse ADME archive (~113 MB) excluded from git via `.gitignore`
- Uploaded via `Code/05-utilities/helper_upload_release_asset.sh` using GitHub Releases API + `GH_TOKEN`
- `.archive/` folder excluded from git (stores old qualification plots for reference)
