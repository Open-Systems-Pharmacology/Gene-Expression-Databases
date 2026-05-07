# Progress — Gene-Expression-Databases

## What Works ✅

### Core Pipeline (Fully Functional)
- **PrepareBioMarts.R**: Downloads Ensembl BioMart annotations for all 18 species into `BioMarts/All_Species_BioMarts.DB`
- **GeneratePKsimDB.R**: Generates ADME-only and full SQLite `.expressionDB` files for all 18 species
- **MakeAllDBs.R**: Orchestrates parallel (PharmaSpecies + AnimalHealthSpecies) + sequential (Human) pipeline
- **Parallel processing**: PSOCK clusters with correct `.export` and `.packages` declarations

### Validation (All 3 Levels Tested ✅)
- **Level 1 — Technical**: 7 X-Y scatter plots by ADME family (Bgee TPM vs OSP sample count)
- **Level 2 — Biological**: 4 composite violin plots (Old RT-PCR vs New Bgee) with gene-symbol mapping fix
- **Level 3 — Cross-species**: 40-gene ADME panel plots across preclinical species

### Data Quality
- `BgeeAges.txt`: All life stage names correct (TAGE → STAGE typo fixed)
- `BioMart dataset IDs`: All 18 species correct (nnemestrina → mnemestrina fixed)
- Species naming: All use underscore format consistently across all files

### Code Quality (PR #7 CodeRabbit Rounds 1 & 2 Complete)
- Resource safety: `on.exit()` guards for all DB connections, setwd, global options
- SQL injection prevention: Parameterized `DBI::dbGetQuery(..., params = list(...))`
- Loop safety: `seq_along()` instead of `1:length()`
- Error signaling: `stop()`/`message()` not `simpleError()`/`simpleMessage()`
- Regex safety: `stringr::str_escape()` for gene name patterns
- Typed NAs: `NA_real_` for numeric columns in dplyr bind operations
- Switch safety: `default = stop()` in all species switch blocks
- Case consistency: `Guineapig` (not `GuineaPig`) matching helper_Species.R

### Documentation
- README: Updated with ASCII pipeline diagram, all species, stage descriptions
- DESCRIPTION + renv.lock: R environment fully specified and pinned
- copilot-instructions.md: Mermaid diagrams properly fenced

### File Organization
```
Code/04-qualification/
├── 01_config/             # proteins-validation.txt, container-mapping.txt, cross-species-genes.txt
├── level1-technical-validation/   # Qualification_BgeeDB_2_PKSimDB.R
├── level2-human-old-vs-new/       # Qualification_PKSimDB.R
├── level3-cross-species/          # Qualification_CrossSpecies.R
├── level1/ level2/ level3/        # Output directories (02_data/, 03_plots/)
├── DESCRIPTION + renv.lock        # R environment specs
└── README.md                      # Framework documentation
```

### Commits on PR Branch (chronological)
| Commit | Description |
|--------|-------------|
| `8f6ad08` | Expand README and fix helper_Relative_Expression.R |
| `8835ad1` | GitHub Release asset strategy for Mouse archive |
| `08b8991` | Reorganize Code/ folder + improved README |
| `a496a01` | Round 2 CodeRabbit feedback (10 fixes) |
| `dc6d47a` | Optimize qualification component (XLSX→TXT, 500→12 plots) |
| `c46e0b7` | Support optional Old DB + network path in qualification |
| `beae1a7` | Enhance visualization with jitter points + gene labels |
| `fc9056b` | Replace simpleMessage with message in GeneratePKsimDB |
| `7e6e41b` | Remaining CodeRabbit fixes (timeout, dataset typo, dir.create) |
| `9dca5e5` | Stage all qualification outputs and configuration updates |

## What's Left to Build 🔲

### Immediate (Blocking PR Merge)
- **Maintainer review** of PR #7 — all code ready, awaiting human approval

### Post-Merge
1. **Create GitHub Release tag** v3.0.2
2. **Upload Mouse release asset** via `helper_upload_release_asset.sh` (requires GH_TOKEN)
3. **(Optional) Olive baboon** — preliminary work on `feature/add-golden-hamster-and-version-docs` branch; needs full integration testing

## Current Status
- **PR #7**: ✅ All CodeRabbit rounds complete, all tests passing, pushed to remote
- **Branch**: `3-technical-validation-workflow-osp-expression-database-v302` up to date with origin
- **HEAD**: `9dca5e5` — working tree clean

## Known Issues

### Resolved This PR
- ~~`nnemestrina_gene_ensembl` typo for Monkey_PigTailed~~ → Fixed in `GeneratePKsimDB.R`
- ~~`40-YEAR-OLD TAGE` typo in BgeeAges.txt~~ → Fixed to "STAGE"
- ~~Timeout option not restored after function call~~ → Added `on.exit()` snapshot/restore
- ~~`dir.create()` using relative paths~~ → Anchored to `PATH` parameter
- ~~`simpleMessage()`/`simpleError()` misuse~~ → Replaced with `message()`/`stop()`
- ~~XLSX config files not version-controlled~~ → Replaced with plain TXT files
- ~~500+ redundant qualification PNG files~~ → Reduced to ~12 publication-ready plots
- ~~`GuineaPig` vs `Guineapig` case mismatch~~ → Standardized to `Guineapig`
- ~~SQL injection risk in helper_Relative_Expression.R~~ → Parameterized queries
- ~~`NA` type mismatch in helper_SQL_Queries.R~~ → Changed to `NA_real_`
- ~~Regex metacharacters in gene names~~ → `stringr::str_escape()` applied
- ~~Helper files with side effects on source()~~ → Wrapped in named functions

### Ongoing Limitations
- Human full DB (~65 GB) requires 64-bit R + large memory when `COMPUTE_IN_RAM=TRUE`
- Mouse ADME archive (~113 MB) always requires separate release asset upload (not in git)
- BgeeDB API changes between package versions — test after any BgeeDB upgrade

## Evolution of Project Decisions

| Decision | Rationale |
|----------|-----------|
| XLSX → TXT config | Version-controllable, diff-readable, no Excel dependency |
| 500 PNG → 12 PNG | Composite family plots more useful; individual gene plots replaced by cross-species level 3 |
| 3-level qualification | Separates technical (L1) vs biological (L2) vs comparative (L3) validation concerns |
| level1/level2/level3 output dirs | Clean separation, avoids cross-contamination of validation artifacts |
| DESCRIPTION + renv.lock | Ensures reproducibility; CodeRabbit suggested formal R environment specification |
| `.archive/` in .gitignore | Old plots preserved locally for reference without cluttering repo history |
| Mouse as release-asset only | GitHub 100 MB limit; release assets have no per-file size restriction |
