# Product Context — Gene-Expression-Databases

## Why This Project Exists
PK-Sim (part of the OSP suite) requires species-specific gene expression data to parameterize PBPK models — specifically, relative expression levels of ADME-relevant genes (CYPs, UGTs, SLCs, ABCs, SULTs, CES) across tissues and age groups.

Previously this data was sourced from RT-PCR measurements. This project replaces and supplements that with bulk RNA-Seq TPM data from Bgee (a curated gene expression database maintained by the Swiss Institute of Bioinformatics).

## Problems It Solves
1. **Scalability**: Manual RT-PCR data collection is slow; Bgee covers thousands of samples automatically
2. **Coverage**: Bgee provides expression data for 18 species vs. limited RT-PCR coverage
3. **Reproducibility**: Scripted pipeline ensures identical results on any machine, for any Bgee release
4. **Traceability**: Git-tracked code + version-pinned R packages (renv.lock) = full audit trail
5. **ADME Focus**: Filters to ADME gene families relevant to pharmacokinetic modeling

## How It Should Work

### User-Facing Workflow
1. Researcher downloads pre-built `.expressionDB.tar.gz` archives from GitHub Releases
2. Archives are imported into PK-Sim as gene expression databases
3. PK-Sim uses these to parameterize tissue-specific enzyme/transporter expression in PBPK models

### Developer Workflow
1. Run `Code/01-biomart-prep/PrepareBioMarts.R` to refresh gene annotations from Ensembl
2. Run `Code/00-pipeline/MakeAllDBs.R` to regenerate all species databases
3. Run qualification scripts to validate data integrity
4. Compress and upload archives via `Code/05-utilities/`
5. Create GitHub Release tag (e.g., v3.0.2)

## User Experience Goals
- **For modelers**: Drop-in `.expressionDB` files, compatible with PK-Sim immediately
- **For developers**: Clean, well-documented R pipeline with clear stage separation
- **For reviewers**: Qualification plots show Old (RT-PCR) vs New (Bgee) comparison + cross-species ADME profiles
- **For auditors**: renv.lock + DESCRIPTION pin exact R package versions; git history shows all changes

## Data Sources
- **Bgee release 15.2** (2024-05-21): RNA-Seq TPM expression data for all species
- **Ensembl BioMart**: Gene annotations (IDs, symbols, synonyms, ADME classification, human orthologs)
