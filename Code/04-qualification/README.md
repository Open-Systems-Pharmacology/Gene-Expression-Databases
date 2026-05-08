# Qualification & Validation Framework

This folder contains the comprehensive qualification framework for OSP Expression Databases across three levels of validation.

## Folder Structure

```
Code/04-qualification/
├── 01_config/                    # Shared configuration (version-controlled)
│   ├── proteins-validation.txt        # 57 genes for validation (one per line)
│   ├── container-mapping.txt          # Tissue-to-container mappings (tab-separated)
│   └── cross-species-genes.txt        # 40 ADME genes for cross-species validation
│
├── level1/                       # LEVEL 1: Technical Validation
│   ├── 02_data/                  # Intermediate data & caches (generated)
│   │   ├── DATA_HUMAN_GSE30611_ERX011211_BgeeDB.csv      (Bgee API cache)
│   │   └── technical_validation_data.csv                 (Bgee → OSP DB comparison)
│   └── 03_plots/                 # Publication-ready plots (generated)
│       └── technical_validation_xy_*.png                 (7 scatter plots by ADME family)
│
├── level2/                       # LEVEL 2: Biological Validation (Human)
│   ├── 02_data/                  # Intermediate data & caches (generated)
│   │   ├── DATA_HUMAN_GSE30611_ERX011211_OSP_DB.csv      (OSP DB extract cache)
│   │   └── validation_summary.csv                        (37 genes × 12 metrics)
│   └── 03_plots/                 # Publication-ready plots (generated)
│       └── human_old_vs_new_*.png                        (4 composite violin plots by family)
│
├── level3/                       # LEVEL 3: Cross-Species Validation
│   ├── 02_data/                  # Intermediate data (generated)
│   │   └── cross_species_data_*.csv                      (Per-gene cross-species data)
│   └── 03_plots/                 # Cross-species plots (generated)
│       └── gene_*.png                                    (40 individual gene plots across species)
│
├── DESCRIPTION                   # R package requirements specification
├── renv.lock                     # R environment lock file (all package versions)
│
├── Qualification_BgeeDB_2_PKSimDB.R    # Wrapper: runs Level 1 validation
├── Qualification_PKSimDB.R             # Wrapper: runs Level 2 validation
├── Qualification_CrossSpecies.R        # Wrapper: runs Level 3 validation
│
└── level{1,2,3}-{name}/          # Actual validation scripts (sourced by wrappers)
    ├── Qualification_BgeeDB_2_PKSimDB.R       (Level 1)
    ├── Qualification_PKSimDB.R                (Level 2)
    └── Qualification_CrossSpecies.R          (Level 3)
```

## Validation Levels

### Level 1: Technical Validation (`level1/`)
**Purpose**: Validate Bgee source data integrity and integration into PKSimDB

- **Input**: Bgee GSE30611 ERX011211 sample (human liver, cached)
- **Output**: 7 X-Y scatter plots + validation CSV
- **Families Tested**: CYP, ABC, UGT, SULT, SLC, CES, Other
- **Method**: Bgee TPM (X-axis) vs OSP sample count (Y-axis)
- **Features**: Gene symbol labels, linear regression overlay, family-specific colors

**Run Command**:
```r
source("Code/04-qualification/Qualification_BgeeDB_2_PKSimDB.R")
# or directly:
source("Code/04-qualification/level1-technical-validation/Qualification_BgeeDB_2_PKSimDB.R")
```

---

### Level 2: Biological Validation (`level2/`)
**Purpose**: Compare ADME gene expression profiles between old (RT-PCR) and new (Bgee-based) databases

- **Input**: Old + New human expressionDB files
- **Output**: 4 composite violin plots + validation summary CSV
- **Families Tested**: CYP, ABC, UGT, SLC
- **Method**: Relative expression violin plots with overlaid jittered points
- **Features**: OSP-branded theme, log10 scale, correlation statistics

**Run Command**:
```r
source("Code/04-qualification/Qualification_PKSimDB.R")
# or directly:
source("Code/04-qualification/level2-human-old-vs-new/Qualification_PKSimDB.R")
```

---

### Level 3: Cross-Species Validation (`level3/`)
**Purpose**: Compare gene expression patterns across 5 species for a 40-gene ADME panel

- **Input**: Expression databases for Human, Macaque, Dog, Rat, Mouse
- **Output**: 40 individual gene plots showing tissue distributions per species
- **Genes**: ADME panel (CYP, ABC, UGT, SLC, SULT, CES genes)
- **Method**: Tissue-level relative expression (normalized per species)
- **Features**: Species-specific colors, tissue distributions, cross-species comparison

**Run Command**:
```r
source("Code/04-qualification/Qualification_CrossSpecies.R")
# or directly:
source("Code/04-qualification/level3-cross-species/Qualification_CrossSpecies.R")
```

---

## R Environment Management

### DESCRIPTION File
- Package name, version, and dependency declarations
- R version requirement: **4.4.1 or higher**
- All dependencies listed with tested versions
- Both CRAN and Bioconductor packages documented

### renv.lock File
- Complete dependency graph with versions and sources
- Lock file for reproducible R environments
- Useful for containerization and CI/CD workflows

### Installation (One-Time Setup)

#### Option 1: Manual Installation
```r
# Install from CRAN
install.packages(c(
  "here", "dplyr", "readr", "ggplot2", "ggrepel", 
  "scales", "tidyr", "DBI", "RSQLite", "stringr"
))

# Install from Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("BgeeDB", "biomaRt"))

# Other packages
install.packages(c("tidyselect", "rlang", "foreach", "doParallel"))
```

#### Option 2: Using renv (Recommended for Reproducibility)
```r
# Install renv if not already present
install.packages("renv")

# From project root:
renv::restore("Code/04-qualification/renv.lock")
```

---

## Quick Start

### Run All Three Validation Levels
```bash
cd /path/to/Gene-Expression-Databases

# Level 1: Technical validation (Bgee → OSP integration)
Rscript Code/04-qualification/Qualification_BgeeDB_2_PKSimDB.R

# Level 2: Biological validation (Old vs New DB)
Rscript Code/04-qualification/Qualification_PKSimDB.R

# Level 3: Cross-species validation
Rscript Code/04-qualification/Qualification_CrossSpecies.R
```

### Run from Within R
```r
setwd("/path/to/Gene-Expression-Databases")

# Level 1
source("Code/04-qualification/Qualification_BgeeDB_2_PKSimDB.R")

# Level 2
source("Code/04-qualification/Qualification_PKSimDB.R")

# Level 3
source("Code/04-qualification/Qualification_CrossSpecies.R")
```

---

## Key Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| **here** | 1.0.1 | Path resolution |
| **dplyr** | 1.2.1 | Data manipulation |
| **ggplot2** | 4.0.0 | Visualization |
| **ggrepel** | 0.9.6 | Text label placement |
| **scales** | 1.4.0 | Color scales |
| **tidyr** | 1.3.1 | Data reshaping |
| **DBI** | 1.2.3 | Database interface |
| **RSQLite** | 2.4.3 | SQLite connector |
| **BgeeDB** | 2.32.0 | Bgee data access |
| **biomaRt** | 2.62.0 | BioMart annotation |
| **stringr** | 1.5.1 | String utilities |
| **foreach** | 1.5.2 | Parallel iteration |
| **doParallel** | 1.0.17 | Parallel backend |

---

## Troubleshooting

### "Cannot find package X"
Ensure all packages are installed (see Installation section above).

### Level 1 fails: "BgeeDB API signature mismatch"
This occurs when BgeeDB package version differs from tested 2.32.0. Update:
```r
BiocManager::install("BgeeDB", version = "latest")
```

### Level 2 fails: "Old DB not found"
The script handles missing Old DB gracefully (Scenario 2 fallback). If intentional, ignore warning.  
To provide Old DB: Set network path in script or place file locally.

### Level 3 runs slowly
Cross-species validation is I/O intensive (5 databases). Expected runtime: 2–5 minutes.

---

## Output Quality

- **Data Exports**: CSV files for downstream analysis/publication supplementary material
- **Plots**: 300 DPI PNG, publication-ready
- **Total Output Size**: ~4 MB for all three levels

---

## Version History

- **v3.0.2** (June 2024): Level-specific output reorganization, simplified folder structure
- **Bgee Release**: 15.2 (2024-05-21)
- **R Version**: 4.4.1 (2024-06-14)

---

## References

- [Bgee Database](http://bgee.org/)
- [BgeeDB R Package](https://bioconductor.org/packages/release/bioc/html/BgeeDB.html)
- [Open Systems Pharmacology](https://github.com/Open-Systems-Pharmacology)
