# Project Brief — Gene-Expression-Databases

## Project Identity
- **Repo**: Open-Systems-Pharmacology/Gene-Expression-Databases
- **Current version**: v3.0.2 (in progress, PR #7)
- **Previous release**: v3.0.1 on master

## Core Purpose
Build and distribute PK-Sim–compatible gene expression databases (SQLite `.expressionDB`) from Bgee RNA-Seq data for 18 species. Used by the Open-Systems-Pharmacology (OSP) suite for physiologically-based pharmacokinetic (PBPK) modeling.

## Supported Species (18 total)

| Category | Species |
|---|---|
| **Human** | *Homo sapiens* |
| **PharmaSpecies (9)** | Mouse, Rat, Rabbit, Guinea pig, Dog, Minipig, Monkey (*M. mulatta*), Monkey (*M. fascicularis*), Monkey (Pig-tailed) |
| **AnimalHealthSpecies (8)** | Cattle, Horse, Cat, Chicken, Goat, Sheep, Turkey, Zebrafish |

## Key Goals
1. Automate DB generation from Bgee release 15.2 RNA-Seq TPM data
2. Map gene expression to PK-Sim SQLite schema (32 tables per DB)
3. Provide ADME-only and full DB variants per species
4. Include a three-level qualification/validation framework
5. Distribute compressed archives via GitHub Releases

## Active Branch
`3-technical-validation-workflow-osp-expression-database-v302`  
→ PR #7: "Technical validation workflow for OSP expression database v3.0.2"  
→ Target: master

## Out of Scope
- Creating new VPCs, IAM roles, or AWS infrastructure
- Modifying PK-Sim application source code
- Gene expression data collection (Bgee provides that)
