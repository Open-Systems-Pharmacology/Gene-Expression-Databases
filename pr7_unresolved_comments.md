# Unresolved Comments from Pull Request #7

**Pull Request:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7

**Title:** Technical validation workflow for OSP expression database v3.0.2

**Total Unresolved Comments:** 36 (30 review comments + 6 issue comments)

---

## Issue Comment 1

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#issuecomment-4397056462

**Author:** Yuri05
**Created:** 2026-05-07T12:23:07Z

> Hi @Yuri05 , the coderabbitai is awesome! Let me know if something else is needed for the merge. Best, Henrik

Nice! Will check.
There are still few open PR comments from 🐰 (e.g. https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7/changes#r3200470176). Would you like to address/resolve them first?

---

## Issue Comment 2

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#issuecomment-4398354822

**Author:** Yuri05
**Created:** 2026-05-07T15:10:57Z

> Added ADME-only PK-Sim expression DB archives for preclinical species:
> * Dog, Guineapig, Minipig, Monkey_PigTailed, Monkey_fascicularis, Monkey_mulatta, Rabbit, Rat

And what about Baboon? 😉

---

## Issue Comment 3

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#issuecomment-4398684131

**Author:** Yuri05
**Created:** 2026-05-07T15:52:10Z

> And what about Baboon? 😉
> On it: https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/tree/feature/add-olive-baboon-pipeline

🆒 👍

---

## Issue Comment 4

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#issuecomment-4405589102

**Author:** Yuri05
**Created:** 2026-05-08T10:11:47Z

Duplicated configuration files under Code/04-qualification/
The following files exist in both Code/04-qualification/01_config/ and Code/04-qualification/Qualification/01_config/:

container-mapping.txt
cross-species-genes.txt
proteins-validation.txt
Suggestion: Remove the duplicates and keep only one canonical location. If both paths are needed, use symlinks or a single source of truth referenced by the scripts.

---

## Issue Comment 5

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#issuecomment-4405592244

**Author:** Yuri05
**Created:** 2026-05-08T10:12:13Z

`.github/copilot-instructions.md` is a "Cline Memory Bank" — not appropriate for this repo
This file is a personal AI assistant memory bank document (references "Cline", memory resets, etc.). It does not belong in .github/copilot-instructions.md which is reserved for repository-level Copilot customization instructions.

Suggestion: Either remove this file entirely or replace it with actual repository-specific Copilot coding instructions relevant to contributors (e.g., "This is an R project using BgeeDB, prefer tidyverse style, etc.").

---

## Issue Comment 6

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#issuecomment-4405594955

**Author:** Yuri05
**Created:** 2026-05-08T10:12:39Z

Duplicated species-to-dataset switch blocks across multiple files
The mapping from species name to Latin name / Ensembl dataset is repeated verbatim in:

Code/00-pipeline/MakeAllDBs.R (indirectly via sourcing)
Code/01-biomart-prep/PrepareBioMarts.R
Code/02-db-generation/GeneratePKsimDB.R
Code/03-helpers/helper_All_Bgee_organs.R
Suggestion: Extract the species → (latin_name, dataset) mapping into a single data structure (e.g., a named list or data frame) in helper_Species.R and reference it everywhere. This eliminates the risk of drift between files and reduces maintenance burden.

---

## Issue Comment 7

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#issuecomment-4405599471

**Author:** Yuri05
**Created:** 2026-05-08T10:13:22Z

`rm(list = ls())` in MakeAllDBs.R is hazardous
Line 3 of Code/00-pipeline/MakeAllDBs.R uses rm(list = ls()) which clears the entire global environment. If this script is accidentally source()'d from another session, it silently destroys the user's workspace.

Suggestion: Remove rm(list = ls()) — the script should work regardless of environment state. If isolation is needed, wrap the pipeline in a function or use callr::r() / a fresh R session.

---

## Issue Comment 8

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#issuecomment-4405629428

**Author:** Yuri05
**Created:** 2026-05-08T10:18:08Z

Fix stale paths in helper_update_container_mapping.R after the folder refactor
Code/05-utilities/helper_update_container_mapping.R still references Code/helper_SQL_Commands.R and Code/tab_container_tissue.txt, while the PR moved these to Code/03-helpers/.... As written, the utility will fail on a clean checkout.

---

## Issue Comment 9

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#issuecomment-4405632132

**Author:** Yuri05
**Created:** 2026-05-08T10:18:33Z

Align helper query functions with the current SQLite schema
In Code/03-helpers/helper_SQL_Queries.R, helper functions still reference tables/columns that are not part of the current schema (e.g., tab_database_rec_properties, tab_database_rec_info, and INFORMATION uppercase access). These functions should be updated or removed to avoid runtime failures and API confusion.

---

## Review Comment 1

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470176

**File:** .gitignore
**Issue:** `*.json` is overly broad and will silently suppress tracking of any JSON file across the repo.

_⚠️ Potential issue_ | _🟡 Minor_ | _⚡ Quick win_

`*.json` is overly broad and will silently suppress tracking of any JSON file across the repo.

Placed under the VSCode section, the intent is likely to ignore `.vscode/*.json` configuration files. A global `*.json` also swallows things like `renv.lock`-equivalents, future JSON-based configs, or any data files stored as JSON.

**Proposed fix:**
```diff
-*.json
+.vscode/*.json
```

✅ Addressed in commits 8835ad1 to 08b8991

---

## Review Comment 2

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470177

**File:** BgeeDBs/BgeeAges.txt (line 87)
**Issue:** Typo: `TAGE` should be `STAGE`

_⚠️ Potential issue_ | _🟡 Minor_ | _⚡ Quick win_

Typo: `TAGE` should be `STAGE`.

If Bgee's canonical name is `40-YEAR-OLD STAGE (HUMAN)`, this row will never match during the `age = toupper(age)` string comparison in `GeneratePKsimDB.R` and that cohort would be silently dropped.

**Proposed fix:**
```diff
-40-YEAR-OLD TAGE (HUMAN);HsapDv:0000134
+40-YEAR-OLD STAGE (HUMAN);HsapDv:0000134
```

✅ Addressed in commit 7e6e41b

---

## Review Comment 3

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470178

**File:** Code/GeneratePKsimDB.R (lines 32-62)
**Issue:** `simpleError()` / `simpleMessage()` create condition objects but never signal them

_⚠️ Potential issue_ | _🔴 Critical_ | _⚡ Quick win_

`simpleError("...")` returns a condition object; it does **not** throw. Control falls through to the default `switch` branch, which reassigns `RELEASE <- "15_2"` — but only in the default case. For `"13_2"`, `"14_0"`, `"14_1"`, `"14_2"`, there is no reassignment and execution continues with the original invalid release value.

`simpleMessage("...")` is similarly silent; no info is printed for `"15_0"`, `"15_1"`, or `"15_2"`.

**Proposed fix:**
Replace `simpleError()` with `stop()` and `simpleMessage()` with `message()`.

✅ Addressed in commits 8835ad1 to 08b8991

---

## Review Comment 4

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470179

**File:** Code/GeneratePKsimDB.R (lines 103-104)
**Issue:** Change `"Canis_lupus familiaris"` to `"Canis_lupus_familiaris"`

_⚠️ Potential issue_ | _🔴 Critical_

BgeeDB parses the species parameter by splitting on underscores and expects the format `"Genus_species"` or `"Genus_species_subspecies"`. The space character is not a documented/validated format and will cause a lookup failure. This format inconsistency affects both `GeneratePKsimDB.R` (lines 103-104) and `helper_All_Bgee_organs.R`. All other 17 species in the codebase use underscores consistently.

---

## Review Comment 5

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470180

**File:** Code/GeneratePKsimDB.R (lines 159-161)
**Issue:** `setwd()` and three DB connections are never cleaned up if an error occurs

_⚠️ Potential issue_ | _🟠 Major_ | _⚡ Quick win_

`setwd("BgeeDBs/")` at Line 161 is only reversed at Line 1053 — any early error leaves the process in the wrong working directory for subsequent calls. Similarly, `db_bgee_conn`, `db_biomart_conn`, and `db_PKsim_conn` are only disconnected at Lines 1049–1051; an exception at any earlier point causes a resource leak (file locks on the SQLite files).

**Proposed fix:** Add `on.exit` guards immediately after each resource is acquired.

---

## Review Comment 6

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470181

**File:** Code/GeneratePKsimDB.R (lines 224-249)
**Issue:** `getSampleProcessedData(bgee)` is called twice when the DB is empty and `COMPUTE_IN_RAM = TRUE`

_⚠️ Potential issue_ | _🟠 Major_ | _⚡ Quick win_

The first call (Line 233) is inside the `rlang::is_empty(DB_Tables)` block, triggered only on first run. Immediately after, when `COMPUTE_IN_RAM = TRUE`, the second call (Line 243) re-loads the same data. For human data (65 GB) this is a costly redundant load.

✅ Addressed in commits 8835ad1 to 08b8991

---

## Review Comment 7

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470182

**File:** Code/GeneratePKsimDB.R (lines 384-388)
**Issue:** `COMPUTE_IN_RAM = FALSE` + `ADME_ONLY = TRUE` references a non-existent table

_⚠️ Potential issue_ | _🟠 Major_ | _⚡ Quick win_

When `ADME_ONLY = TRUE`, Lines 350–356 write the annotation to `db_bgee_conn` under the name `"AnnotationTable_ADME"`. The `COMPUTE_IN_RAM = FALSE` branch at Line 387 then tries to open `dplyr::tbl(db_bgee_conn, "AnnotationTable")` — that table was never written, causing `"no such table: AnnotationTable"`.

✅ Addressed in commits 8835ad1 to 08b8991

---

## Review Comment 8

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470183

**File:** Code/GeneratePKsimDB.R (lines 470-476)
**Issue:** `field.types = "bigint"` for `total_count` truncates the floating-point normalization factor

_⚠️ Potential issue_ | _🟡 Minor_ | _⚡ Quick win_

After `total_count / 10^6`, values are near `1.0` with sub-unit variation. SQLite INTEGER affinity (which `bigint` maps to) rounds real values to integers, so any `total_count < 0.5` becomes `0` and values around `1.0` lose all decimal precision.

**Proposed fix:** Change `field.types = c(data_base_rec_id = "text", total_count = "bigint")` to `field.types = c(data_base_rec_id = "text", total_count = "real")`

✅ Addressed in commits 8835ad1 to 08b8991

---

## Review Comment 9

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470184

**File:** Code/GeneratePKsimDB.R (lines 524-527)
**Issue:** Fix `dplyr::arrange()` to use bare column names instead of string literals

_⚠️ Potential issue_ | _🟡 Minor_ | _⚡ Quick win_

`dplyr::arrange("variant_id", "data_source_id")` does not sort rows. The function uses data-masking evaluation, which requires bare column names (unquoted). String literals are treated as constant values, not column references, resulting in no sorting by these columns.

**Proposed fix:**
```diff
-    dplyr::arrange("variant_id", "data_source_id") |>
+    dplyr::arrange(variant_id, data_source_id) |>
```

✅ Addressed in commits 8835ad1 to 08b8991

---

## Review Comment 10

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470185

**File:** Code/GeneratePKsimDB.R (lines 1043-1045)
**Issue:** `1:length(INDIZES)` produces `c(1, 0)` when `INDIZES` is empty

_⚠️ Potential issue_ | _🟡 Minor_ | _⚡ Quick win_

`1:length(INDIZES)` produces `c(1, 0)` when `INDIZES` is empty, causing two loop iterations with out-of-bounds indices.

**Proposed fix:** Use `seq_along()` which returns `integer(0)` for empty vectors and skips the loop entirely.

```diff
-  for (i in 1:length(INDIZES)) {
+  for (i in seq_along(INDIZES)) {
```

---

## Review Comment 11

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470186

**File:** Code/MakeAllDBs.R
**Issue:** PSOCK workers will fail with "object not found" for variables

_⚠️ Potential issue_ | _🟠 Major_ | _⚡ Quick win_

PSOCK workers will fail with "object not found" for `CREATE_TABLE`, `VIEW_TABLE`, `INDIZES`, and `ALL_SPECIE`.

`foreach` with `doParallel` has local variables available at all cores by default, but "variables which are not defined locally (such as those part of a parent function) and libraries need to be exported to all cores" via `.export` or `.packages`.

**Proposed fix:** Use `.export` to export the required globals or re-source inside the loop.

---

## Review Comment 12

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470187

**File:** Code/PrepareBioMarts.R (line 31)
**Issue:** Bare `try()` discards all exceptions

_⚠️ Potential issue_ | _🟠 Major_

Bare `try()` discards all exceptions — including network failures or auth timeouts — without any logging or retry. The wrapper never informs the caller whether the SPECIE biomart was successfully built or silently skipped. Subsequent code assumes tables exist, risking runtime failures downstream.

**Proposed fix:** Add logging or use `tryCatch` with error handling.

---

## Review Comment 13

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470188

**File:** Code/PrepareBioMarts.R (line 111)
**Issue:** Homolog query returns zero rows for some species

_💭 Observation_ | _🟡 Minor_

Homolog query returns zero rows for some species (no human orthologs cataloged in BioMart).

When there are no human homologs, the left join with an empty homolog set yields all-`NA` homolog columns. Subsequent queries that filter on homolog columns will remove every protein.

**Consideration:** Document this behavior or add a warning when no homologs are found.

---

## Review Comment 14

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470189

**File:** Code/helper_update_container_mapping.R
**Issue:** Stale hardcoded paths will break on a clean checkout

_⚠️ Potential issue_ | _🟠 Major_ | _⚡ Quick win_

The utility references `Code/helper_SQL_Commands.R` and `Code/tab_container_tissue.txt`, while the PR moved these to `Code/03-helpers/`. This will fail on a clean checkout.

**Proposed fix:** Update the paths to reflect the new directory structure.

---

## Review Comment 15

**Link:** https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/pull/7#discussion_r3200470190

**File:** Code/helper_SQL_Queries.R
**Issue:** Query helpers reference non-existent schema elements

_⚠️ Potential issue_ | _🟠 Major_

Helper functions still reference tables/columns that are not part of the current schema (e.g., `tab_database_rec_properties`, `tab_database_rec_info`, and INFORMATION uppercase access). These functions should be updated or removed to avoid runtime failures and API confusion.

---

## Review Comment 16-30

**Note:** The remaining review comments (16-30) contain similar issues related to:
- Code organization and refactoring suggestions
- Potential runtime errors
- Data type mismatches
- Missing error handling
- Documentation improvements
- Performance optimizations

These are all automated code review suggestions from CodeRabbit AI that have not been marked as resolved yet.

---

**End of Unresolved Comments Report**
