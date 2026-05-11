# Preparations to run the script ###
# free workspace
rm(list = ls())

# set working and export directory
setwd(here::here())
PATH <- getwd()

# define release version
RELEASE <- "15_2"

# load dependency functions ####
# support function that builds the Bgee organ/age lookup tables (definition only)
source(paste0(PATH, "/Code/03-helpers/helper_All_Bgee_organs.R"))
source(paste0(PATH, "/Code/03-helpers/helper_SQL_Commands.R"))
source(paste0(PATH, "/Code/03-helpers/helper_Species.R"))
source(paste0(PATH, "/Code/01-biomart-prep/PrepareBioMarts.R"))
source(paste0(PATH, "/Code/02-db-generation/GeneratePKsimDB.R")) # Main code

# set timeout to enable download of large files ####
# is needed to allow download of large data sets
options(timeout = 60 * 60 * 2)  # 2 hours — sufficient for large Bgee downloads

# depending on your system you might need to set a proxy to enable data download
# Sys.setenv("http_proxy" = "http://PROXY:PORT")
# Sys.setenv("ftp_proxy" = "http://PROXY:PORT")

# Build Bgee organ/age lookup tables (BgeeOrgans.txt, BgeeAges.txt) ####
# These are cached in BgeeDBs/ and only need rebuilding when the Bgee release changes.
build_bgee_lookup_tables(PATH)

# Get gene annotation information ####
# Human information needs to be added first,
# is the basis for gene homology of other species
biomart_db <- file.path(PATH, "BioMarts", "All_Species_BioMarts.DB")
biomart_cache_metadata_table <- "tab_biomart_cache_metadata"
biomart_cache_key_name <- "pinned_config_md5"

compute_biomart_cache_key <- function(path) {
  config_file <- file.path(path, "Code", "01-biomart-prep", "PrepareBioMarts.R")
  if (!file.exists(config_file)) {
    return("")
  }

  config_lines <- grep(
    "dataset\\s*=|hosts\\s*=\\s*c\\(|expected_version_pattern\\s*=",
    readLines(config_file, warn = FALSE),
    value = TRUE
  )
  if (length(config_lines) == 0) {
    return("")
  }

  tmp <- tempfile("biomart-cache-key-")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(config_lines, tmp, useBytes = TRUE)
  as.character(unname(tools::md5sum(tmp)))
}

read_biomart_cache_key <- function(conn, table_name, key_name) {
  if (!DBI::dbExistsTable(conn, table_name)) {
    return("")
  }

  result <- tryCatch(
    DBI::dbGetQuery(
      conn,
      paste0("SELECT cache_key FROM ", table_name, " WHERE key_name = ? LIMIT 1"),
      params = list(key_name)
    ),
    error = function(e) data.frame(cache_key = character(0))
  )

  if (nrow(result) == 0) {
    return("")
  }
  as.character(result$cache_key[1])
}

write_biomart_cache_key <- function(conn, table_name, key_name, cache_key) {
  DBI::dbExecute(
    conn,
    paste0(
      "CREATE TABLE IF NOT EXISTS ", table_name,
      " (key_name TEXT PRIMARY KEY, cache_key TEXT NOT NULL, updated_at TEXT NOT NULL)"
    )
  )

  DBI::dbExecute(
    conn,
    paste0(
      "INSERT INTO ", table_name, " (key_name, cache_key, updated_at) VALUES (?, ?, ?) ",
      "ON CONFLICT(key_name) DO UPDATE SET cache_key = excluded.cache_key, updated_at = excluded.updated_at"
    ),
    params = list(key_name, cache_key, as.character(Sys.time()))
  )
}

expected_tables <- c(
  paste0(ALL_SPECIE, "_Annotations"),
  paste0(ALL_SPECIE, "_ADME")
)
biomart_cache_complete <- FALSE
existing_tables <- character(0)
current_biomart_cache_key <- compute_biomart_cache_key(PATH)
stored_biomart_cache_key <- ""
has_expected_tables <- FALSE

if (file.exists(biomart_db)) {
  biomart_conn <- DBI::dbConnect(RSQLite::SQLite(), biomart_db, synchronous = NULL)
  tryCatch({
    existing_tables <- DBI::dbListTables(biomart_conn)
    has_expected_tables <- all(expected_tables %in% existing_tables)
    if (has_expected_tables) {
      stored_biomart_cache_key <- read_biomart_cache_key(
        biomart_conn,
        biomart_cache_metadata_table,
        biomart_cache_key_name
      )
    }
  }, finally = try(DBI::dbDisconnect(biomart_conn), silent = TRUE))
}

biomart_cache_complete <- has_expected_tables &&
  nzchar(current_biomart_cache_key) &&
  identical(stored_biomart_cache_key, current_biomart_cache_key)

if (biomart_cache_complete) {
  message("BioMart cache already complete for all species; skipping PrepareBioMarts().")
} else {
  target_species <- c("Human", PharmaSpecies, AnimalHealthSpecies)
  missing_species <- target_species[vapply(target_species, function(species) {
    annotation_table <- paste0(species, "_Annotations")
    adme_table <- paste0(species, "_ADME")
    !(annotation_table %in% existing_tables && adme_table %in% existing_tables)
  }, logical(1))]

  species_to_refresh <- if (has_expected_tables) {
    target_species
  } else {
    missing_species
  }

  for (Specie in species_to_refresh) {
    PrepareBioMarts(SPECIE = Specie)
  }

  if (nzchar(current_biomart_cache_key) && file.exists(biomart_db)) {
    biomart_conn <- DBI::dbConnect(RSQLite::SQLite(), biomart_db, synchronous = NULL)
    tryCatch({
      write_biomart_cache_key(
        biomart_conn,
        biomart_cache_metadata_table,
        biomart_cache_key_name,
        current_biomart_cache_key
      )
    }, finally = try(DBI::dbDisconnect(biomart_conn), silent = TRUE))
  }
}

# make PKsimDB for pharmacological species and their ADME genes ####
# run with foreach in parallel
library(foreach)
cl <- parallel::makePSOCKcluster(length(PharmaSpecies))
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
doParallel::registerDoParallel(cl, cores = length(PharmaSpecies))
foreach::foreach(
  mol = seq_along(PharmaSpecies),
  .export = c("PharmaSpecies", "RELEASE", "PATH",
              "GeneratePKsimDB", "ALL_SPECIE",
              "CREATE_TABLE", "VIEW_TABLE", "INDIZES"),
  .packages = c("DBI", "RSQLite", "dplyr", "BgeeDB",
                "biomaRt", "tibble", "readr",
                "tidyr", "tidyselect", "stringr",
                "rlang", "here")
) %dopar% {
  Specie <- PharmaSpecies[mol]
  # for (Specie in PharmaSpecies) {
  GeneratePKsimDB(
    SPECIE = Specie,
    COMPUTE_IN_RAM = TRUE,
    ADME_ONLY = TRUE,
    RELEASE = RELEASE
  )
  GeneratePKsimDB(
    SPECIE = Specie,
    COMPUTE_IN_RAM = TRUE,
    ADME_ONLY = FALSE,
    RELEASE = RELEASE
  )
}
# run with foreach in parallel
library(foreach)
# make PKsimDB for animal health species and their ADME genes ####
cl <- parallel::makePSOCKcluster(length(AnimalHealthSpecies))
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
doParallel::registerDoParallel(cl, cores = length(AnimalHealthSpecies))
foreach::foreach(
  mol = seq_along(AnimalHealthSpecies),
  .export = c("AnimalHealthSpecies", "RELEASE", "PATH",
              "GeneratePKsimDB", "ALL_SPECIE",
              "CREATE_TABLE", "VIEW_TABLE", "INDIZES"),
  .packages = c("DBI", "RSQLite", "dplyr", "BgeeDB",
                "biomaRt", "tibble", "readr",
                "tidyr", "tidyselect", "stringr",
                "rlang", "here")
) %dopar% {
  Specie <- AnimalHealthSpecies[mol]
  # for (Specie in AnimalHealthSpecies) {
  GeneratePKsimDB(
    SPECIE = Specie,
    COMPUTE_IN_RAM = TRUE,
    ADME_ONLY = TRUE,
    RELEASE = RELEASE
  )
  GeneratePKsimDB(
    SPECIE = Specie,
    COMPUTE_IN_RAM = TRUE,
    ADME_ONLY = FALSE,
    RELEASE = RELEASE
  )
}

# Make databases for human  ####
GeneratePKsimDB(SPECIE = "Human", COMPUTE_IN_RAM = TRUE, ADME_ONLY = TRUE,  RELEASE = RELEASE)
GeneratePKsimDB(SPECIE = "Human", COMPUTE_IN_RAM = TRUE, ADME_ONLY = FALSE, RELEASE = RELEASE)

# Compress ADME Databases for upload
system("bash Code/05-utilities/helper_compress_DBs.sh")

# Run technical validation that data from Bgee is correctly transfered to PK-Sim DB
source(paste0(PATH, "/Code/04-qualification/level1-technical-validation/Qualification_BgeeDB_2_PKSimDB.R"))

# Run comparison between previous and new expression profiles of ADME genes in OSP-Model-Library
source(paste0(PATH, "/Code/04-qualification/level2-human-old-vs-new/Qualification_PKSimDB.R"))

# Run cross-species comparison for selected ADME genes in key preclinical species
source(paste0(PATH, "/Code/04-qualification/level3-cross-species/Qualification_CrossSpecies.R"))

# Clean up ####
mem_info <- gc(reset = TRUE)
cat("memory used: ", round(sum(mem_info[, 2]), 1), " Mb\n")
