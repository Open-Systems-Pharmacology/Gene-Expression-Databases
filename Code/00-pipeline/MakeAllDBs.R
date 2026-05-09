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
expected_tables <- c(
  paste0(ALL_SPECIE, "_Annotations"),
  paste0(ALL_SPECIE, "_ADME")
)
biomart_cache_complete <- FALSE

if (file.exists(biomart_db)) {
  biomart_conn <- DBI::dbConnect(RSQLite::SQLite(), biomart_db, synchronous = NULL)
  on.exit(try(DBI::dbDisconnect(biomart_conn), silent = TRUE), add = TRUE)
  existing_tables <- DBI::dbListTables(biomart_conn)
  biomart_cache_complete <- all(expected_tables %in% existing_tables)
}

if (biomart_cache_complete) {
  message("BioMart cache already complete for all species; skipping PrepareBioMarts().")
} else {
  PrepareBioMarts(SPECIE = "Human")
  for (Specie in PharmaSpecies) {
    PrepareBioMarts(SPECIE = Specie)
  }
  for (Specie in AnimalHealthSpecies) {
    PrepareBioMarts(SPECIE = Specie)
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
