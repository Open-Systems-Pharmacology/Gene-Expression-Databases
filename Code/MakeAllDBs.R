# Preparations to run the script ###
# free workspace
rm(list = ls())

# set working and export directory
setwd(here::here())
PATH <- getwd()

# define release version
RELEASE <- "15_2"

# load dependency functions ####
# support function that builds the gene information look up tables
source(paste0(PATH, "/Code/helper_All_Bgee_organs.R"))
source(paste0(PATH, "/Code/helper_SQL_Commands.R"))
source(paste0(PATH, "/Code/helper_Species.R"))
source(paste0(PATH, "/Code/PrepareBioMarts.R"))
source(paste0(PATH, "/Code/GeneratePKsimDB.R")) # Main code

# set timeout to enable download of large files ####
# is needed to allow download of large data sets
options(timeout = 60 * 60 * 60)

# depending on your system you might need to set a proxy to enable data download
# Sys.setenv("http_proxy" = "http://PROXY:PORT")
# Sys.setenv("ftp_proxy" = "http://PROXY:PORT")

# define species for which to build the gene expression databases ####
source(paste0(PATH, "/Code/helper_Species.R"))

# Get gene annotation information ####
# Human information needs to be added first,
# is the basis for gene homology of other species
PrepareBioMarts(SPECIE = "Human")
for (Specie in PharmaSpecies) {
  PrepareBioMarts(SPECIE = Specie)
}
for (Specie in AnimalHealthSpecies) {
  PrepareBioMarts(SPECIE = Specie)
}

# make PKsimDB for pharmacological species and their ADME genes ####
# run with foreach in parallel
library(foreach)
cl <- parallel::makePSOCKcluster(length(PharmaSpecies))
doParallel::registerDoParallel(cl, cores = length(PharmaSpecies))
foreach::foreach(mol = 1:length(PharmaSpecies)) %dopar% {
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
# doParallel::stopImplicitCluster() #
parallel::stopCluster(cl)
# run with foreach in parallel
library(foreach)
# make PKsimDB for animal health species and their ADME genes ####
cl <- parallel::makePSOCKcluster(length(AnimalHealthSpecies))
doParallel::registerDoParallel(cl, cores = length(AnimalHealthSpecies))
foreach::foreach(mol = 1:length(AnimalHealthSpecies)) %dopar% {
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
# doParallel::stopImplicitCluster()
parallel::stopCluster(cl)

# Make databases for human  ####
GeneratePKsimDB(SPECIE = "Human", COMPUTE_IN_RAM = TRUE, ADME_ONLY = TRUE)
GeneratePKsimDB(SPECIE = "Human", COMPUTE_IN_RAM = TRUE, ADME_ONLY = FALSE)

# Compress ADME Databases for upload
system("bash Code/helper_compress_DBs.sh")

# Run technical validation that data from Bgee is correctly transfered to PK-Sim DB
source("Code/Qualification_BgeeDB_2_PKSimDB.R")

# Run comparison between previous and new expression profiles of ADME genes in OSP-Model-Library
source("Code/Qualification_PKSimDB.R")

# Clean up ####
cat(
  "memory in storage: ", memory.size(),
  " available memory: ", memory.size(max = TRUE), "\n"
)
gc()
