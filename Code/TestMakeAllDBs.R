# Preparations to run the script ###
# set working and export directory
# LOCATION <- dirname(rstudioapi::getSourceEditorContext()$path)
PATH <- getwd()

# initiate R environment in which this code was developed
# renv::restore()

# load dependency functions ####
# support function that builds the gene information look up tables
source(paste0(PATH, "/Code/PrepareBioMarts.R"))
source(paste0(PATH, "/Code/GeneratePKsimDB.R")) # Main code

# set timeout to enable download of large files ####
# is needed to allow download of large data sets
options(timeout = 60 * 60 * 15)

# depending on your system you might need to set a proxy to enable data download
# Sys.setenv("http_proxy" = "http://PROXY:PORT")
# Sys.setenv("ftp_proxy" = "http://PROXY:PORT")

# define species for which to build the gene expression databases ####
PharmaSpecies <- c(
  "Mouse", "Rat", "Rabbit", "GuineaPig", "Dog", "Minipig",
  "Monkey_mulatta", "Monkey_fascicularis", "Monkey_PigTailed"
)

AnimalHealthSpecies <- c("Cattle", "Horse", "Cat", "Chicken",
                         "Goat", "Sheep", "Turkey", "Zebrafish")

# Get gene annotation information ####
# Human information needs to be added first,
# is the basis for gene homology of other species
PrepareBioMarts("Human")
for (Specie in PharmaSpecies) {
  PrepareBioMarts(Specie)
}

# Make databases for pre-clinical species all data ####
# performance tests
# profvis::profvis(GeneratePKsimDB(SPECIE = "Dog", COMPUTE_IN_RAM = FALSE))

# define release version
RELEASE <- "15_2"

# make PKsimDB for pharmacological species and their ADME genes ####
for (Specie in PharmaSpecies) {
  GeneratePKsimDB(
    SPECIE = Specie,
    COMPUTE_IN_RAM = TRUE,
    ADME_ONLY = TRUE,
    RELEASE = RELEASE
  )
}

# make PKsimDB for animal health species and their ADME genes ####
for (Specie in AnimalHealthSpecies) {
  GeneratePKsimDB(
    SPECIE = Specie,
    COMPUTE_IN_RAM = TRUE,
    ADME_ONLY = TRUE,
    RELEASE = RELEASE
  )
}

# Make databases for human  ####
GeneratePKsimDB(SPECIE = "Human", COMPUTE_IN_RAM = FALSE, ADME_ONLY = TRUE)

# Clean up ####
cat("memory in storage: ", memory.size(), 
    " available memory: ", memory.size(max = TRUE), "\n")
gc()
