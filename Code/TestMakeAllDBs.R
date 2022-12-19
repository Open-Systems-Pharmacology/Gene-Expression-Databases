#### Preparations to run the script ####
PATH <- getwd()
source(paste0(PATH,"/Code/1.PrepareBioMarts.R")) # supply function
source(paste0(PATH,"/Code/2.GeneratePKsimDB.R")) # Main code
#source(paste0(PATH,"/Code/DownloadAll_BioMarts.R")) # supply function
options(timeout = 60*60*3) # is needed to allow download of human data (65 GB take longer thatnk 1 min ;-) )
LOCATION = dirname(rstudioapi::getSourceEditorContext()$path)
# depending on your system you might need to set a proxy to enable data download
# Sys.setenv("http_proxy" = "http://PROXY:PORT")
# Sys.setenv("ftp_proxy" = "http://PROXY:PORT")

ALL_SPECIE <- c("Mouse","Rat","Rabbit","Dog","Minipig","Monkey","Monkey_CrabEatingMacaque",
                "GuineaPig","Cat","Chicken","Turkey","Goat","Sheep","Cattle","Horse","Zebrafish")

#### Get gene annotation information ####
# Human information needs to be added first, is the basis for gene homology of other species
PrepareBioMarts("Human")
for (Specie in ALL_SPECIE) {
  PrepareBioMarts(Specie)
}

#### Make databases for pre-clinical species ####
# performance tests
# profvis::profvis(GeneratePKsimDB(SPECIE = "Dog", COMPUTE_IN_RAM = T_No))

# Dog causes error, Variant names can be found from withing PK-Sim  but no date is retrieved; organ maaping might be wrong
RELEASE = "15_0"
GeneratePKsimDB(SPECIE = "Dog"    ,   COMPUTE_IN_RAM = F, ADME_ONLY = F)
GeneratePKsimDB(SPECIE = "Minipig",   COMPUTE_IN_RAM = F, ADME_ONLY = F)
GeneratePKsimDB(SPECIE = "Monkey" ,   COMPUTE_IN_RAM = F, ADME_ONLY = F)
GeneratePKsimDB(SPECIE = "Monkey_CrabEatingMacaque", COMPUTE_IN_RAM = F, ADME_ONLY = F)
GeneratePKsimDB(SPECIE = "Rat"    ,   COMPUTE_IN_RAM = F, ADME_ONLY = F)
GeneratePKsimDB(SPECIE = "Mouse"  ,   COMPUTE_IN_RAM = F, ADME_ONLY = F)
GeneratePKsimDB(SPECIE = "Rabbit" ,   COMPUTE_IN_RAM = F, ADME_ONLY = F)

#### Make databases for animal-health species ####
GeneratePKsimDB(SPECIE = "Chicken",   COMPUTE_IN_RAM = F, ADME_ONLY = T)
GeneratePKsimDB(SPECIE = "Zebrafish", COMPUTE_IN_RAM = F, ADME_ONLY = T)
GeneratePKsimDB(SPECIE = "Horse",     COMPUTE_IN_RAM = F, ADME_ONLY = T)
GeneratePKsimDB(SPECIE = "Cat",       COMPUTE_IN_RAM = F, ADME_ONLY = T)
GeneratePKsimDB(SPECIE = "GuineaPig", COMPUTE_IN_RAM = F, ADME_ONLY = T)
GeneratePKsimDB(SPECIE = "Cattle",    COMPUTE_IN_RAM = F, ADME_ONLY = T)
GeneratePKsimDB(SPECIE = "Turkey",    COMPUTE_IN_RAM = F, ADME_ONLY = T)
GeneratePKsimDB(SPECIE = "Sheep",     COMPUTE_IN_RAM = F, ADME_ONLY = T)
GeneratePKsimDB(SPECIE = "Goat",      COMPUTE_IN_RAM = F, ADME_ONLY = T)

#### Make databases for human  ####
GeneratePKsimDB(SPECIE = "Human", COMPUTE_IN_RAM = F, ADME_ONLY = T)

# #### Make all expression databases ####
# # try human ADME ONLY first, since it is a HUGE data set!
# for (Specie in ALL_SPECIE) {
#   if (Specie == "Human") {
#     GeneratePKsimDB(SPECIE = Specie, USE_PROXY = TRUE, ADME_ONLY = TRUE) 
#   } else {
#     GeneratePKsimDB(SPECIE = Specie, USE_PROXY = TRUE, ADME_ONLY = FALSE) 
#   }
# }


#### Clean up ####
cat("memory in storage: ", memory.size()," available memory: ", memory.size(max = TRUE), "\n")
gc()
