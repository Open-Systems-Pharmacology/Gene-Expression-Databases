#### Get information on annotated organs from all species
## This information is needed to map the organ information from data source to the 
## PK-Sim internal organ compartment nomenclature (found in TAB_CONTAINER_TISSUE.txt)
## If names are not matched the expression data can not be linked to the PK-Sim organ compartment!
require("BgeeDB", quietly = T)
require("dplyr", quietly = T)
ALL_SPECIE <- c("Dog","Mouse","Rat","Rabbit","Monkey","Minipig","Zebrafish","Cattle","Horse","Cat","GuineaPig","Chicken","Human")

Organs <- tibble()
Ages <- Organs

for (SPECIE in ALL_SPECIE) {
  switch(SPECIE,
         #### Clinical specie ####
         Human = {SPECIE_LAT <- "Homo_sapiens"
         DATASET = "hsapiens_gene_ensembl"},
         #### Pre-clinical species ####
         Monkey = {SPECIE_LAT <- "Macaca_mulatta"
         DATASET = "mmulatta_gene_ensembl"},
         Minipig = {SPECIE_LAT <- "Sus_scrofa"
         DATASET = "sscrofa_gene_ensembl"}, 
         Dog = {SPECIE_LAT <- "Canis_lupus familiaris"
         DATASET = "cfamiliaris_gene_ensembl" },
         Mouse = { SPECIE_LAT <- "Mus_musculus"
         DATASET = "mmusculus_gene_ensembl"},
         Rat = {SPECIE_LAT <- "Rattus_norvegicus"
         DATASET = "rnorvegicus_gene_ensembl"},
         Rabbit = {SPECIE_LAT <- "Oryctolagus_cuniculus"
         DATASET = "ocuniculus_gene_ensembl"},
         #### Animal health relevant species ####
         Zebrafish = {SPECIE_LAT <- "Danio_rerio" 
         DATASET = "drerio_gene_ensembl"},
         Cattle = {SPECIE_LAT <- "Bos_taurus"
         DATASET = "btaurus_gene_ensembl"},
         Horse = {SPECIE_LAT <- "Equus_caballus"
         DATASET = "ecaballus_gene_ensembl"},
         Cat = {SPECIE_LAT <- "Felis_catus"
         DATASET = "fcatus_gene_ensembl"},
         GuineaPig = {SPECIE_LAT <- "Cavia_porcellus"
         DATASET = "cporcellus_gene_ensembl"},
         Chicken = {SPECIE_LAT <- "Gallus_gallus" 
         DATASET = "ggallus_gene_ensembl"},
         # potentially new species with Bgee 15_0
         Turkey = { SPECIE_LAT <- "Meleagris_gallopavo"
         DATASET = "mgallopavo_gene_ensembl"}, 
         Goat = {SPECIE_LAT <- "Capra_hircus"
         DATASET = "chircus_gene_ensembl"},
         Sheep = {SPECIE_LAT <- "Ovis_aries"
         DATASET = "oaries_gene_ensembl"}, # oarambouillet_gene_ensembl
         Monkey_CrabEatingMacaque = {SPECIE_LAT <- "Macaca_fascicularis"
         DATASET = "mfascicularis_gene_ensembl"}
  )
  # listBgeeSpecies(release = "14.1")
  bgee <- Bgee$new(species = SPECIE_LAT, dataType = "rna_seq") # "affymetrix", "est", "in_situ"
  print(paste0("Collecting data from Bgee, load experimental data set No:", as.character(1)))
  # Combine RNA-seq data from all given experiments
  #DataTableBgee <- data.table::data.table(getData(bgee))
  
  annotation_bgee <- getAnnotation(bgee)
  TMP1 <- tibble(annotation_bgee[[1]]) %>% dplyr::select(Anatomical.entity.name, Anatomical.entity.ID) %>% distinct()
  TMP2 <- tibble(annotation_bgee[[1]]) %>% dplyr::select(Stage.name, Stage.ID) %>% distinct()
  # TMP <- as.data.table(unique(annotation_bgee[[1]]))#[, c("Anatomical.entity.ID","Anatomical.entity.name","Stage.ID","Stage.name")]))
  # TMP[, Specie := as.factor(SPECIE)]
  # unique(annotation_bgee[[2]][, c("Experiment.ID","Experiment.name","Data.source.URL")])
  Organs <- rbind(Organs, TMP1)    
  Ages <- rbind(Ages, TMP2)  
  
}
#write.table(sort(toupper(unique(Organs[, c("Anatomical.entity.ID","Anatomical.entity.name","Stage.ID","Stage.name")]))), file = "BgeeOrgans.txt", row.names = F)
write.table(Organs %>% distinct() %>% dplyr::mutate(Anatomical.entity.name = toupper(Anatomical.entity.name)), file = "BgeeOrgans.txt", row.names = F)
write.table(Ages %>% distinct() %>% dplyr::mutate(Stage.name = toupper(Stage.name)), file = "BgeeAges.txt", row.names = F)

#write.table(Organs, file = "BgeeAnnotations.txt", row.names = F)
# 2.0 Use tissue container information of getAnnotation(bgee) for UBERON: annotations
## paste0(PATH_DB,"TAB_CONTAINER_TISSUE.txt")
## unique(DataTableFinal$Anatomical.entity.name)
