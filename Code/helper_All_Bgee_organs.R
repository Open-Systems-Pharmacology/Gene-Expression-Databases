#### Get information on annotated organs from all species
## This information is needed to map the organ information from data source to the
## PK-Sim internal organ compartment nomenclature (found in TAB_CONTAINER_TISSUE.txt)
## If names are not matched the expression data can not be linked to the PK-Sim organ compartment!
source(paste0(PATH, "/Code/helper_Species.R"))
# ALL_SPECIE <- c(
#  "Dog", "Mouse", "Rat", "Rabbit", "Monkey_mulatta", "Monkey_fascicularis", "Monkey_PigTailed", "Minipig",
#  "Zebrafish", "Cattle", "Horse", "Cat", "GuineaPig", "Chicken", "Sheep", "Goat", "Human"
# )
ALL_SPECIE <- c("Human", PharmaSpecies, AnimalHealthSpecies)
dir.create("BgeeDBs")
setwd("BgeeDBs/")
Organs <- tibble::tibble()
Ages <- Organs

for (SPECIE in ALL_SPECIE) {
  switch(SPECIE,
    #### Clinical specie ####
    Human = {
      SPECIE_LAT <- "Homo_sapiens"
      DATASET <- "hsapiens_gene_ensembl"
    },
    #### Pre-clinical species ####
    Monkey_mulatta = {
      SPECIE_LAT <- "Macaca_mulatta"
      DATASET <- "mmulatta_gene_ensembl"
    },
    Monkey_fascicularis = {
      SPECIE_LAT <- "Macaca_fascicularis"
      DATASET <- "mfascicularis_gene_ensembl"
    },
    Monkey_PigTailed = {
      SPECIE_LAT <- "Macaca_nemestrina"
      DATASET <- "nnemestrina_gene_ensembl"
    },
    Minipig = {
      SPECIE_LAT <- "Sus_scrofa"
      DATASET <- "sscrofa_gene_ensembl"
    },
    Dog = {
      SPECIE_LAT <- "Canis_lupus familiaris"
      DATASET <- "cfamiliaris_gene_ensembl"
    },
    Mouse = {
      SPECIE_LAT <- "Mus_musculus"
      DATASET <- "mmusculus_gene_ensembl"
    },
    Rat = {
      SPECIE_LAT <- "Rattus_norvegicus"
      DATASET <- "rnorvegicus_gene_ensembl"
    },
    Rabbit = {
      SPECIE_LAT <- "Oryctolagus_cuniculus"
      DATASET <- "ocuniculus_gene_ensembl"
    },
    #### Animal health relevant species ####
    Zebrafish = {
      SPECIE_LAT <- "Danio_rerio"
      DATASET <- "drerio_gene_ensembl"
    },
    Cattle = {
      SPECIE_LAT <- "Bos_taurus"
      DATASET <- "btaurus_gene_ensembl"
    },
    Horse = {
      SPECIE_LAT <- "Equus_caballus"
      DATASET <- "ecaballus_gene_ensembl"
    },
    Cat = {
      SPECIE_LAT <- "Felis_catus"
      DATASET <- "fcatus_gene_ensembl"
    },
    GuineaPig = {
      SPECIE_LAT <- "Cavia_porcellus"
      DATASET <- "cporcellus_gene_ensembl"
    },
    Chicken = {
      SPECIE_LAT <- "Gallus_gallus"
      DATASET <- "ggallus_gene_ensembl"
    },
    # potentially new species with Bgee 15_0
    Turkey = {
      SPECIE_LAT <- "Meleagris_gallopavo"
      DATASET <- "mgallopavo_gene_ensembl"
    },
    Goat = {
      SPECIE_LAT <- "Capra_hircus"
      DATASET <- "chircus_gene_ensembl"
    },
    Sheep = {
      SPECIE_LAT <- "Ovis_aries"
      DATASET <- "oaries_gene_ensembl"
    } # oarambouillet_gene_ensembl
  )
  # listBgeeSpecies(release = "14.1")
  bgee <- BgeeDB::Bgee$new(
    species = SPECIE_LAT,
    dataType = "rna_seq"
  ) # "affymetrix", "est", "in_situ"
  print(paste0("Collecting data from Bgee, load experimental data set No:", as.character(1)))
  # Combine RNA-seq data from all given experiments
  # DataTableBgee <- data.table::data.table(getData(bgee))

  annotation_bgee <- BgeeDB::getAnnotation(bgee)
  TMP1 <- tibble::tibble(annotation_bgee[[1]]) |>
    dplyr::select(Anatomical.entity.name, Anatomical.entity.ID)

  TMP2 <- tibble::tibble(annotation_bgee[[1]]) |>
    dplyr::select(Stage.name, Stage.ID)
  # TMP <- as.data.table(unique(annotation_bgee[[1]]))#[, c("Anatomical.entity.ID","Anatomical.entity.name","Stage.ID","Stage.name")]))
  # TMP[, Specie := as.factor(SPECIE)]
  # unique(annotation_bgee[[2]][, c("Experiment.ID","Experiment.name","Data.source.URL")])
  Organs <- rbind(Organs, TMP1) |>
    dplyr::distinct()
  Ages <- rbind(Ages, TMP2) |>
    dplyr::distinct()
}
# write.table(sort(toupper(unique(Organs[, c("Anatomical.entity.ID","Anatomical.entity.name","Stage.ID","Stage.name")]))), file = "BgeeOrgans.txt", row.names = F)
utils::write.table(
  Organs |>
    dplyr::mutate(Anatomical.entity.name = toupper(Anatomical.entity.name)) |>
    dplyr::arrange(Anatomical.entity.name) |>
    dplyr::distinct(),
  file = "BgeeOrgans.txt",
  row.names = FALSE,
  sep = ";",
  quote = FALSE
)
utils::write.table(
  Ages |>
    dplyr::mutate(Stage.name = toupper(Stage.name)) |>
    dplyr::arrange(Stage.ID) |>
    dplyr::distinct(),
  file = "BgeeAges.txt",
  row.names = FALSE,
  sep = ";",
  quote = FALSE
)
setwd("../")
# write.table(Organs, file = "BgeeAnnotations.txt", row.names = F)
# 2.0 Use tissue container information of getAnnotation(bgee) for UBERON: annotations
## paste0(PATH_DB,"TAB_CONTAINER_TISSUE.txt")
## unique(DataTableFinal$Anatomical.entity.name)
