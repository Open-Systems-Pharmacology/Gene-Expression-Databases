#### Get information on annotated organs from all species
## This information is needed to map the organ information from data source to the
## PK-Sim internal organ compartment nomenclature (found in TAB_CONTAINER_TISSUE.txt)
## If names are not matched the expression data can not be linked to the PK-Sim organ compartment!

#' Build Bgee organ and age lookup tables for all species
#'
#' Downloads Bgee RNA-Seq annotation for all species and writes
#' BgeeOrgans.txt and BgeeAges.txt to BgeeDBs/. These files must exist
#' before GeneratePKsimDB can map tissues to PK-Sim container names.
#' Results are cached; re-running overwrites the files.
#'
#' @param PATH Character. Absolute path to the repository root.
#' @return Invisibly NULL. Side effect: writes BgeeOrgans.txt and BgeeAges.txt.
build_bgee_lookup_tables <- function(PATH) {
  source(paste0(PATH, "/Code/03-helpers/helper_Species.R"))

  ALL_SPECIE_LOCAL <- c("Human", PharmaSpecies, AnimalHealthSpecies)

  dir.create(file.path(PATH, "BgeeDBs"), recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd(file.path(PATH, "BgeeDBs"))
  on.exit(setwd(old_wd), add = TRUE)


  Ages <- Organs

  for (SPECIE in ALL_SPECIE_LOCAL) {
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
        DATASET <- "mnemestrina_gene_ensembl"
      },
      Minipig = {
        SPECIE_LAT <- "Sus_scrofa"
        DATASET <- "sscrofa_gene_ensembl"
      },
      Dog = {
        SPECIE_LAT <- "Canis_lupus_familiaris"
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
      Guineapig = {
        SPECIE_LAT <- "Cavia_porcellus"
        DATASET <- "cporcellus_gene_ensembl"
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
      Chicken = {
        SPECIE_LAT <- "Gallus_gallus"
        DATASET <- "ggallus_gene_ensembl"
      },
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
      }, # oarambouillet_gene_ensembl
      stop(paste("Unknown SPECIE:", SPECIE))
    )

    bgee <- BgeeDB::Bgee$new(species = SPECIE_LAT, dataType = "rna_seq")
    print(paste0("Collecting Bgee annotation for: ", SPECIE))

    annotation_bgee <- BgeeDB::getAnnotation(bgee)
    TMP1 <- tibble::tibble(annotation_bgee[[1]]) |>
      dplyr::select(Anatomical.entity.name, Anatomical.entity.ID)

    TMP2 <- tibble::tibble(annotation_bgee[[1]]) |>
      dplyr::select(Stage.name, Stage.ID)

    Organs <- rbind(Organs, TMP1) |> dplyr::distinct()
    Ages   <- rbind(Ages,   TMP2) |> dplyr::distinct()
  }

  utils::write.table(
    Organs |>
      dplyr::mutate(Anatomical.entity.name = toupper(Anatomical.entity.name)) |>
      dplyr::arrange(Anatomical.entity.name) |>
      dplyr::distinct(),
    file      = "BgeeOrgans.txt",
    row.names = FALSE,
    sep       = ";",
    quote     = FALSE,
    append    = FALSE
  )
  utils::write.table(
    Ages |>
      dplyr::mutate(Stage.name = toupper(Stage.name)) |>
      dplyr::arrange(Stage.ID) |>
      dplyr::distinct(),
    file      = "BgeeAges.txt",
    row.names = FALSE,
    sep       = ";",
    quote     = FALSE,
    append    = FALSE
  )
  invisible(NULL)
}

