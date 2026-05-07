# Suppress dplyr global variable binding warnings
utils::globalVariables(c(
  "ensembl_gene_id", "entrezgene_id", "external_gene_name", "uniprot_gn_symbol", "description", "external_synonym", "wikigene_name",
  "hsapiens_homolog_associated_gene_name", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_perc_id_r1",
  "VARIANT_NAME", "ENTREZID", "PREFFERED_NAME", "SYMBOL", "OFFICIAL_FULL_NAME", "SYNONYM", "OTHER_NAME",
  "HOMOLOG_SYMBOL", "HOMOLOG", "HOMOLOGY_PERCENT"
))
#' PrepareBioMarts
#'
#' Description:
#'   This function downloads and prepares gene annotation data for a given species
#'   using Ensembl BioMart and other sources. It stores the information as RDS and
#'   database files, enabling downstream RNA-Seq analysis, gene identifier mapping,
#'   and ortholog mapping between species. The function also filters for ADME-relevant
#'   genes and supports lazy tibble operations for efficient data handling.
#'
#' Parameters:
#'   SPECIE  -  Character. Species name (e.g., "Rat", "Human", "Dog").
#'              Must be one of the supported species in the function.
#'
#' Returns:
#'   No return value. Writes annotation tables to a local SQLite database.
#'
#' Example:
#'   PrepareBioMarts(SPECIE = "Human")
#'   PrepareBioMarts(SPECIE = "Rat")
#'
#' Details:
#'   - Estimates gene length from chromosome coordinates.
#'   - Maps multiple gene identifiers (Ensembl, NCBI, symbol, etc.).
#'   - Maps orthologs to human genes (where available).
#'   - Filters for ADME genes using curated lists from multiple sources.
#'   - Stores results in a database for efficient downstream analysis.
#'
PrepareBioMarts <- function(SPECIE = "Rat") {
  # Function is designed to download species gene annotation and store the information (as RDS & DB files)
  # The information is needed to:
  #  1. Estimate RPKM from other RNAseq data (gene length needed; https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)
  #  2. Annotate a gene to enable search over multiple identifiers
  #     (e.g. NCBI gene id: 1543, ensembl id: ENSG00000140465, symbol: CYP1A1, etc.)
  #  3. Map ortholog genes (homologs between species)
  #
  # Data storage ideally in an data base to enable lazy tibble operations
  # List of human to ortholog animal ensembl id mapping (downloaded from: https://www.ensembl.org/biomart/martview)
  # Length of a gene is estimated form start and stop region on chromosome
  # A subset of genes is prepared including only AMDE relevant genes (reduce data load)
  #
  # ADME gene subset is defined via gene symbols compiled from different sources:
  #  1. from: https://www.physiome.org/Course/2015winBioen498e/wk4/KellDobson08.pdf
  #  2. pharmaadme.org
  #  3. http://www.membranetransport.org/
  #  4. http://www.tcdb.org/
  #
  #
  # Notably Sheep homologous genes are not available from oaries_gene_ensembl
  # Alternative homology mappings via: OMA # https://bioconductor.org/packages/release/bioc/html/OmaDB.html
  #### Preparations ####
  # adding the various identifier and information is extremely memory
  # intensive and might only work with a 64-bit R version!
  source(paste0(PATH, "/Code/03-helpers/helper_Species.R"))
  # ALL_SPECIE <- c(
  #    "Mouse", "Rat", "Rabbit", "Guineapig", "Dog", "Minipig",
  #    "Monkey_mulatta", "Monkey_fascicularis", "Monkey_PigTailed",
  #    "Cattle", "Horse", "Cat", "Chicken", "Goat", "Sheep", "Turkey",
  #    "Zebrafish", "Human"
  # )

  if (!SPECIE %in% ALL_SPECIE) {
    stop(paste0(
      "Given SPECIE: '",
      SPECIE, "' not in list. \n Must be one of: ",
      paste(ALL_SPECIE, collapse = " , ")
    ))
  }

  # Create annotation database, if exist only connection is made
  dir.create("./BioMarts/", showWarnings = FALSE)

  Conn <- DBI::dbConnect(
    drv = RSQLite::SQLite(),
    paste0("./BioMarts/All_Species_BioMarts.DB"),
    synchronous = NULL
  )
  if (!DBI::dbExistsTable(conn = Conn, name = "Human_ADME") && SPECIE != "Human") {
    print(paste0("Human annotation do not exist in 'All_Species_BioMarts.DB' \n Human data is loaded first."))
    PrepareBioMarts("Human")
  }

  #### Load biomaRt data and store ####
  print(paste0("Loading ", SPECIE, " gene inforamtion from biomart mirror."))
  # List of available species:
  # ensembl <- biomaRt::useMart("ensembl")
  # ListOfEnsemblSpecies <- biomaRt::listDatasets(ensembl)
  # ListOfEnsemblArchives <- biomaRt::listEnsemblArchives()
  #
  switch(SPECIE,
    Cat = {
      ensembl <- biomaRt::useMart("ensembl",
        dataset = "fcatus_gene_ensembl",
        host = "https://may2021.archive.ensembl.org"
      )
    },
    Cattle = {
      ensembl <- biomaRt::useMart("ensembl",
        dataset = "btaurus_gene_ensembl",
        host = "https://may2021.archive.ensembl.org"
      )
    },
    Chicken = {
      ensembl <- biomaRt::useMart("ensembl",
        dataset = "ggallus_gene_ensembl",
        host = "https://apr2022.archive.ensembl.org"
      )
    },
    Dog = {
      ensembl <- biomaRt::useMart("ensembl",
        dataset = "clfamiliaris_gene_ensembl",
        host = "https://may2021.archive.ensembl.org"
      )
    },
    Goat = {
      ensembl <- biomaRt::useMart("ensembl",
        dataset = "chircus_gene_ensembl"
      )
    },
    Guineapig = {
      ensembl <- biomaRt::useMart("ensembl",
        dataset = "cporcellus_gene_ensembl",
        host = "https://may2025.archive.ensembl.org"
      )
    },
    Horse = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "ecaballus_gene_ensembl")
    },
    Human = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    },
    Minipig = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "sscrofa_gene_ensembl")
    },
    Monkey_fascicularis = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "mfascicularis_gene_ensembl")
    },
    Monkey_mulatta = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "mmulatta_gene_ensembl")
    },
    Monkey_PigTailed = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "mnemestrina_gene_ensembl")
    },
    Mouse = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    },
    Rabbit = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "ocuniculus_gene_ensembl")
    },
    Rat = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
    },
    Sheep = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "oaries_gene_ensembl")
    },
    Turkey = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "mgallopavo_gene_ensembl")
    },
    Zebrafish = {
      ensembl <- biomaRt::useMart("ensembl", dataset = "drerio_gene_ensembl", host = "https://oct2024.archive.ensembl.org")
    }
  )

  #### Select annotation information and identifiers; mapping ortholog genes  ####
  print(paste0("Fetiching ", SPECIE, " gene identifiers and annotations."))

  # Get SPECIE specific gene annotations
  SPECIE_ANNOTATION <- biomaRt::getBM(attributes = c(
    "ensembl_gene_id", "description", "entrezgene_id",
    "external_synonym", "external_gene_name", "wikigene_name",
    "uniprot_gn_symbol"
  ), mart = ensembl) # hgnc_symbol

  # Add SPECIE specific gene lengths and map human ortholog onto animal genes
  SPECIE_ANNOTATION <- dplyr::left_join(SPECIE_ANNOTATION,
    biomaRt::getBM(
      attributes = c("ensembl_gene_id", "start_position", "end_position"),
      mart = ensembl
    ),
    by = "ensembl_gene_id"
  )
  SPECIE_ANNOTATION <-
    dplyr::rename(
      .data = SPECIE_ANNOTATION,
      VARIANT_NAME = ensembl_gene_id,
      ENTREZID = entrezgene_id,
      PREFFERED_NAME = external_gene_name,
      SYMBOL = uniprot_gn_symbol,
      OFFICIAL_FULL_NAME = description,
      SYNONYM = external_synonym,
      OTHER_NAME = wikigene_name
    )

  if (SPECIE != "Human" && SPECIE != "Sheep") {
    print(paste0("Link ", SPECIE, " gene identifiers to human orthologs based on ensembl IDs."))
    SPECIE_ANNOTATION <- dplyr::left_join(SPECIE_ANNOTATION,
      biomaRt::getBM(
        attributes = c(
          "ensembl_gene_id", "hsapiens_homolog_associated_gene_name",
          "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_perc_id_r1"
        ),
        mart = ensembl
      ),
      by = c("VARIANT_NAME" = "ensembl_gene_id")
    ) # "hsapiens_homolog_orthology_type"

    SPECIE_ANNOTATION <- dplyr::rename(
      .data = SPECIE_ANNOTATION,
      HOMOLOG_SYMBOL = hsapiens_homolog_associated_gene_name,
      HOMOLOG = hsapiens_homolog_ensembl_gene,
      HOMOLOGY_PERCENT = hsapiens_homolog_perc_id_r1
    )
  }
  # String cleaning
  SPECIE_ANNOTATION <- tibble::tibble(SPECIE_ANNOTATION) |>
    dplyr::arrange(VARIANT_NAME)
  SPECIE_ANNOTATION <- SPECIE_ANNOTATION |>
    dplyr::mutate(
      OFFICIAL_FULL_NAME =
        stringr::str_remove(string = OFFICIAL_FULL_NAME, pattern = "\\[Source:.*")
    )

  # Write annotation table to database
  DBI::dbWriteTable(
    conn = Conn,
    name = paste0(SPECIE, "_Annotations"),
    value = tibble::as_tibble(SPECIE_ANNOTATION),
    overwrite = TRUE
  )

  #### Only ADME genes ####
  # List of ADME gene abbreviations
  ADMEgene <- c(
    "^ABC", "^ADH", "^AHR", "^ALD", "^ALDH", "^ANXA", "^AOX", "^ARN", "^ARNT",
    "^ARS", "^ATP", "^CACN", "^CAT", "^CBR", "^CDA", "^CES", "^CFT", "^CHS",
    "^CHST", "^CYB", "^CYP", "^DDO", "^DHR", "^DHRS", "^DPE", "^DPY", "^DPYD",
    "^EPH", "^EPHX", "^FMO", "^GPX", "^GSR", "^GSS", "^GST", "^HAG", "^HNF",
    "^HNM", "^HNMT", "^HSD", "^IAP", "^KCN", "^MAT", "^MET", "^MGS", "^MGST",
    "^MPO", "^NAT", "^NNM", "^NNMT", "^NOS", "^NR1", "^NUDT", "^PDE", "^PIAS",
    "^PLG", "^PNM", "^PNMT", "^PON", "^POR", "^PPA", "^PPAR", "^PPARA",
    "^RXR", "^SER", "^SLC", "^SCN", "^SOD", "^SUL", "^SULF", "^SULT", "^SQST",
    "^TAP", "^TPM", "^TPMT", "^UGT", "^UROC", "^VKOR", "^XDH", "^XRC" # ,"^LOC"
  )

  # Select only ADME genes
  if (SPECIE != "Human" && SPECIE != "Sheep") {
    SPECIE_ANNOTATION <- SPECIE_ANNOTATION |>
      dplyr::filter(grepl(paste(ADMEgene, collapse = "|"), SYNONYM) |
        grepl(paste(ADMEgene, collapse = "|"), SYMBOL) |
        grepl(paste(ADMEgene, collapse = "|"), HOMOLOG_SYMBOL) |
        grepl(paste(ADMEgene, collapse = "|"), PREFFERED_NAME) |
        grepl(paste(ADMEgene, collapse = "|"), OTHER_NAME)) |>
      dplyr::distinct() |>
      dplyr::arrange(VARIANT_NAME)
  } else {
    SPECIE_ANNOTATION <- SPECIE_ANNOTATION |>
      dplyr::filter(grepl(paste(ADMEgene, collapse = "|"), SYNONYM) |
        grepl(paste(ADMEgene, collapse = "|"), SYMBOL) |
        grepl(paste(ADMEgene, collapse = "|"), PREFFERED_NAME) |
        grepl(paste(ADMEgene, collapse = "|"), OTHER_NAME)) |>
      dplyr::distinct() |>
      dplyr::arrange(VARIANT_NAME)
  }
  # Save to database, compress database and close connection
  DBI::dbWriteTable(
    conn = Conn,
    name = paste0(SPECIE, "_ADME"),
    value = tibble::as_tibble(SPECIE_ANNOTATION),
    overwrite = TRUE
  )
  DBI::dbExecute(conn = Conn, statement = "VACUUM") # reduce storage space
  DBI::dbDisconnect(Conn)
}
