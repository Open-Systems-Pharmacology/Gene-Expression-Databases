###############################################################
# GeneratePKsimDB
#
# Description:
#   This function generates a PK-Sim compatible gene expression database
#   for a given species using Bgee RNA-Seq data. It supports filtering for
#   ADME genes, different Bgee releases, and normalization options.
#   The output is stored in the specified path and can be computed in RAM
#   for performance.
#
# Parameters:
#   SPECIE        - Character. Species name (e.g., "Dog", "Human").
#   PATH          - Character. Output directory path. Default: current working directory.
#   ADME_ONLY     - Logical. If TRUE, only ADME genes are included. Default: TRUE.
#   RELEASE       - Character. Bgee release version (e.g., "15_2"). Default: "15_2".
#   COMPUTE_IN_RAM- Logical. If TRUE, computation is performed in RAM. Default: FALSE.
#   TPM_ONLY      - Logical. If TRUE, only TPM normalized data is used. Default: TRUE.
#   INCLUDE_RPKM  - Logical. If TRUE, RPKM values are estimated. Default: FALSE.
#
# Returns:
#   No return value. Writes output files to disk.
#
# Example:
#   GeneratePKsimDB(SPECIE = "Human", PATH = "./output", ADME_ONLY = TRUE)
###############################################################
GeneratePKsimDB <- function(
    SPECIE = "Rabbit",
    PATH = getwd(),
    ADME_ONLY = TRUE,
    RELEASE = "15_2",
    COMPUTE_IN_RAM = FALSE,
    TPM_ONLY = TRUE,
    INCLUDE_RPKM = FALSE) {
  ### Dependencies ####
  switch(RELEASE,
    "13_2" = {
      simpleError("Code was discontiniued for Bgee release 13_2, try >= 15_0")
    },
    "14_0" = {
      simpleError("Code was discontiniued for Bgee release 14_0, try >= 15_0")
    },
    "14_1" = {
      simpleError("Code was discontiniued for Bgee release 14_1, try >= 15_0")
    },
    "14_2" = {
      simpleError("Code was discontiniued for Bgee release 14_2, try >= 15_0")
    },
    "15_0" = {
      simpleMessage("Code was tested for for Bgee release 15_0")
    },
    "15_1" = {
      simpleMessage("Code was not tested for for Bgee release 15_1")
    },
    "15_2" = {
      simpleMessage("Code was tested for for Bgee release 15_2")
    },
    {
      simpleMessage(paste0("Input Bgee release: ", RELEASE, " was not recognized. Default: 15_2 is used instead"))
      RELEASE <- "15_2"
    }
  )

  # Test if input species is valid
  ALL_SPECIE <- c(
    "Human", "Monkey_mulatta", "Minipig", "Dog", "Mouse", "Rat", "Rabbit",
    "Zebrafish", "Cattle", "Horse", "Cat", "GuineaPig", "Chicken",
    "Goat", "Sheep", "Turkey", "Monkey_fascicularis", "Monkey_PigTailed"
  )

  if (!SPECIE %in% ALL_SPECIE) {
    stop(paste0(
      "Given SPECIE: '", SPECIE, "' not in list. \n Must be one of: ",
      paste(ALL_SPECIE, collapse = " , ")
    ))
  }

  # load SQL commands to set views and column types
  source(paste0(PATH, "/Code/helper_SQL_Commands.R"))

  # is needed to allow download of human data
  # (65 GB take longer than 1 min ;-) )
  options(timeout = 60 * 60 * 15)

  #### Load experimental data from bgee ####
  # listBgeeSpecies(ordering = 1)
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
    # new species with Bgee > 15_0
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
    }
  )
  print("Fetch Bgee expression data sets")
  # depending on RELEASE
  if (is.na(RELEASE)) { # if specific release is desired
    bgee <- BgeeDB::Bgee$new(
      species = SPECIE_LAT,
      dataType = "rna_seq",
      sendStats = FALSE
    ) # other available data domains "affymetrix", "est", "in_situ"
  } else {
    bgee <- BgeeDB::Bgee$new(
      species = SPECIE_LAT,
      dataType = "rna_seq",
      sendStats = FALSE,
      release = RELEASE
    )
  }

  # Set naming of data bases; listBgeeRelease()[1,"release"]
  BegeeRelease <- bgee$release
  DB_Bgee <- paste0(bgee[["species"]], bgee[["sqlite_extension"]])
  DB_PKsim <- paste0("GENEDB_", tolower(SPECIE))
  if (ADME_ONLY) {
    DB_PKsim <- paste0(DB_PKsim, "_ADME_ONLY")
  }
  if (TPM_ONLY) {
    DB_PKsim <- paste0(DB_PKsim, "_TPM_ONLY")
  }
  DB_PKsim <- paste0(
    DB_PKsim, "_BgeeRelease_",
    BegeeRelease, "_", Sys.Date(), ".expressionDB"
  )

  PATH2PKsim <- paste0(PATH, "/", SPECIE, "/")
  PATH_DB_Bgee <- bgee[["pathToData"]]

  dir.create(PATH2PKsim, showWarnings = FALSE)

  #### Connections to local data bases ####
  print("Build local data bases")
  # local bgee DB
  db_bgee_conn <-
    DBI::dbConnect(RSQLite::SQLite(),
      paste0(PATH_DB_Bgee, "/", gsub(x = DB_Bgee, pattern = " ", replacement = "_")),
      synchronous = "off",
      cache_size = -1000
    )

  DB_Tables <- DBI::dbListTables(conn = db_bgee_conn)

  # local biomart DB
  db_biomart_conn <-
    DBI::dbConnect(RSQLite::SQLite(),
      "BioMarts/All_Species_BioMarts.DB",
      synchronous = "off",
      cache_size = -1000
    )

  # local pkSim DB
  db_PKsim_conn <- DBI::dbConnect(RSQLite::SQLite(),
    paste0(PATH2PKsim, DB_PKsim),
    synchronous = "off",
    cache_size = -1000
  )

  # If tables are empty they have been created but not filled with data
  if (rlang::is_empty(DB_Tables)) {
    # If data is not yet locally stored execute getSampleProcessedData()
    # to initiate download, currently the latest version
    # getSampleProcessedData() switched to a forced load of whole
    # db content into ram. This can cause to errors when the function
    # GeneratePKsimDB.R is called the first time, re-executing the function
    # usually resolves the issue (then DB is already downloaded).
    if (SPECIE_LAT == "Homo_sapiens") {
      # Access expression data table, either in local storage or as remote
      # lazy table (connection to database). Human data is implemented to
      # not include SRP012682 "Genotype-Tissue Expression (GTEx) Common Fund
      # Project. single data source with > 60 GB!
      # code can currently not process this amount of data... :-(
      rna_seq_selected <- BgeeDB::getSampleProcessedData(bgee,
        experimentId = # SRP058036
          c(
            "GSE57344", "GSE62098", "GSE43520", "GSE30352",
            "GSE58387", "GSE58608", "GSE64283", "GSE30611"
          )
        # RNA-Seq experiments sorted by number of organs/tissue count 75 -> 4
        # c("SRP012682", "ERP003613", "SRP163252", "GSE30611", "SRP043364",
        #   "ERP109002", "ERP006650", "SRP111096", "GSE30352", "GSE64283",
        #   "SRP058036", "SRP007359", "SRP028336", "SRP102989", "GSE43520",
        #   "ERP120553", "SRP145002", "SRP262331", "SRP058740", "SRP136499")
      )
    } else {
      rna_seq_selected <- BgeeDB::getSampleProcessedData(bgee)
      # The default Bgee databases contain quotes causing error during later
      # database generation, remove to reduce free disc space and for clear naming
      DBI::dbExecute(
        conn = db_bgee_conn,
        statement = c("UPDATE rna_seq SET \"Anatomical.entity.name\" = REPLACE(\"Anatomical.entity.name\", '\"' , '' )")
      )
      DBI::dbExecute(
        conn = db_bgee_conn,
        statement = c("UPDATE rna_seq SET \"Stage.name\" = REPLACE(\"Stage.name\",'\"' , '' )")
      )
      DBI::dbExecute(
        conn = db_bgee_conn,
        statement = c("UPDATE rna_seq SET \"Stage.name\" = REPLACE(\"Stage.name\",' (human)' , '' )")
      )
      DBI::dbExecute(
        conn = db_bgee_conn,
        statement = c("UPDATE rna_seq SET \"Strain\" = REPLACE(\"Strain\", '\"' , '' )")
      )
      # (conn = db_bgee_conn, statement = c("UPDATE rna_seq SET \"Anatomical.entity.ID\" = REPLACE(\"Anatomical.entity.ID\", 'UBERON:' , '' )" )   )
      # DBI::dbExecute(conn = db_bgee_conn, statement = c("UPDATE rna_seq SET \"Anatomical.entity.ID\" = REPLACE(\"Anatomical.entity.ID\", 'CL:' , '' )" )   )
      # DBI::dbExecute(conn = db_bgee_conn, statement = c("UPDATE rna_seq SET \"Stage.ID\" = REPLACE(\"Stage.ID\", 'UBERON:' , '' )" )   )
      # DBI::dbExecute(conn = db_bgee_conn, statement = c("UPDATE rna_seq SET \"Stage.ID\" = REPLACE(\"Stage.ID\", 'HsapDv:' , '' )" )   )
      # DBI::dbExecute(conn = db_bgee_conn, statement = c("VACUUM" )   )
      if (readr::parse_number(BegeeRelease) >= 14) {
        # strain column only include in bgee > 14
        DBI::dbExecute(
          conn = db_bgee_conn,
          statement = c("UPDATE rna_seq SET \"Strain\" = REPLACE(\"Strain\", '\"' , '' )")
        )
      }
    }
  }

  print("Re-organize data from BgeeDB for PK-Sim compatibility")
  if (COMPUTE_IN_RAM) {
    # as local table in RAM
    print("Fetch expression data from local BgeeDB.")
    rna_seq_selected <- BgeeDB::getSampleProcessedData(bgee)
  } else {
    # as remote table
    rna_seq_selected <- dplyr::tbl(db_bgee_conn, "rna_seq")
    # Human data sets can be quite large and cause the code to crash on
    # systems with limited RAM
  }

  # Subset only high quality and as present declared genes form
  # (inconsistent over different Bgee releases)
  if (readr::parse_number(BegeeRelease) >= 15) {
    rna_seq_selected <- rna_seq_selected |>
      dplyr::filter(
        Detection.flag == "present",
        Anatomical.entity.name != "multicellular organism"
        # , Stage.name != "life cycle"
      ) |>
      dplyr::select(
        !tidyselect::any_of(c(
          "Anatomical.entity.ID",
          "Library.type",
          "Stage.ID"
        ))
      ) |>
      dplyr::select(
        !tidyselect::any_of(c(
          "Detection.flag",
          "State.in.Bgee",
          "pValue",
          "Rank",
          "FPKM",
          "Read.count"
        ))
      )
  } else {
    rna_seq_selected <- rna_seq_selected |>
      dplyr::filter(Detection.flag == "present") |>
      dplyr::filter(Detection.quality == "high quality") |>
      dplyr::select(!tidyselect::any_of(c(
        "Anatomical.entity.ID",
        "Library.type", "Stage.ID"
      ))) |>
      dplyr::select(!tidyselect::any_of(c(
        "Detection.flag",
        "Detection.quality",
        "State.in.Bgee"
      )))
  }

  if (readr::parse_number(BegeeRelease) < 14) {
    rna_seq_selected <- rna_seq_selected |>
      mutate(Sex = "UNSPECIFIED") |>
      dplyr::mutate(Strain = SPECIE) |>
      dplyr::mutate(FPKM = NA) |>
      dplyr::mutate(TPM = NA)
  }

  # Rename columns to match PK-Sim required namespace
  rna_seq_selected <- dplyr::rename(
    .data = rna_seq_selected,
    data_base = Experiment.ID,
    data_base_rec_id = Library.ID,
    variant_name = Gene.ID,
    tissue = Anatomical.entity.name,
    gender = Sex,
    strain = Strain,
    age = Stage.name
  ) #|> dplyr::compute()

  # Further format, rename and prepare data ####
  rna_seq_selected <- rna_seq_selected |>
    dplyr::mutate(gender = tolower(gender)) |>
    dplyr::mutate(age = tolower(age)) |>
    dplyr::mutate(strain = tolower(strain)) |>
    dplyr::mutate(gender = dplyr::if_else(gender == tolower("NA"),
      tolower("UNSPECIFIED"), gender
    )) |>
    dplyr::mutate(strain = dplyr::if_else(strain == tolower("NA"),
      tolower("UNSPECIFIED"), strain
    )) |>
    dplyr::mutate(tissue = dplyr::if_else(tissue == tolower("NA"),
      tolower("UNSPECIFIED"), tissue
    )) |>
    dplyr::mutate(tissue = tolower(tissue)) |>
    dplyr::mutate(state = tolower("NORMAL")) |>
    dplyr::mutate(health_state = tolower(paste(strain, state, sep = " "))) |>
    dplyr::compute(name = "rna_seq_selected", temporary = FALSE)

  ##### Add Gene IDs / Annotations needed for allow broad reach ####
  print("Annotete expression data with gene & protein identifiers.")
  # ADME only genes
  AnnotationTable <-
    DBI::dbReadTable(db_biomart_conn, paste0(SPECIE, "_ADME")) |>
    dplyr::select(-c("start_position", "end_position"))
  # For only homologous genes, add "HOMOLOG" prefix in Gene Symbol;
  # This is necessary since PK-SIM gene names are based on Symbols
  AnnotationTable <- AnnotationTable |>
    dplyr::collect() |>
    readr::type_convert() |>
    dplyr::mutate(OFFICIAL_FULL_NAME = trimws(OFFICIAL_FULL_NAME))
  colnames(AnnotationTable) <- tolower(colnames(AnnotationTable))

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "AnnotationTable_ADME",
    value = AnnotationTable,
    row.names = FALSE,
    overwrite = TRUE
  )

  AnnotationTable <-
    DBI::dbReadTable(db_biomart_conn, paste0(SPECIE, "_Annotations")) |>
    dplyr::select(-c("start_position", "end_position"))
  AnnotationTable <- AnnotationTable |>
    dplyr::collect() |>
    readr::type_convert() |>
    dplyr::mutate(OFFICIAL_FULL_NAME = trimws(OFFICIAL_FULL_NAME))
  colnames(AnnotationTable) <- tolower(colnames(AnnotationTable))

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "AnnotationTable",
    value = AnnotationTable,
    row.names = FALSE,
    overwrite = TRUE
  )

  #### generate TAB tables from data base data frame ####
  print("Write data to TAB files")
  # Currently realizing the query in local RAM (collect) needed for
  # copying accross database is the ratelimiting step.
  # Try sqlpipe for data transfer between local DBs without detour through local RAM:
  # https://github.com/sqlpipe/sqlpipe?tab=readme-ov-file
  # icacls sqlpipe /grant USER:RX

  # Combine expression data and annotation information
  if (ADME_ONLY) {
    AnnotationTable <- dplyr::tbl(db_bgee_conn, "AnnotationTable_ADME")
    rna_seq_selected <- dplyr::left_join(
      AnnotationTable |>
        dplyr::select(variant_name) |>
        dplyr::distinct(),
      by = dplyr::join_by(variant_name),
      rna_seq_selected
    ) |>
      dplyr::filter(!is.na(age)) |>
      dplyr::mutate(age = as.character(age)) |>
      dplyr::compute(name = "rna_seq_selected_ADME", temporary = FALSE)
    rna_seq_selected <- dplyr::tbl(db_bgee_conn, "rna_seq_selected_ADME")
  } else {
    AnnotationTable <- dplyr::tbl(db_bgee_conn, "AnnotationTable")
    rna_seq_selected <- dplyr::left_join(
      AnnotationTable |>
        dplyr::select(variant_name) |>
        dplyr::distinct(),
      by = dplyr::join_by(variant_name),
      rna_seq_selected
    ) |>
      dplyr::filter(!is.na(age)) |>
      dplyr::mutate(age = as.character(age)) |>
      dplyr::compute(name = "rna_seq_selected", temporary = FALSE)
    rna_seq_selected <- dplyr::tbl(db_bgee_conn, "rna_seq_selected")
  }

  # Remove rna_seq_selected entries where age is NA
  rna_seq_selected <- rna_seq_selected |>
    dplyr::filter(!is.na(age)) |>
    dplyr::mutate(age = as.character(age))

  # Add gene IDs
  print("Build tab_gene_variants and write to database")
  tab_gene_variants <- rna_seq_selected |>
    dplyr::select(variant_name) |>
    dplyr::distinct() |>
    # dplyr::collect() |>
    dplyr::mutate(variant_id = dplyr::row_number()) |>
    dplyr::mutate(gene_id = variant_id) |>
    dplyr::relocate(gene_id) |>
    dplyr::arrange(variant_name) |>
    dplyr::compute(
      name = "tab_gene_variants",
      temporary = FALSE,
      overwrite = TRUE
    )

  # Add data source IDs
  print("Build tab_expression_data_records_tmp and write to database")
  tab_expression_data_records_tmp <- rna_seq_selected |>
    dplyr::select(data_base_rec_id) |>
    dplyr::distinct() |>
    # dplyr::collect() |>
    dplyr::mutate(data_source_id = dplyr::row_number()) |>
    dplyr::compute(
      name = "tab_expression_data_records_tmp",
      temporary = FALSE,
      overwrite = TRUE
    )

  # Add data age IDs
  print("Build tab_expression_data_ages_tmp and write to database")
  tab_expression_data_ages_tmp <- rna_seq_selected |>
    dplyr::select(age) |>
    dplyr::distinct() |>
    # dplyr::collect() |>
    dplyr::mutate(age = age) |>
    dplyr::mutate(age_id = dplyr::row_number()) |>
    dplyr::arrange(age_id, age) |>
    dplyr::relocate(age_id) |>
    dplyr::compute(
      name = "tab_expression_data_ages_tmp",
      temporary = FALSE,
      overwrite = TRUE
    )

  #### Merge different expression measures (later Units in PK-Sim DB) ####
  # Here cut offs could be considered
  print("Merge expression data from multiple units.")
  tpm_table <- rna_seq_selected |>
    dplyr::select(!tidyselect::any_of(c("FPKM", "Read.count", "RPKM"))) |> # dplyr::distinct() |>
    # dplyr::mutate(unit = "TPM") |>
    # dplyr::rename(sample_count = TPM) |>
    # dplyr::filter(!is.na(sample_count)) |>
    # filter(sample_count > CUT_OFF[["TPM"]]) |>
    # dplyr::mutate(data_base_rec_id = paste0(data_base_rec_id, "_TPM")) |>
    dplyr::group_by(data_base_rec_id) |>
    dplyr::summarise(total_count = sum(TPM, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    # dplyr::collect() |>
    dplyr::mutate(total_count = total_count / 10^6) |>
    dplyr::compute(
      name = "tpm_table",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_properties #
  print("Build tab_expression_data_properties and write to database")
  KEYS <- c(
    "data_source_id", "data_base", "tissue", "health_state",
    "gender", "age"
  )
  TMP <- dplyr::tbl(db_bgee_conn, "tab_expression_data_records_tmp")
  tab_expression_data_properties_tmp <-
    dplyr::left_join(rna_seq_selected, TMP, copy = TRUE) |>
    dplyr::select(tidyselect::all_of(KEYS)) |>
    dplyr::distinct() |>
    dplyr::compute(
      name = "tab_expression_data_properties_tmp",
      temporary = FALSE,
      overwrite = TRUE
    )

  tab_expression_data_properties <-
    tab_expression_data_properties_tmp |>
    tidyr::pivot_longer(!data_source_id,
      names_to = "property",
      values_to = "property_value"
    ) |>
    dplyr::distinct() |>
    dplyr::compute(
      name = "tab_expression_data_properties",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_values #
  print("Build tab_expression_data_values and write to database")
  KEYS <- c("variant_id", "data_source_id", "sample_count", "total_count", "unit")
  tab_expression_data_values <- dplyr::left_join(
    rna_seq_selected, # |> dplyr::mutate(data_base_rec_id = paste0(data_base_rec_id,"_TPM"))
    dplyr::tbl(db_bgee_conn, "tpm_table")
  ) |>
    dplyr::mutate(unit = "TPM") |>
    dplyr::rename(sample_count = TPM) |>
    dplyr::left_join(dplyr::tbl(db_bgee_conn, "tab_gene_variants")) |> # , copy = TRUE) |>
    dplyr::left_join(dplyr::tbl(db_bgee_conn, "tab_expression_data_records_tmp")) |> # , copy = TRUE) |> # |> dplyr::mutate(data_base_rec_id = paste0(data_base_rec_id,"_TPM")
    dplyr::select(tidyselect::all_of(KEYS)) |>
    dplyr::distinct() |>
    dplyr::compute(
      name = "tab_expression_data_values",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_age_properties #
  print("Build tab_expression_data_age_properties and write to database")
  # optional age is extracted form the data sets, in bgee strings like "2-3 weeks old rat" are given
  tab_expression_data_age_properties <- tab_expression_data_ages_tmp |>
    dplyr::select(age_id, age) |>
    dplyr::distinct() |>
    tidyr::pivot_longer(!age_id,
      names_to = "property",
      values_to = "property_value"
    ) |>
    dplyr::compute(
      name = "tab_expression_data_age_properties",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_ages #
  print("Build tab_expression_data_ages and write to database")
  # for mice also Theiler Stage could be included: https://www.emouseatlas.org/emap/ema/theiler_stages/StageDefinition/stagedefinition.html#dpc
  #  Numextract <- function(string){
  #    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  #  }
  tab_expression_data_ages <- tab_expression_data_ages_tmp |>
    dplyr::select(age_id, age) |>
    dplyr::distinct() |>
    dplyr::mutate(age = tolower(age)) |>
    dplyr::mutate(age_min = 0) |>
    dplyr::mutate(age_max = 0) |>
    dplyr::mutate(age_median = 0) |>
    dplyr::mutate(age_mean = 0) |>
    dplyr::collect() |>
    dplyr::mutate(
      age_max = rowSums(tibble::tibble(
        a = as.numeric(stringr::str_detect(age, "first decade")) * 10,
        b = as.numeric(stringr::str_detect(age, "second decade")) * 20,
        c = as.numeric(stringr::str_detect(age, "third decade")) * 30,
        d = as.numeric(stringr::str_detect(age, "fourth decade")) * 40,
        e = as.numeric(stringr::str_detect(age, "fifth decade")) * 50,
        f = as.numeric(stringr::str_detect(age, "sixth decade")) * 60,
        g = as.numeric(stringr::str_detect(age, "seventh decade")) * 70,
        h = as.numeric(stringr::str_detect(age, "eighth decade")) * 80,
        i = as.numeric(stringr::str_detect(age, "ninth decade")) * 90,
        j = as.numeric(stringr::str_detect(age, "tenth decade")) * 100
      ), na.rm = TRUE),
      age_min = rowSums(tibble::tibble(
        a = readr::parse_number(stringr::str_extract(age, "^[\\d].*day")) / 365,
        b = readr::parse_number(stringr::str_extract(age, "^[\\d].*week")) / 52,
        c = readr::parse_number(stringr::str_extract(age, "^[\\d].*month")) / 12,
        d = readr::parse_number(stringr::str_extract(age, "^[\\d].*year"))
      ), na.rm = TRUE),
      age_max = dplyr::if_else(!stringr::str_detect(age, "decade"), age_min, age_max),
      age_min = dplyr::if_else(stringr::str_detect(age, "fertilization"), age_min * -1, age_min),
      age_max = dplyr::if_else(stringr::str_detect(age, "fertilization"), age_max * -1, age_max),
      age_mean = dplyr::if_else(stringr::str_detect(age, "fertilization"), age_mean * -1, age_mean),
      age_min = dplyr::if_else(stringr::str_detect(age, "decade"), age_max - 10, age_min),
      age_max = dplyr::if_else(stringr::str_detect(age, "decade"), age_max - 1, age_max),
      age_mean = dplyr::if_else(stringr::str_detect(age, "decade"), (age_max + age_min) / 2, age_mean)
    ) |>
    dplyr::mutate(age_median = age_mean) |>
    dplyr::compute(
      name = "tab_expression_data_ages",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_bases #
  print("Build tab_expression_data_bases and write to database")
  annotation_bgee <- BgeeDB::getAnnotation(bgee)
  DataAnnotation <- tibble::as_tibble(unique(annotation_bgee[[2]][, c("Experiment.ID", "Data.source.URL")]))
  DataAnnotation <- DataAnnotation |>
    dplyr::rename(data_base = Experiment.ID) |>
    dplyr::rename(url = Data.source.URL) |>
    dplyr::compute(
      name = "DataAnnotation",
      temporary = FALSE,
      overwrite = TRUE
    )

  tab_expression_data_bases_tmp <- rna_seq_selected |>
    dplyr::select(data_base) |>
    dplyr::distinct() |>
    dplyr::compute(
      name = "tab_expression_data_bases_tmp",
      temporary = FALSE,
      overwrite = TRUE
    )

  tab_expression_data_bases <-
    dplyr::left_join(tab_expression_data_bases_tmp, DataAnnotation, copy = TRUE) |>
    dplyr::compute(
      name = "tab_expression_data_bases",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_gender_properties #
  print("Build tab_expression_data_gender_properties and write to database")
  tab_expression_data_gender_properties <- rna_seq_selected |>
    dplyr::select(gender, age, tissue) |>
    dplyr::distinct() |>
    dplyr::mutate(developmental_stage_source = "-") |>
    tidyr::pivot_longer(!gender,
      names_to = "property",
      values_to = "property_value"
    ) |>
    dplyr::distinct() |>
    dplyr::compute(
      name = "tab_expression_data_gender_properties",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_genders #
  print("Build tab_expression_data_genders and write to database")
  tab_expression_data_genders <- rna_seq_selected |>
    dplyr::select(gender) |>
    dplyr::distinct() |>
    dplyr::mutate(information = paste0("A ", gender, " individual or population")) |>
    dplyr::compute(
      name = "tab_expression_data_genders",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_health_state #
  print("Build tab_expression_data_health_state and write to database")
  tab_expression_data_health_state <- rna_seq_selected |>
    dplyr::select(health_state) |>
    dplyr::distinct() |>
    dplyr::mutate(information = paste("This refers to from an ", health_state, " individual.", sep = "")) |>
    dplyr::compute(
      name = "tab_expression_data_health_state",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_health_state_properties #
  print("Build tab_expression_data_health_state_properties and write to database")
  tab_expression_data_health_state_properties <- rna_seq_selected |>
    dplyr::select(health_state) |>
    dplyr::distinct() |>
    dplyr::mutate(property = "health_state") |>
    dplyr::mutate(property_value = health_state) |>
    dplyr::compute(
      name = "tab_expression_data_health_state_properties",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_records #
  print("Build tab_expression_data_records and write to database")
  KEYS <- c(
    "data_source_id", "data_base_rec_id", "tissue", "health_state",
    "gender", "age_id"
  )
  assign("TODAY", as.character(Sys.Date()))
  tab_expression_data_records <-
    dplyr::left_join(
      tab_expression_data_records_tmp,
      rna_seq_selected |>
        dplyr::select(
          data_base, data_base_rec_id,
          tissue, health_state, gender, age
        ) |>
        dplyr::distinct()
    ) |>
    dplyr::left_join(dplyr::tbl(db_bgee_conn, "tab_expression_data_ages") |>
      dplyr::select(age_id, age)) |>
    dplyr::select(tidyselect::all_of(KEYS)) |>
    dplyr::distinct() |>
    dplyr::mutate(sample_source = "tissue") |>
    dplyr::mutate(data_base = "RNAseq") |>
    dplyr::mutate(last_refresh_date = TODAY) |>
    # mutate(last_refresh_date = as.character(as.POSIXlt(strptime(Sys.Date(), format = "%Y-%m-%d")))) |>
    # collect() |>
    dplyr::relocate(age_id, .after = tidyselect::last_col()) |>
    dplyr::compute(
      name = "tab_expression_data_records",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_sample_source_properties #
  print("Build tab_expression_data_sample_source_properties and write to database")
  tab_expression_data_sample_source_properties <- rna_seq_selected |>
    dplyr::select(tissue) |>
    dplyr::distinct() |>
    dplyr::rename(property_value = tissue) |>
    dplyr::mutate(property = "tissue_SOURCE") |>
    dplyr::mutate(sample_source = "tissue") |>
    dplyr::arrange(property_value) |>
    dplyr::compute(
      name = "tab_expression_data_sample_source_properties",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_sample_sources #
  print("Build tab_expression_data_sample_sources and write to database")
  tab_expression_data_sample_sources <-
    tibble::tibble(
      sample_source = c("CELL LINE", "PRIMARY CULTURE", "TISSUE", "UNSPECIFIED"),
      information = c(
        "The sample source is a cell line.",
        "The sample source is a tissue culture started from cells, tissues, or organs taken directly from the organism.",
        "The sample source is tissue.",
        "The sample source is unknown."
      )
    ) |>
    dplyr::compute(
      name = "tab_expression_data_sample_sources",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_tissue_properties #
  print("Build tab_expression_data_tissue_properties and write to database")
  tab_expression_data_tissue_properties <- rna_seq_selected |>
    dplyr::select(tissue) |>
    dplyr::distinct() |>
    dplyr::mutate(property_value = tissue) |>
    dplyr::mutate(property = "tissue") |>
    dplyr::arrange(property_value) |>
    dplyr::compute(
      name = "tab_expression_data_tissue_properties",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_tissues #
  print("Build tab_expression_data_tissues and write to database")
  tab_expression_data_tissues <- rna_seq_selected |>
    dplyr::select(tissue) |>
    dplyr::distinct() |>
    dplyr::mutate(information = "Organ") |>
    dplyr::arrange(tissue) |>
    dplyr::compute(
      name = "tab_expression_data_tissues",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_expression_data_units #
  print("Build tab_expression_data_units and write to database")
  # tab_expression_data_units <-
  #     tab_expression_data_values |>
  #     dplyr::select(unit) |>
  #     dplyr::distinct() |>
  #     dplyr::mutate(information = "NGS data")
  tab_expression_data_units <- tibble::tibble(
    unit = c("TPM", "FPKM", "RPKM", "Read.count"),
    information = c(
      "Transcript per million",
      "Fragments per kilobase million - Normalized to library size",
      "Reads per kilobase million",
      "Raw counts of transcripts"
    )
  ) |>
    dplyr::compute(
      name = "tab_expression_data_units",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_gene_names #
  print("Build tab_gene_names and write to database")
  # linking the unique variant ID to the different identifier
  KEYS <- c(
    "variant_id", "variant_name", "symbol", "official_full_name",
    "entrezid", "synonym", "protein_id", "preferred_name", "other_name"
  )
  if (SPECIE != "Human" & SPECIE != "Sheep") {
    KEYS <- c(KEYS, "homolog", "homolog_symbol")
  }

  tab_gene_variants <- dplyr::tbl(db_bgee_conn, "tab_gene_variants")
  tab_gene_names <-
    dplyr::left_join(tab_gene_variants, AnnotationTable) |>
    dplyr::select(tidyselect::any_of(KEYS)) |>
    dplyr::distinct() |>
    dplyr::rename(gene_id = entrezid) |>
    dplyr::mutate(gene_id = as.character(gene_id)) |>
    dplyr::distinct()

  tab_gene_names <- tab_gene_names |>
    tidyr::pivot_longer(!variant_id,
      names_to = "name_type",
      values_to = "gene_name"
    ) |>
    dplyr::filter(gene_name != "") |>
    dplyr::rename(gene_id = variant_id) |>
    dplyr::distinct() |>
    dplyr::arrange(gene_id) |>
    dplyr::compute(
      name = "tab_gene_names",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_gene_name_types #
  print("Build tab_gene_name_types and write to database")
  tab_gene_name_types <- tab_gene_names |>
    dplyr::select(name_type) |>
    dplyr::distinct() |>
    dplyr::mutate(information = paste0("The identifier '", name_type, "' is based on the Cran R biomaRt package.")) |>
    dplyr::compute(
      name = "tab_gene_name_types",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_gene_variants #

  ### tab_genes #
  print("Build tab_genes and write to database")
  tab_genes <- tab_gene_variants |>
    dplyr::select(gene_id) |>
    dplyr::distinct() |>
    dplyr::compute(
      name = "tab_genes",
      temporary = FALSE,
      overwrite = TRUE
    )

  ### tab_global_statistics #
  print("Build tab_global_statistics and write to database")
  # function in original Access DB == avg: 10^mean(Logarithmus([sample_count]/[total_count])/Logarithmus(10))
  tab_global_statistics <- tab_expression_data_values |>
    dplyr::select(unit, sample_count, total_count) |>
    dplyr::group_by(unit) |>
    dplyr::mutate(avg = dplyr::if_else(condition = unit %in% c("RPKM", "FPKM"),
      true = 10^mean(log(sample_count / total_count) / log(10), na.rm = TRUE),
      false = 10^mean(log(sample_count) / log(10), na.rm = TRUE)
    )) |>
    dplyr::select(unit, avg) |>
    dplyr::distinct() |>
    dplyr::compute(
      name = "tab_global_statistics",
      temporary = FALSE,
      overwrite = TRUE
    )

  ##### Write data into PK-Sim expression database #####
  print("Write data into PK-Sim expression database")
  assign(
    x = "tab_container_tissue",
    value = tibble::tibble(read.table("Code/tab_container_tissue.txt",
      header = 1, sep = "\t"
    ))
  )
  assign(
    x = "tab_dts_properties",
    value = tibble::tibble(read.table("Code/tab_dts_properties.txt",
      header = 1, sep = "\t"
    ))
  )
  colnames(tab_dts_properties) <- colnames(tab_dts_properties) |> tolower()
  Files <- ls(pattern = "^tab_.{4}")
  Files <- Files[-c(grep(x = Files, pattern = "_tmp"))]

  ## Connect to SQLite db
  # DB_Tables <- DBI::dbListTables(db_PKsim_conn)
  ## Remove old tables
  # for (TAB in DB_Tables) {
  #  base::ifelse(rlang::is_empty(base::grep("QRY", TAB)),
  #    DBI::dbRemoveTable(db_PKsim_conn, TAB),
  #    DBI::dbExecute(conn = db_PKsim_conn, paste0("DROP VIEW ", TAB, ";"))
  #  )
  # }
  ## Add empty tables SQLite data base.
  # CREATE_TABLE_COMMAND <- CREATE_TABLE$COMMAND
  # for (i in 1:length(CREATE_TABLE_COMMAND)) {
  #  DBI::dbExecute(conn = db_PKsim_conn, CREATE_TABLE_COMMAND[i])
  # }
  # Load modified tab files and write to db
  for (i in 1:length(c("tab_container_tissue", "tab_dts_properties"))) {
    TAB <- c("tab_container_tissue", "tab_dts_properties")[i]
    # info output to screen
    print(paste0("Write \'", tolower(TAB), "\' to \'", SPECIE, "\' PK-Sim expression database"))
    tmp_table_names <- TYPE[[which(TAB_s %in% tolower(TAB))]]
    names(tmp_table_names) <- colnames(TMP)

    DBI::dbWriteTable(
      conn = db_bgee_conn,
      name = TAB,
      value = data.frame(TMP),
      row.names = FALSE # field.types = tmp_table_names
    )
  }
  # }

  for (i in 1:length(VIEW)) { # if error occurs generate DB and show warning
    DBI::dbExecute(db_bgee_conn, VIEW[i])
  }

  for (i in 1:length(INDIZES)) { # if error occurs generate DB and show warning
    DBI::dbExecute(db_bgee_conn, INDIZES[i])
  }

  # CREATE_TABLE$ALTER_COMMAND
  # for (i in 1:length(CREATE_TABLE$ALTER_COMMAND)) { # if error occurs generate DB and show warning
  #   DBI::dbExecute(db_PKsim_conn, CREATE_TABLE$ALTER_COMMAND[i])
  # }
  print("Compress PK-Sim expression database")
  DBI::dbExecute(conn = db_bgee_conn, statement = "VACUUM") # clears pre-allocated disc space for db
  DBI::dbDisconnect(db_bgee_conn)
  gc() # clears used memory
}
