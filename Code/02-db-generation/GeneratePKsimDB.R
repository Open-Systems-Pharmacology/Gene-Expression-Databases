#' GeneratePKsimDB
#'
#' Description:
#'   This function generates a PK-Sim compatible gene expression database
#'   for a given species using Bgee RNA-Seq data. It supports filtering for
#'   ADME genes, different Bgee releases, and normalization options.
#'   The output is stored in the specified path and can be computed in RAM
#'   for performance.
#'
#' Parameters:
#'   SPECIE        - Character. Species name (e.g., "Dog", "Human").
#'   PATH          - Character. Output directory path. Default: current working directory.
#'   ADME_ONLY     - Logical. If TRUE, only ADME genes are included. Default: TRUE.
#'   RELEASE       - Character. Bgee release version (e.g., "15_2"). Default: "15_2".
#'   COMPUTE_IN_RAM- Logical. If TRUE, computation is performed in RAM. Default: FALSE.
#'   (INCLUDE_RPKM  - Logical. If TRUE, RPKM values are estimated. Default: FALSE.)
#'
#' Returns:
#'   No return value. Writes output files to disk.
#'
#' Example:
#'   GeneratePKsimDB(SPECIE = "Human", PATH = "./output", ADME_ONLY = TRUE)
#'
GeneratePKsimDB <- function(
    SPECIE = "Rabbit",
    PATH = here::here(),
    ADME_ONLY = TRUE,
    RELEASE = "15_2",
    COMPUTE_IN_RAM = TRUE) {
  ### Dependencies ####
  switch(RELEASE,
    "13_2" = {
      stop("Code was discontinued for Bgee release 13_2, try 15_2")
    },
    "14_0" = {
      stop("Code was discontinued for Bgee release 14_0, try 15_2")
    },
    "14_1" = {
      stop("Code was discontinued for Bgee release 14_1, try 15_2")
    },
    "14_2" = {
      stop("Code was discontinued for Bgee release 14_2, try 15_2")
    },
    "15_0" = {
      message("Code was not tested for Bgee release 15_0")
    },
    "15_1" = {
      message("Code was not tested for Bgee release 15_1")
    },
    "15_2" = {
      message("Code was tested for Bgee release 15_2")
    },
    {
      message(
        paste0(
          "Input Bgee release: ", RELEASE,
          " was not recognized. Default: 15_2 is used instead"
        )
      )
      RELEASE <- "15_2"
    }
  )

  # Test if input species is valid
  source(paste0(PATH, "/Code/03-helpers/helper_Species.R"))

  if (!SPECIE %in% ALL_SPECIE) {
    stop(paste0(
      "Given SPECIE: '", SPECIE, "' not in list. \n Must be one of: ",
      paste(ALL_SPECIE, collapse = " , ")
    ))
  }

  # is needed to allow download of human data (65 GB takes some time)
  old_timeout <- getOption("timeout")
  options(timeout = 60 * 60 * 60)
  on.exit(options(timeout = old_timeout), add = TRUE)

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
      DATASET <- "mnemestrina_gene_ensembl"
    },
    Minipig = {
      SPECIE_LAT <- "Sus_scrofa"
      DATASET <- "sscrofa_gene_ensembl"
    },
    Dog = {
      SPECIE_LAT <- "Canis_lupus_familiaris"
      DATASET <- "clfamiliaris_gene_ensembl"
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
    Guineapig = {
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

  print(getwd())
  dir.create(file.path(PATH, "BgeeDBs"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PATH, "PK-Sim DBs"), recursive = TRUE, showWarnings = FALSE)
  old_wd <- setwd("BgeeDBs/")
  on.exit(setwd(old_wd), add = TRUE)
  print(paste0("Fetch Bgee expression data sets for ", SPECIE))

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
  DB_PKsim <- paste0(
    DB_PKsim, "_BgeeRelease_",
    BegeeRelease, ".expressionDB" # "_", Sys.Date(), ".expressionDB"
  )
  PATH2PKsim <- paste0(PATH, "/PK-Sim DBs/", SPECIE, "/")
  PATH_DB_Bgee <- bgee[["pathToData"]]
  dir.create(PATH2PKsim, recursive = TRUE, showWarnings = FALSE)

  #### Connections to local data bases ####
  print(paste0("Connect to local ", SPECIE, " bgee data base"))
  # local bgee DB
  db_bgee_conn <-
    DBI::dbConnect(RSQLite::SQLite(),
      paste0(PATH_DB_Bgee, "/", gsub(x = DB_Bgee, pattern = " ", replacement = "_")),
      synchronous = "off",
      cache_size = -1000
    )
  on.exit(try(DBI::dbDisconnect(db_bgee_conn), silent = TRUE), add = TRUE)
  DB_Tables <- DBI::dbListTables(conn = db_bgee_conn)

  # local biomart DB
  print(paste0("Connect to local ", SPECIE, " BioMarts"))
  db_biomart_conn <-
    DBI::dbConnect(RSQLite::SQLite(),
      paste0(PATH, "/BioMarts/All_Species_BioMarts.DB"),
      synchronous = "off",
      cache_size = -1000
    )
  on.exit(try(DBI::dbDisconnect(db_biomart_conn), silent = TRUE), add = TRUE)

  # local pkSim DB
  print(paste0("Connect to local ", SPECIE, " PK-Sim DB"))
  db_PKsim_conn <- DBI::dbConnect(RSQLite::SQLite(),
    paste0(PATH2PKsim, DB_PKsim),
    synchronous = "off",
    cache_size = -1000
  )
  on.exit(try(DBI::dbDisconnect(db_PKsim_conn), silent = TRUE), add = TRUE)

  # If tables are empty they have been created but not filled with data
  if (rlang::is_empty(DB_Tables)) {
    print(paste0("Build local ", SPECIE, " BgeeDB data bases"))
    # If data is not yet locally stored execute getSampleProcessedData()
    # to initiate download, currently the latest version
    # getSampleProcessedData() switched to a forced load of whole
    # db content into ram. This can cause to errors when the function
    # GeneratePKsimDB.R is called the first time, re-executing the function
    # usually resolves the issue (then DB is already downloaded).
    print(paste0("Load ", SPECIE, " expression data into RAM"))
    rna_seq_selected <- BgeeDB::getSampleProcessedData(bgee)
    DBI::dbExecute(db_bgee_conn, 'UPDATE rna_seq SET "Anatomical.entity.name" = REPLACE("Anatomical.entity.name", \'"\', \'\')')
    DBI::dbExecute(db_bgee_conn, 'UPDATE rna_seq SET "Stage.name" = REPLACE("Stage.name", \'"\', \'\')')
    DBI::dbExecute(db_bgee_conn, 'UPDATE rna_seq SET Strain = REPLACE(Strain, \'"\', \'\')')
    already_loaded <- TRUE
  } else {
    already_loaded <- FALSE
  }

  print("Re-organize data from BgeeDB for PK-Sim compatibility")
  if (COMPUTE_IN_RAM) {
    # as local table in RAM
    if (!already_loaded) {
      print("Fetch expression data from local BgeeDB.")
      rna_seq_selected <- BgeeDB::getSampleProcessedData(bgee)
    }
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
      dplyr::mutate(Sex = "UNSPECIFIED") |>
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
  )

  # Further format, rename and prepare data ####
  rna_seq_selected <- rna_seq_selected |>
    dplyr::mutate(gender = toupper(gender)) |>
    dplyr::mutate(age = toupper(age)) |>
    dplyr::mutate(strain = toupper(strain)) |>
    dplyr::mutate(gender = dplyr::if_else(gender == toupper("NA"),
      toupper("UNSPECIFIED"), gender
    )) |>
    dplyr::mutate(strain = dplyr::if_else(strain == toupper("NA"),
      toupper("UNSPECIFIED"), strain
    )) |>
    dplyr::mutate(tissue = dplyr::if_else(tissue == toupper("NA"),
      toupper("UNSPECIFIED"), tissue
    )) |>
    dplyr::mutate(tissue = toupper(tissue)) |>
    dplyr::mutate(state = toupper("NORMAL")) |>
    dplyr::mutate(health_state = toupper(paste(strain, state, sep = " ")))

  if (COMPUTE_IN_RAM) {
    rna_seq_selected <- rna_seq_selected |>
      dplyr::compute(
        name = "rna_seq_selected",
        temporary = FALSE,
        overwrite = TRUE
      )
  }

  ##### Add Gene IDs / Annotations needed for allow broad reach ####
  print("Annotete expression data with gene & protein identifiers.")
  if (ADME_ONLY) {
    AnnotationTable <-
      DBI::dbReadTable(db_biomart_conn, paste0(SPECIE, "_ADME")) |>
      dplyr::select(-c("start_position", "end_position")) |>
      dplyr::mutate(OFFICIAL_FULL_NAME = trimws(OFFICIAL_FULL_NAME)) |>
      readr::type_convert()

    colnames(AnnotationTable) <- tolower(colnames(AnnotationTable))

    DBI::dbWriteTable(
      conn = db_bgee_conn,
      name = "AnnotationTable_ADME",
      value = data.frame(AnnotationTable),
      row.names = FALSE,
      overwrite = TRUE
    )
  } else {
    AnnotationTable <-
      DBI::dbReadTable(db_biomart_conn, paste0(SPECIE, "_Annotations")) |>
      dplyr::select(-c("start_position", "end_position")) |>
      dplyr::mutate(OFFICIAL_FULL_NAME = trimws(OFFICIAL_FULL_NAME)) |>
      readr::type_convert()

    colnames(AnnotationTable) <- tolower(colnames(AnnotationTable))

    DBI::dbWriteTable(
      conn = db_bgee_conn,
      name = "AnnotationTable",
      value = data.frame(AnnotationTable),
      row.names = FALSE,
      overwrite = TRUE
    )
  }
  # For only homologous genes, add "HOMOLOG" prefix in Gene Symbol;
  # This is necessary since PK-SIM gene names are based on Symbols

  #### generate TAB tables from data base data frame ####
  print("Write data to TAB files")
  # Currently realizing the query in local RAM (collect) needed for
  # copying accross database is the ratelimiting step.
  # Try sqlpipe for data transfer between local DBs without detour through local RAM:
  # https://github.com/sqlpipe/sqlpipe?tab=readme-ov-file
  # icacls sqlpipe /grant USER:RX
  if (COMPUTE_IN_RAM) {
    AnnotationTable <- dplyr::collect(AnnotationTable)
  } else {
    table_name <- if (ADME_ONLY) "AnnotationTable_ADME" else "AnnotationTable"
    AnnotationTable <- dplyr::tbl(db_bgee_conn, table_name)
  }

  # Combine expression data and annotation information
  rna_seq_selected <- dplyr::left_join(
    AnnotationTable |>
      dplyr::select(variant_name) |>
      dplyr::distinct(),
    by = dplyr::join_by(variant_name),
    rna_seq_selected
  )

  # Remove rna_seq_selected entries where age is NA
  rna_seq_selected <- rna_seq_selected |>
    dplyr::filter(!is.na(age)) |>
    dplyr::mutate(age = as.character(age))

  # Add gene IDs
  print("Build tab_gene_variants")
  tab_gene_variants <- rna_seq_selected |>
    dplyr::select(variant_name) |>
    dplyr::distinct() |>
    dplyr::collect() |>
    dplyr::mutate(variant_id = dplyr::row_number()) |>
    dplyr::mutate(gene_id = variant_id) |>
    dplyr::relocate(gene_id) |>
    dplyr::arrange(variant_name)

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_gene_variants",
    value = data.frame(tab_gene_variants),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_gene_variants"][[1]]
  )
  # Add data source IDs
  print("Build tab_expression_data_records_tmp")
  tab_expression_data_records_tmp <- rna_seq_selected |>
    dplyr::select(data_base_rec_id) |>
    dplyr::distinct() |>
    dplyr::collect() |>
    dplyr::mutate(data_source_id = dplyr::row_number()) #|>
  # dplyr::mutate(AGE = tolower(age))

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_records_tmp",
    value = data.frame(tab_expression_data_records_tmp),
    overwrite = TRUE,
    field.types = c(data_base_rec_id = "text", data_source_id = "bigint")
  )

  # Add data age IDs
  print("Build tab_expression_data_ages_tmp")
  tab_expression_data_ages_tmp <- rna_seq_selected |>
    dplyr::select(age) |>
    dplyr::distinct() |>
    dplyr::collect() |>
    dplyr::mutate(age = toupper(age)) |>
    dplyr::mutate(age_id = dplyr::row_number()) |>
    dplyr::arrange(age_id, age) |>
    dplyr::relocate(age_id)

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_ages_tmp",
    value = data.frame(tab_expression_data_ages_tmp),
    overwrite = TRUE,
    field.types = c(age_id = "bigint", age = "text")
  )

  #### Merge different expression measures (later Units in PK-Sim DB) ####
  # Here cut offs could be considered
  print("Merge expression data from multiple units.")
  tpm_table <- rna_seq_selected |>
    dplyr::select(!tidyselect::any_of(c("FPKM", "Read.count", "RPKM"))) |>
    dplyr::group_by(data_base_rec_id) |>
    dplyr::summarise(total_count = sum(TPM, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::collect() |>
    dplyr::mutate(total_count = total_count / 10^6)

  print("Build tpm_table")
  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tpm_table",
    value = data.frame(tpm_table),
    overwrite = TRUE,
    field.types = c(data_base_rec_id = "text", total_count = "real")
  )

  ### tab_expression_data_properties #
  print("Build tab_expression_data_properties")
  KEYS <- c(
    "data_source_id", "data_base", "tissue", "health_state",
    "gender", "age"
  )
  TMP <- dplyr::tbl(db_bgee_conn, "tab_expression_data_records_tmp")
  tab_expression_data_properties_tmp <-
    dplyr::left_join(rna_seq_selected, TMP, copy = TRUE) |>
    dplyr::select(tidyselect::all_of(KEYS)) |>
    dplyr::distinct()

  tab_expression_data_properties <-
    tab_expression_data_properties_tmp |>
    tidyr::pivot_longer(!data_source_id,
      names_to = "property",
      values_to = "property_value"
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(property = toupper(property))

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_properties",
    value = as.data.frame(tab_expression_data_properties),
    overwrite = TRUE,
    field.types = c(data_source_id = "bigint", property = "text", property_value = "text")
  )
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_properties SET property_value = REPLACE(property_value, \'"\', \'\')')

  if (COMPUTE_IN_RAM) {
    tab_expression_data_properties <- DBI::dbReadTable(db_bgee_conn, "tab_expression_data_properties")
  } else {
    tab_expression_data_properties <- dplyr::tbl(db_bgee_conn, "tab_expression_data_properties")
  }
  ### tab_expression_data_values #
  print("Build tab_expression_data_values")
  KEYS <- c("variant_id", "data_source_id", "sample_count", "total_count", "unit")
  tab_expression_data_values <- dplyr::left_join(rna_seq_selected,
    dplyr::tbl(db_bgee_conn, "tpm_table"),
    copy = TRUE
  ) |>
    dplyr::mutate(unit = "TPM") |>
    dplyr::rename(sample_count = TPM) |>
    dplyr::left_join(dplyr::tbl(db_bgee_conn, "tab_gene_variants"), copy = TRUE) |>
    dplyr::left_join(dplyr::tbl(db_bgee_conn, "tab_expression_data_records_tmp"), copy = TRUE) |>
    dplyr::select(tidyselect::all_of(KEYS)) |>
    dplyr::arrange(variant_id, data_source_id) |>
    dplyr::distinct() |>
    dplyr::collapse()

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_values",
    value = as.data.frame(tab_expression_data_values),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_values"][[1]]
  )

  ### tab_expression_data_age_properties #
  print("Build tab_expression_data_age_properties")
  # optional age is extracted form the data sets, in bgee strings like "2-3 weeks old rat" are given
  tab_expression_data_age_properties <- tab_expression_data_ages_tmp |>
    dplyr::select(age_id, age) |>
    dplyr::distinct() |>
    tidyr::pivot_longer(!age_id,
      names_to = "property",
      values_to = "property_value"
    ) |>
    dplyr::mutate(property = toupper(property))

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_age_properties",
    value = as.data.frame(tab_expression_data_age_properties),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_age_properties"][[1]]
  )

  ### tab_expression_data_ages #
  print("Build tab_expression_data_ages")
  # for mice also Theiler Stage could be included: https://www.emouseatlas.org/emap/ema/theiler_stages/StageDefinition/stagedefinition.html#dpc
  # for human also Carnegie Stages could be included: https://embryology.ch/de/embryogenese/periode-embryonnaire/carnegie-stadien/
  tab_expression_data_ages <- tab_expression_data_ages_tmp |>
    dplyr::select(age_id, age) |>
    dplyr::distinct() |>
    dplyr::collect() |>
    dplyr::mutate(age = toupper(age)) |>
    dplyr::mutate(age_min = 0) |>
    dplyr::mutate(age_max = 0) |>
    dplyr::mutate(age_max = base::rowSums(tibble::tibble(
      a = matrix(stringr::str_detect(string = age, pattern = toupper("FIRST DECADE")) * 10, ncol = 1),
      b = matrix(stringr::str_detect(string = age, pattern = toupper("SECOND DECADE")) * 20, ncol = 1),
      c = matrix(stringr::str_detect(string = age, pattern = toupper("THIRD DECADE")) * 30, ncol = 1),
      d = matrix(stringr::str_detect(string = age, pattern = toupper("FOURTH DECADE")) * 40, ncol = 1),
      e = matrix(stringr::str_detect(string = age, pattern = toupper("FIFTH DECADE")) * 50, ncol = 1),
      f = matrix(stringr::str_detect(string = age, pattern = toupper("SIXTH DECADE")) * 60, ncol = 1),
      g = matrix(stringr::str_detect(string = age, pattern = toupper("SEVENTH DECADE")) * 70, ncol = 1),
      h = matrix(stringr::str_detect(string = age, pattern = toupper("EIGHTH DECADE")) * 80, ncol = 1),
      i = matrix(stringr::str_detect(string = age, pattern = toupper("NINTH DECADE")) * 90, ncol = 1),
      j = matrix(stringr::str_detect(string = age, pattern = toupper("TENTH DECADE")) * 100, ncol = 1)
    ), na.rm = TRUE)) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "FIRST MONTH", replacement = "1")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "SECOND MONTH", replacement = "2")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "THIRD MONTH", replacement = "3")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "FOURTH MONTH", replacement = "4")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "FIFTH MONTH", replacement = "5")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "SIXTH MONTH", replacement = "6")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "SEVENTH MONTH", replacement = "7")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "EIGHTH MONTH", replacement = "8")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "NINTH MONTH", replacement = "9")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "DAY 12", replacement = "12 DAY")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "DAY 14", replacement = "14 DAY")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "DAY 16", replacement = "16 DAY")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "DAY 20", replacement = "20 DAY")) |>
    dplyr::mutate(age = stringr::str_replace(string = age, pattern = "DAY 21", replacement = "21 DAY"))

  tab_expression_data_ages <- tab_expression_data_ages |>
    dplyr::mutate(age_min = base::rowSums(tibble::tibble(
      a = matrix(readr::parse_number(stringr::str_extract(string = age, pattern = ("^[\\d].*DAY"))) / 365, ncol = 1),
      b = matrix(readr::parse_number(stringr::str_extract(string = age, pattern = ("^[\\d].*WEEK"))) / 52, ncol = 1),
      c = matrix(readr::parse_number(stringr::str_extract(string = age, pattern = ("^[\\d].*MONTH"))) / 12, ncol = 1),
      d = matrix(readr::parse_number(stringr::str_extract(string = age, pattern = ("^[\\d].*YEAR"))), ncol = 1)
    ), na.rm = TRUE))
  tab_expression_data_ages <- tab_expression_data_ages |>
    dplyr::mutate(age_median = 0) |>
    dplyr::mutate(age_mean = 0) |>
    dplyr::mutate(age_max = dplyr::if_else(!(stringr::str_detect(string = age, pattern = toupper("DECADE"))), age_min, age_max)) |>
    dplyr::mutate(age_min = dplyr::if_else(stringr::str_detect(string = age, pattern = toupper("FERTILIZATION")), -9 / 12 + age_min, age_min)) |>
    dplyr::mutate(age_max = dplyr::if_else(stringr::str_detect(string = age, pattern = toupper("FERTILIZATION")), -9 / 12 + age_max, age_max)) |>
    dplyr::mutate(age_min = dplyr::if_else(stringr::str_detect(string = age, pattern = toupper("GESTATION")), -9 / 12 + age_min, age_min)) |>
    dplyr::mutate(age_max = dplyr::if_else(stringr::str_detect(string = age, pattern = toupper("GESTATION")), -9 / 12 + age_max, age_max)) |>
    dplyr::mutate(age_min = dplyr::if_else(stringr::str_detect(string = age, pattern = toupper("EMBRYO")), -9 / 12 + age_min, age_min)) |>
    dplyr::mutate(age_max = dplyr::if_else(stringr::str_detect(string = age, pattern = toupper("EMBRYO")), -9 / 12 + age_max, age_max)) |>
    dplyr::mutate(age_min = dplyr::if_else(stringr::str_detect(string = age, pattern = toupper("DECADE")), age_max - 10, age_min)) |>
    dplyr::mutate(age_max = dplyr::if_else(stringr::str_detect(string = age, pattern = toupper("DECADE")), age_max - 1, age_max))

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_ages",
    value = as.data.frame(tab_expression_data_ages),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_ages"][[1]]
  )

  ### tab_expression_data_bases #
  print("Build tab_expression_data_bases")
  annotation_bgee <- BgeeDB::getAnnotation(bgee)
  DataAnnotation <- tibble::as_tibble(unique(annotation_bgee[[2]][, c("Experiment.ID", "Data.source.URL")]))
  DataAnnotation <- DataAnnotation |>
    dplyr::rename(data_base = Experiment.ID) |>
    dplyr::rename(url = Data.source.URL) # |> rename(NAME = Experiment.name)

  tab_expression_data_bases <- rna_seq_selected |>
    dplyr::select(data_base) |>
    dplyr::distinct()

  if (!COMPUTE_IN_RAM) {
    DBI::dbWriteTable(
      conn = db_bgee_conn,
      name = "DataAnnotation",
      value = data.frame(DataAnnotation),
      overwrite = TRUE
    )
    DataAnnotation <- dplyr::tbl(db_bgee_conn, "DataAnnotation")
  }
  tab_expression_data_bases <-
    dplyr::left_join(tab_expression_data_bases, DataAnnotation, copy = TRUE)

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_bases",
    value = as.data.frame(tab_expression_data_bases),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_bases"][[1]]
  )

  ### tab_expression_data_gender_properties #
  print("Build tab_expression_data_gender_properties")
  tab_expression_data_gender_properties <- rna_seq_selected |>
    dplyr::select(gender, age, tissue) |>
    dplyr::distinct() |>
    dplyr::mutate(developmental_stage_source = "-") |>
    tidyr::pivot_longer(!gender,
      names_to = "property",
      values_to = "property_value"
    ) |>
    dplyr::mutate(property = toupper(property)) |>
    dplyr::distinct()

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_gender_properties",
    value = as.data.frame(tab_expression_data_gender_properties),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_gender_properties"][[1]]
  )
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_gender_properties SET property_value = REPLACE(property_value, \'"\', \'\')')

  ### tab_expression_data_genders #
  print("Build tab_expression_data_genders")
  tab_expression_data_genders <- rna_seq_selected |>
    dplyr::select(gender) |>
    dplyr::distinct() |>
    dplyr::mutate(information = paste0("A ", gender, " individual or population"))

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_genders",
    value = as.data.frame(tab_expression_data_genders),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_genders"][[1]]
  )

  ### tab_expression_data_health_state #
  print("Build tab_expression_data_health_state")
  tab_expression_data_health_state <- rna_seq_selected |>
    dplyr::select(health_state) |>
    dplyr::distinct() |>
    dplyr::mutate(information = paste("This refers to data from '", health_state, "' individual.", sep = ""))

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_health_state",
    value = as.data.frame(tab_expression_data_health_state),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_health_state"][[1]]
  )
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_health_state SET health_state = REPLACE(health_state, \'"\', \'\')')
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_health_state SET information = REPLACE(information, \'"\', \'\')')

  ### tab_expression_data_health_state_properties #
  print("Build tab_expression_data_health_state_properties")
  tab_expression_data_health_state_properties <- rna_seq_selected |>
    dplyr::select(health_state) |>
    dplyr::distinct() |>
    dplyr::mutate(property = "HEALTH_STATE") |>
    dplyr::mutate(property_value = health_state)

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_health_state_properties",
    value = as.data.frame(tab_expression_data_health_state_properties),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_health_state_properties"][[1]]
  )
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_health_state_properties SET health_state = REPLACE(health_state, \'"\', \'\')')
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_health_state_properties SET property_value = REPLACE(property_value, \'"\', \'\')')

  ### tab_expression_data_records #
  print("Build tab_expression_data_records")
  KEYS <- c(
    "data_source_id", "data_base_rec_id", "tissue", "health_state",
    "gender", "age_id"
  )
  assign("TODAY", as.character(Sys.Date()))
  tab_expression_data_records <-
    dplyr::left_join(tab_expression_data_records_tmp,
      rna_seq_selected |>
        dplyr::select(
          data_base, data_base_rec_id,
          tissue, health_state, gender, age
        ) |>
        dplyr::distinct(),
      copy = TRUE
    ) |>
    dplyr::left_join(
      dplyr::tbl(db_bgee_conn, "tab_expression_data_ages") |>
        dplyr::select(age_id, age),
      copy = TRUE
    ) |>
    dplyr::select(tidyselect::all_of(KEYS)) |>
    dplyr::distinct() |>
    dplyr::mutate(sample_source = "tissue") |>
    dplyr::mutate(data_base = "RNAseq") |>
    dplyr::mutate(last_refresh_date = TODAY) |>
    # mutate(last_refresh_date = as.character(as.POSIXlt(strptime(Sys.Date(), format = "%Y-%m-%d")))) |>
    # collect() |>
    dplyr::relocate(age_id, .after = tidyselect::last_col())

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_records",
    value = as.data.frame(tab_expression_data_records),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_records"][[1]]
  )
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_records SET tissue = REPLACE(tissue, \'"\', \'\')')
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_records SET health_state = REPLACE(health_state, \'"\', \'\')')

  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_age_properties SET property_value = REPLACE(property_value, \'"\', \'\')')
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_ages_tmp SET age = REPLACE(age, \'"\', \'\')')
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_ages SET age = REPLACE(age, \'"\', \'\')')

  ### tab_expression_data_sample_source_properties #
  print("Build tab_expression_data_sample_source_properties")
  tab_expression_data_sample_source_properties <- rna_seq_selected |>
    dplyr::select(tissue) |>
    dplyr::distinct() |>
    dplyr::rename(property_value = tissue) |>
    dplyr::mutate(property = "TISSUE_SOURCE") |>
    dplyr::mutate(sample_source = "TISSUE") |>
    dplyr::arrange(property_value)

  DBI::dbWriteTable(
    conn = db_bgee_conn, name = "tab_expression_data_sample_source_properties",
    value = as.data.frame(tab_expression_data_sample_source_properties),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_sample_source_properties"][[1]]
  )
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_sample_source_properties SET property_value = REPLACE(property_value, \'"\', \'\')')

  ### tab_expression_data_sample_sources #
  print("Build tab_expression_data_sample_sources")
  tab_expression_data_sample_sources <-
    tibble::tibble(
      sample_source = c("CELL LINE", "PRIMARY CULTURE", "TISSUE", "UNSPECIFIED"),
      information = c(
        "The sample source is a cell line.",
        "The sample source is a tissue culture started from cells, tissues, or organs taken directly from the organism.",
        "The sample source is tissue.",
        "The sample source is unknown."
      )
    )

  DBI::dbWriteTable(
    conn = db_bgee_conn, name = "tab_expression_data_sample_sources",
    value = as.data.frame(tab_expression_data_sample_sources),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_sample_sources"][[1]]
  )

  ### tab_expression_data_tissue_properties #
  print("Build tab_expression_data_tissue_properties")
  tab_expression_data_tissue_properties <- rna_seq_selected |>
    dplyr::select(tissue) |>
    dplyr::distinct() |>
    dplyr::mutate(property_value = tissue) |>
    dplyr::mutate(property = "TISSUE") |>
    dplyr::arrange(property_value)

  DBI::dbWriteTable(
    conn = db_bgee_conn, name = "tab_expression_data_tissue_properties",
    value = as.data.frame(tab_expression_data_tissue_properties),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_tissue_properties"][[1]]
  )
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_tissue_properties SET property_value = REPLACE(property_value, \'"\', \'\')')
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_tissue_properties SET tissue = REPLACE(tissue, \'"\', \'\')')

  ### tab_expression_data_tissues #
  print("Build tab_expression_data_tissues")
  tab_expression_data_tissues <- rna_seq_selected |>
    dplyr::select(tissue) |>
    dplyr::distinct() |>
    dplyr::mutate(information = "Organ") |>
    dplyr::arrange(tissue)

  DBI::dbWriteTable(
    conn = db_bgee_conn, name = "tab_expression_data_tissues",
    value = as.data.frame(tab_expression_data_tissues),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_tissues"][[1]]
  )
  DBI::dbExecute(db_bgee_conn, 'UPDATE tab_expression_data_tissues SET tissue = REPLACE(tissue, \'"\', \'\')')

  ### tab_expression_data_units #
  print("Build tab_expression_data_units")
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
  )

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_expression_data_units",
    value = as.data.frame(tab_expression_data_units),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_expression_data_units"][[1]]
  )

  ### tab_gene_names #
  print("Build tab_gene_names")
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
    dplyr::left_join(tab_gene_variants, AnnotationTable, copy = TRUE) |>
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
    dplyr::mutate(name_type = toupper(name_type)) |>
    dplyr::arrange(gene_id)

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_gene_names",
    value = as.data.frame(tab_gene_names),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_gene_names"][[1]]
  )

  ### tab_gene_name_types #
  print("Build tab_gene_name_types")
  tab_gene_name_types <- tab_gene_names |>
    dplyr::select(name_type) |>
    dplyr::distinct() |>
    dplyr::mutate(information = paste0("The identifier '", name_type, "' is based on the Cran R biomaRt package.")) #|>
  # dplyr::compute()

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_gene_name_types",
    value = as.data.frame(tab_gene_name_types),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_gene_name_types"][[1]]
  )

  ### tab_gene_variants #

  ### tab_genes #
  print("Build tab_genes")
  tab_genes <- tab_gene_variants |>
    dplyr::select(gene_id) |>
    dplyr::distinct()

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_genes",
    value = as.data.frame(tab_genes),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_genes"][[1]]
  )

  ### tab_global_statistics #
  print("Build tab_global_statistics")
  # function in original Access DB == avg: 10^mean(Logarithmus([sample_count]/[total_count])/Logarithmus(10))
  tab_global_statistics <- tab_expression_data_values |>
    dplyr::select(unit, sample_count, total_count) |>
    dplyr::group_by(unit) |>
    dplyr::mutate(avg = dplyr::if_else(condition = unit %in% c("RPKM", "FPKM"),
      true = 10^mean(log(sample_count / total_count) / log(10), na.rm = TRUE),
      false = 10^mean(log(sample_count) / log(10), na.rm = TRUE)
      # false = geometric_mean(sample_count, na.rm = TRUE)
    )) |>
    dplyr::select(unit, avg) |>
    dplyr::distinct()

  DBI::dbWriteTable(
    conn = db_bgee_conn,
    name = "tab_global_statistics",
    value = as.data.frame(tab_global_statistics),
    overwrite = TRUE,
    field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == "tab_global_statistics"][[1]]
  )

  ##### Write data into PK-Sim expression database #####
  print(paste0("Write '", SPECIE, "' data into PK-Sim expression database"))
  assign(
    x = "tab_container_tissue",
    value = tibble::tibble(read.table(paste0(PATH, "/Code/03-helpers/tab_container_tissue.txt"),
      header = 1, sep = "\t", quote = '"'
    ))
  )
  assign(
    x = "tab_dts_properties",
    value = tibble::tibble(read.table(paste0(PATH, "/Code/03-helpers/tab_dts_properties.txt"),
      header = 1, sep = "\t"
    ))
  )
  #  tab_dts_properties <- tab_dts_properties |>
  #  dplyr::mutate(table_name = tolower(table_name),
  #                column_name = tolower(column_name))
  # Files <- ls(pattern = "^tab_.{4}")
  # Files <- Files[-c(grep(x = Files, pattern = "_tmp"))]
  # Files <- Files[-c(grep(x = Files, pattern = "_TMP"))]
  # Files <- DBI::dbListTables(db_bgee_conn)
  # Files <- Files[grep(x = Files, pattern = "^TAB_.{4}")]

  # Connect to SQLite db
  DB_Tables <- DBI::dbListTables(db_PKsim_conn)
  DB_Object <- DBI::dbListObjects(db_PKsim_conn)

  # Remove old views
  for (view in names(VIEW_TABLE)) {
    # DBI::dbExecute(conn = db_PKsim_conn, paste0("DROP VIEW ", toupper(view), ";"))
    DBI::dbExecute(conn = db_PKsim_conn, paste0("DROP VIEW IF EXISTS ", view, ";"))
  }
  # Remove index
  for (tab in CREATE_TABLE$TAB) {
    DBI::dbExecute(conn = db_PKsim_conn, paste0("DROP INDEX IF EXISTS ", tab, ";"))
  }
  # Remove tables
  for (tab in CREATE_TABLE$TAB) {
    DBI::dbExecute(conn = db_PKsim_conn, paste0("DROP TABLE IF EXISTS ", tab, ";"))
  }

  # Add empty tables SQLite data base.
  for (command in CREATE_TABLE$COMMAND) {
    DBI::dbExecute(conn = db_PKsim_conn, command)
  }

  # Load modified tab files and write to db
  for (tab in CREATE_TABLE$TAB) {
    # info output to screen
    print(paste0("Write \'", tolower(tab), "\' to \'", SPECIE, "\' PK-Sim expression database"))

    if (sum(tab == c("tab_container_tissue", "tab_dts_properties")) == 1) {
      DBI::dbWriteTable(
        conn = db_PKsim_conn,
        name = tolower(tab),
        value = get(tab),
        row.names = FALSE,
        append = TRUE # ,
        # overwrite = FALSE,
        # field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == tab][[1]]
      )
    } else {
      TMP <- dplyr::tbl(db_bgee_conn, tab)
      DBI::dbWriteTable(
        conn = db_PKsim_conn,
        name = tolower(tab),
        value = data.frame(TMP) |> dplyr::collect(),
        row.names = FALSE,
        append = TRUE # ,
        # overwrite = FALSE,
        # field.types = CREATE_TABLE$TYPE[CREATE_TABLE$TAB == tab][[1]]
      )
    }
  }

  print(paste0("Create views for '", SPECIE, "\' PK-Sim expression database"))
  for (view in names(VIEW_TABLE)) {
    DBI::dbExecute(db_PKsim_conn, VIEW_TABLE[[view]])
  }

  print(paste0("Set table indizes for '", SPECIE, "\' PK-Sim expression database"))
  for (i in seq_along(INDIZES)) { # if error occurs generate DB and show warning
    DBI::dbExecute(db_PKsim_conn, INDIZES[i])
  }

  print(paste0("Compress '", SPECIE, "' PK-Sim expression database"))
  DBI::dbExecute(conn = db_PKsim_conn, statement = "VACUUM") # clears pre-allocated disc space for db
  gc() # clears used memory
  # DB connections and working directory are restored by on.exit() handlers
}
