GeneratePKsimDB <- function(SPECIE = "Rabbit", 
                            PATH = getwd(),
                            ADME_ONLY = F,
                            RELEASE = "15_0", #RELEASE = "14_0" # "14_2" # "15_0"
                            COMPUTE_IN_RAM = F, #CUT_OFF = c(Read.count = 10, RPKM = 1, FPKM = 1, TPM = 1)
                            TPM_ONLY = T, # use only TPM normalized expression data (reduces workload and memory usage)
                            INCLUDE_RPKM = F # indicate if RPKM values should be estimated (gene length needed)
) {
  ### Dependencies ####
  # Test if input species is valid
  ALL_SPECIE <- c("Human","Monkey","Minipig","Dog","Mouse","Rat","Rabbit","Zebrafish","Cattle","Horse","Cat","GuineaPig","Chicken","Goat","Sheep","Turkey","Monkey_CrabEatingMacaque")
  
  if (!SPECIE %in% ALL_SPECIE) {
    stop(paste0("Given SPECIE: '", SPECIE, "' not in list. \n Must be one of: ", paste(ALL_SPECIE, collapse = " , ") ))
  }
  
  # attach needed packages; 
  require("BgeeDB") 
  require("biomaRt") 
  require("RSQLite") 
  require("dplyr") 
  require("readr") 
  # require("reshape2") 
  # require("data.table") 
  # require("tidyverse") 
  # require("biomartr") 
  # require("dm") 
  
  source(paste0(PATH, "/Code/helper_SQL_Commands.R")) # load SQL commands to set views and column types
  
  options(timeout = 60*60*3) # is needed to allow download of human data (65 GB take longer than 1 min ;-) )
  
  #### Load experimental data from bgee ####
  # listBgeeSpecies(ordering = 1)
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
         # new species with Bgee 15_0
         Turkey = { SPECIE_LAT <- "Meleagris_gallopavo"
         DATASET = "mgallopavo_gene_ensembl"}, 
         Goat = {SPECIE_LAT <- "Capra_hircus"
         DATASET = "chircus_gene_ensembl"},
         Sheep = {SPECIE_LAT <- "Ovis_aries"
         DATASET = "oaries_gene_ensembl"}, 
         Monkey_CrabEatingMacaque = {SPECIE_LAT <- "Macaca_fascicularis"
         DATASET = "mfascicularis_gene_ensembl"}
  )
  if (is.na(RELEASE)) { # if specific release is desired
    bgee <- Bgee$new(species = SPECIE_LAT, dataType = "rna_seq", sendStats = F) # other available data domains "affymetrix", "est", "in_situ"
  } else {
    bgee <- Bgee$new(species = SPECIE_LAT, dataType = "rna_seq", sendStats = F, release = RELEASE) 
  }
  
  # Set naming of data bases; listBgeeRelease()[1,"release"]
  BegeeRelease <- bgee$release 
  DB_Bgee <- paste0(bgee[["species"]], bgee[["sqlite_extension"]])
  DB_PKsim <- paste0("GENEDB_", tolower(SPECIE))
  if (ADME_ONLY) { DB_PKsim <- paste0(DB_PKsim,"_ADME_ONLY") } 
  if (TPM_ONLY) { DB_PKsim <- paste0(DB_PKsim,"_TPM_ONLY") }
  DB_PKsim <- paste0(DB_PKsim,"_BgeeRelease_", BegeeRelease,"_",Sys.Date(),".expressionDB") 
  
  PATH2PKsim <- paste0(PATH,"/", SPECIE,"/")
  PATH_DB_Bgee <- bgee[["pathToData"]]
  
  dir.create(PATH2PKsim, showWarnings = F)
  
  #### Connect to bgee data base and do formatting ####
  db_bgee_conn <- 
    dbConnect(RSQLite::SQLite(), 
              paste0(PATH_DB_Bgee,"/", gsub(x = DB_Bgee, pattern = " ", replacement = "_") ), 
              synchronous = 'off', cache_size = -1000)
  
  DB_Tables <- dbListTables(conn = db_bgee_conn)
  
  db_biomart_conn <- 
    dbConnect(RSQLite::SQLite(), 
              "BioMarts/All_Species_BioMarts.DB", 
              synchronous = 'off', cache_size = -1000)
  
  
  if (rlang::is_empty(DB_Tables)) { 
    # If data is not yet locally stored execute getData() to initiate download, 
    # currently the latest version getData() switched to a forced load of whole db content into ram. 
    # This can cause to errors when the function GeneratePKsimDB.R is called the first time,
    # re-excecuting the funciton usually resolves the issue (then DB is already downloaded).
    if (SPECIE_LAT == "Homo_sapiens") {
      # Access expression data table, either in local storage or as remote lazy table (connection to database)
      # Human data is implemented to not include SRP012682 "Genotype-Tissue Expression (GTEx) Common Fund Project"
      # >> single data source with > 60 GB! code can currently not process this amount of data... :-(
      rna_seq_selected <- getData(bgee, experimentId = 
                                    c("GSE57344","GSE62098","GSE43520","GSE30352","GSE58387","GSE58608","GSE64283","GSE30611"))
    } else { rna_seq_selected <- getData(bgee) } 
    
    # The default Bgee databases contain quotes causing erroer during the database generation and are removed to free up disc space and for clear naming
    dbExecute(conn = db_bgee_conn, statement = c("UPDATE rna_seq SET \"Stage.name\" = REPLACE(\"Stage.name\",'\"' , '' )" )   )
    dbExecute(conn = db_bgee_conn, statement = c("UPDATE rna_seq SET \"Anatomical.entity.name\" = REPLACE(\"Anatomical.entity.name\", '\"' , '' )" )   )
    dbExecute(conn = db_bgee_conn, statement = c("UPDATE rna_seq SET \"Anatomical.entity.ID\" = REPLACE(\"Anatomical.entity.ID\", 'UBERON:' , '' )" )   )
    dbExecute(conn = db_bgee_conn, statement = c("UPDATE rna_seq SET \"Stage.ID\" = REPLACE(\"Stage.ID\", 'UBERON:' , '' )" )   )
    if ( parse_number(BegeeRelease) >= 14 ) { # strain column only include in bgee > 14
      dbExecute(conn = db_bgee_conn, statement = c("UPDATE rna_seq SET \"Strain\" = REPLACE(\"Strain\", '\"' , '' )" )   )
    }
    
  } 
  
  if (COMPUTE_IN_RAM) { 
    # as local table in RAM
    print("Fetch expression data from local BgeeDB.")
    rna_seq_selected <- getData(bgee) 
  } else { 
    # as remote table
    rna_seq_selected <- tbl(db_bgee_conn, "rna_seq") 
    if (SPECIE_LAT == "Homo_sapiens") {
      rna_seq_selected <- rna_seq_selected %>% dplyr::filter(Experiment.ID != "SRP012682") %>% dplyr::compute() 
    } 
  }
  
  # Subset only high quality and as present declared genes form (inconsistent over different Bgee releases)
  if ( parse_number(BegeeRelease) >= 15 ) {
    rna_seq_selected <- rna_seq_selected %>% # bgee$
      filter(Detection.flag == "present") %>%
      dplyr::select(!c("Detection.flag","State.in.Bgee"))
  } else {
    rna_seq_selected <- rna_seq_selected %>% # bgee$
      filter(Detection.flag == "present") %>%
      filter(Detection.quality == "high quality") %>%
      dplyr::select(!c("Detection.flag","Detection.quality","State.in.Bgee"))
  }
  if ( parse_number(BegeeRelease) < 14 ) {
    rna_seq_selected <- rna_seq_selected %>% 
      mutate(Sex = "UNSPECIFIED") %>% mutate(Strain = SPECIE) %>% mutate(FPKM = NA) %>% mutate(TPM = NA)
  }
  
  if (INCLUDE_RPKM) {
    # Add RPKMs to the data set; for this gene length is needed, extracted form biomart annotations
    # How to derive RPKM: https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
    # https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
    universe <- as.data.frame(readRDS(file = paste0(PATH,"/BioMarts/",SPECIE,"_Annotations.rds")))
    
    # Get sequence length for each gen
    AnnotationTable  <- dbReadTable(db_biomart_conn, paste0(SPECIE,"_Annotations")) %>% 
      mutate(length = end_position - start_position) %>% 
      dplyr::select(VARIANT_NAME, length) %>% rename(Gene.ID = VARIANT_NAME)
    
    # add the annotation data to the data table
    rna_seq_selected <- rna_seq_selected %>% 
      dplyr::left_join(., AnnotationTable %>% distinct(), by = "Gene.ID") 
    
    # Calculate total counts per experiment
    TcountExp <- rna_seq_selected %>% 
      group_by(Library.ID) %>%
      summarise(SumedReadCount = sum(Read.count)) %>%
      dplyr::mutate(SumedReadCount_per_Milion = SumedReadCount / 10^6) %>%
      ungroup()
    
    # Add total counts per experiment to data set
    rna_seq_selected <- rna_seq_selected %>% 
      dplyr::left_join(TcountExp, by = "Library.ID", copy = T) %>% dplyr::compute()
    
    # Calculate RPKM for each experiment
    rna_seq_selected  <- rna_seq_selected %>% 
      dplyr::mutate(RPKM =  (Read.count / SumedReadCount_per_Milion) / (length / 10^3))
    
    rna_seq_selected <- rna_seq_selected %>%
      dplyr::select(!length) %>%
      dplyr::select(!SumedReadCount_per_Milion) %>%
      dplyr::select(!SumedReadCount)
    print("Read count were translated to RPKM values and added to Bgee data.")
  }
  
  # Rename columns to match PK-Sim requiered namespace
  rna_seq_selected <- rna_seq_selected %>%
    dplyr::rename(DATA_BASE = Experiment.ID,
                  DATA_BASE_REC_ID = Library.ID,
                  VARIANT_NAME = Gene.ID,
                  TISSUE = Anatomical.entity.name,
                  GENDER = Sex,
                  STRAIN = Strain,
                  AGE = Stage.name) 
  
  # Further format, rename and prepare data ####
  rna_seq_selected <- rna_seq_selected %>%
    dplyr::mutate(GENDER = toupper(GENDER)) %>%
    dplyr::mutate(GENDER = if_else(GENDER == "NA", "UNSPECIFIED", GENDER)) %>%
    dplyr::mutate(STRAIN = if_else(STRAIN == "NA", "UNSPECIFIED", STRAIN)) %>%
    dplyr::mutate(TISSUE = if_else(TISSUE == "NA", "UNSPECIFIED", TISSUE)) %>%
    dplyr::mutate(TISSUE = toupper(TISSUE)) %>%
    dplyr::mutate(STATE = toupper("NORMAL")) %>%
    dplyr::mutate(HEALTH_STATE = toupper(paste(STRAIN, STATE, sep = " ")) ) %>%
    dplyr::select(!c("Anatomical.entity.ID","Library.type","Stage.ID"))  %>% dplyr::compute()
  
  # Add IDs for database structre
  rna_seq_selected_annotated <- dplyr::left_join(rna_seq_selected  %>% 
                                                   dplyr::select(VARIANT_NAME) %>% distinct() %>%
                                                   dplyr::mutate(VARIANT_ID = row_number() ), 
                                                 rna_seq_selected )  %>% dplyr::compute()
  rna_seq_selected_annotated <- dplyr::left_join(rna_seq_selected_annotated  %>% 
                                                   dplyr::select(DATA_BASE_REC_ID) %>% distinct() %>%
                                                   dplyr::mutate(DATA_SOURCE_ID = row_number() ), 
                                                 rna_seq_selected_annotated ) %>% dplyr::compute()
  rna_seq_selected_annotated <- dplyr::left_join(rna_seq_selected_annotated %>%
                                                   dplyr::select(AGE) %>% distinct() %>%
                                                   dplyr::mutate(AGE_ID = row_number() ),
                                                 rna_seq_selected_annotated )  %>%
    dplyr::mutate(AGE = toupper(AGE)) %>%
    dplyr::mutate(STRAIN = toupper(STRAIN)) %>% dplyr::compute()
  rm(rna_seq_selected)
  
  #### Merge different expression measures (later Units in PK-Sim DB) ####
  # Here cut offs could be considered
  print("Merge expression data from multiple units.")
  TPM_Table <- rna_seq_selected_annotated %>% 
    dplyr::select(!one_of(c("FPKM","Read.count","RPKM")) ) %>% distinct() %>%
    dplyr::mutate(UNIT = "TPM") %>%
    dplyr::rename(SAMPLE_COUNT = TPM) %>%
    dplyr::filter(!is.na(SAMPLE_COUNT)) %>%
    #filter(SAMPLE_COUNT > CUT_OFF[["TPM"]]) %>%
    dplyr::mutate(DATA_BASE_REC_ID = paste0(DATA_BASE_REC_ID, "_", UNIT)) %>%
    group_by(DATA_BASE_REC_ID ) %>%
    dplyr::mutate(TOTAL_COUNT = sum(SAMPLE_COUNT)) %>%
    ungroup() %>% dplyr::compute()
  
  if (!TPM_ONLY) {
    FPKM_Table <- rna_seq_selected_annotated %>% 
      dplyr::select(!one_of(c("TPM","Read.count","RPKM")) ) %>% distinct() %>%
      dplyr::mutate(UNIT = "FPKM") %>%
      dplyr::rename(SAMPLE_COUNT = FPKM) %>%
      dplyr::filter(!is.na(SAMPLE_COUNT)) %>%
      #filter(SAMPLE_COUNT > CUT_OFF[["FPKM"]]) %>%
      dplyr::mutate(DATA_BASE_REC_ID = paste0(DATA_BASE_REC_ID, "_", UNIT)) %>%
      group_by(DATA_BASE_REC_ID ) %>%
      dplyr::mutate(TOTAL_COUNT = sum(SAMPLE_COUNT)) %>%
      ungroup()
    
    Read.count_Table <- rna_seq_selected_annotated %>% 
      dplyr::select(!one_of(c("FPKM","TPM","RPKM")) ) %>% distinct() %>% 
      dplyr::mutate(UNIT = "Read.count") %>%
      dplyr::rename(SAMPLE_COUNT = Read.count) %>%
      dplyr::filter(!is.na(SAMPLE_COUNT)) %>%
      # filter(SAMPLE_COUNT > CUT_OFF[["Read.count"]]) %>%
      dplyr::mutate(DATA_BASE_REC_ID = paste0(DATA_BASE_REC_ID, "_", UNIT)) %>%
      group_by(DATA_BASE_REC_ID ) %>%
      dplyr::mutate(TOTAL_COUNT = sum(SAMPLE_COUNT)) %>%
      ungroup()
    
    if (INCLUDE_RPKM) {
      RPKM_Table <- rna_seq_selected_annotated %>% 
        dplyr::select(!one_of(c("FPKM","TPM","Read.count")) ) %>% distinct() %>% 
        dplyr::mutate(UNIT = "RPKM") %>%
        dplyr::rename(SAMPLE_COUNT = RPKM) %>%
        dplyr::filter(!is.na(SAMPLE_COUNT)) %>%
        #filter(SAMPLE_COUNT > CUT_OFF[["RPKM"]]) %>%
        dplyr::mutate(DATA_BASE_REC_ID = paste0(DATA_BASE_REC_ID, "_", UNIT)) %>%
        group_by(DATA_BASE_REC_ID ) %>%
        dplyr::mutate(TOTAL_COUNT = sum(SAMPLE_COUNT)) %>%
        ungroup()
    }
    
    DataTableCombined <- dplyr::full_join(TPM_Table, FPKM_Table)
    
    DataTableCombined <- dplyr::full_join(DataTableCombined, Read.count_Table)
    
    if (INCLUDE_RPKM) {DataTableCombined <- dplyr::full_join(DataTableCombined, RPKM_Table) }
    if (!COMPUTE_IN_RAM) {DataTableCombined <- DataTableCombined %>% dplyr::compute() } 
  } else { 
    DataTableCombined <- TPM_Table
    if (!COMPUTE_IN_RAM) {DataTableCombined <- DataTableCombined %>% dplyr::compute() } 
  }
  
  ##### Add Gene IDs / Annotations needed for allow broad reach ####
  print("Annotete expression data with gene & protein identifiers.")
  if (ADME_ONLY) { 
    AnnotationTable  <- dbReadTable(db_biomart_conn, paste0(SPECIE,"_ADME")) %>% dplyr::select(-c("start_position","end_position")) 
  } else { 
    AnnotationTable  <- dbReadTable(db_biomart_conn, paste0(SPECIE,"_Annotations")) %>% dplyr::select(-c("start_position","end_position")) 
  }
  
  # For only homologous genes, add "HOMOLOG" prefix in Gene Symbol; This is necessary since PK-SIM gene names are based on Symbols
  if (COMPUTE_IN_RAM == F) {
    dbWriteTable(conn = db_bgee_conn, name = "AnnotationTable", value = AnnotationTable, overwrite = T)
    AnnotationTable <- tbl(db_bgee_conn, "AnnotationTable")
  } else {
    AnnotationTable <- AnnotationTable %>% dplyr::collect()
  }
  # Combine expression data and annotation information
  rm(TPM_Table)
  rm(rna_seq_selected_annotated)
  
  if (ADME_ONLY) {
    DataTableCombined <- dplyr::left_join(AnnotationTable %>% dplyr::select(VARIANT_NAME) %>% distinct(), DataTableCombined) %>% dplyr::filter(!is.na(DATA_SOURCE_ID)) %>% dplyr::arrange(DATA_SOURCE_ID,VARIANT_ID) %>% dplyr::compute()
  } 
  DataTableCombined <- dplyr::left_join(DataTableCombined, AnnotationTable, by = "VARIANT_NAME") %>% dplyr::arrange(DATA_SOURCE_ID,VARIANT_ID) %>% dplyr::compute()
  
  #### generate TAB tables from data base data frame ####
  print("Write data to TAB files")
  
  ### TAB_EXPRESSION_DATA_PROPERTIES #
  KEYS <- c("DATA_SOURCE_ID","DATA_BASE","TISSUE","HEALTH_STATE","GENDER","AGE") 
  TAB_EXPRESSION_DATA_PROPERTIES <- 
    DataTableCombined %>% dplyr::select(all_of(KEYS)) %>% distinct() 
  
  TAB_EXPRESSION_DATA_PROPERTIES <- TAB_EXPRESSION_DATA_PROPERTIES %>% 
    pivot_longer(!DATA_SOURCE_ID, names_to = "PROPERTY", values_to = "PROPERTY_VALUE") %>% distinct() 
  
  ###TAB_EXPRESSION_DATA_VALUES #
  KEYS <- c("VARIANT_ID","DATA_SOURCE_ID","SAMPLE_COUNT","TOTAL_COUNT","UNIT") 
  TAB_EXPRESSION_DATA_VALUES <- DataTableCombined %>% dplyr::select(KEYS) %>% distinct() %>% dplyr::arrange(VARIANT_ID, DATA_SOURCE_ID) 
  
  ### TAB_EXPRESSION_DATA_AGE_PROPERTIES #
  # optional age is extracted form the data sets, in bgee strings like "2-3 weeks old rat" are given
  TAB_EXPRESSION_DATA_AGE_PROPERTIES <- DataTableCombined %>% dplyr::select(AGE_ID, AGE) %>%
    distinct() %>% pivot_longer(!AGE_ID , names_to = "PROPERTY", values_to = "PROPERTY_VALUE")
  
  ### TAB_EXPRESSION_DATA_AGES #
  # for mice also Theiler Stage could be included: https://www.emouseatlas.org/emap/ema/theiler_stages/StageDefinition/stagedefinition.html#dpc
#  Numextract <- function(string){
#    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
#  }
  TAB_EXPRESSION_DATA_AGES <- DataTableCombined %>% dplyr::select(AGE_ID, AGE) %>% distinct() %>% collect() %>%
    dplyr::mutate(AGE_MIN = 0) %>%
    dplyr::mutate(AGE_MAX = 0) %>%
    # dplyr::mutate(AGE_MIN = if_else(stringr::str_extract(string = AGE, pattern = "^[\\d].*MONTH"), AGE_MIN / 12,
    #                                 if_else(stringr::str_extract(string = AGE_MIN, pattern = "^[\\d].*WEEK"), AGE_MAX / 52 * -1, AGE_MIN )
    # ) ) %>%
    dplyr::mutate(AGE_MEDIAN = 0) %>%
    dplyr::mutate(AGE_MEAN = 0) %>%
    dplyr::mutate(AGE_MEAN = rowSums(tibble(a = matrix(readr::parse_number(stringr::str_extract(string = AGE, pattern = "^[\\d].*DAY")) / 365, ncol = 1),
                                            b = matrix(readr::parse_number(stringr::str_extract(string = AGE, pattern = "^[\\d].*WEEK")) / 52, ncol = 1),
                                            c = matrix(readr::parse_number(stringr::str_extract(string = AGE, pattern = "^[\\d].*MONTH")) / 12, ncol = 1),
                                            d = matrix(readr::parse_number(stringr::str_extract(string = AGE, pattern = "^[\\d].*YEAR")) , ncol = 1) ), na.rm = T)
    ) %>%  
    dplyr::mutate(AGE_MIN = rowSums(tibble(a = matrix(stringr::str_detect(string = AGE, pattern = "FIRST DECADE") * 0, ncol = 1),
                                           b = matrix(stringr::str_detect(string = AGE, pattern = "SECOND DECADE") * 11, ncol = 1),
                                           c = matrix(stringr::str_detect(string = AGE, pattern = "THIRD DECADE") * 20 , ncol = 1),
                                           d = matrix(stringr::str_detect(string = AGE, pattern = "FOURTH DECADE") * 30 , ncol = 1),
                                           e = matrix(stringr::str_detect(string = AGE, pattern = "FIFTH DECADE") * 40 , ncol = 1),
                                           f = matrix(stringr::str_detect(string = AGE, pattern = "SIXTH DECADE") * 50 , ncol = 1),
                                           g = matrix(stringr::str_detect(string = AGE, pattern = "SEVENTH DECADE") * 60 , ncol = 1),
                                           h = matrix(stringr::str_detect(string = AGE, pattern = "EIGHTH DECADE") * 70 , ncol = 1),
                                           i = matrix(stringr::str_detect(string = AGE, pattern = "NINTH DECADE") * 80 , ncol = 1),
                                           j = matrix(stringr::str_detect(string = AGE, pattern = "TENTH DECADE") * 90 , ncol = 1)), na.rm = T)
    ) %>%
    dplyr::mutate(AGE_MAX = rowSums(tibble(a = matrix(stringr::str_detect(string = AGE, pattern = "FIRST DECADE") * 10, ncol = 1),
                                           b = matrix(stringr::str_detect(string = AGE, pattern = "SECOND DECADE") * 20, ncol = 1),
                                           c = matrix(stringr::str_detect(string = AGE, pattern = "THIRD DECADE") * 30 , ncol = 1),
                                           d = matrix(stringr::str_detect(string = AGE, pattern = "FOURTH DECADE") * 40 , ncol = 1),
                                           e = matrix(stringr::str_detect(string = AGE, pattern = "FIFTH DECADE") * 50 , ncol = 1),
                                           f = matrix(stringr::str_detect(string = AGE, pattern = "SIXTH DECADE") * 60 , ncol = 1),
                                           g = matrix(stringr::str_detect(string = AGE, pattern = "SEVENTH DECADE") * 70 , ncol = 1),
                                           h = matrix(stringr::str_detect(string = AGE, pattern = "EIGHTH DECADE") * 80 , ncol = 1),
                                           i = matrix(stringr::str_detect(string = AGE, pattern = "NINTH DECADE") * 90 , ncol = 1),
                                           j = matrix(stringr::str_detect(string = AGE, pattern = "TENTH DECADE") * 100 , ncol = 1)), na.rm = T)
    ) %>%
    dplyr::mutate(AGE_MIN = if_else(stringr::str_detect(string = AGE, pattern = "FERTILIZATION"), AGE_MIN*-1, AGE_MIN) ) %>% 
    dplyr::mutate(AGE_MAX = if_else(stringr::str_detect(string = AGE, pattern = "FERTILIZATION"), AGE_MAX*-1, AGE_MAX) ) %>%
    dplyr::mutate(AGE_MEDIAN = if_else(stringr::str_detect(string = AGE, pattern = "FERTILIZATION"), AGE_MEDIAN*-1, AGE_MEDIAN) ) %>%
    dplyr::mutate(AGE_MEAN = if_else(stringr::str_detect(string = AGE, pattern = "FERTILIZATION"), AGE_MEAN*-1, AGE_MEAN) ) 
  # dplyr::mutate(across(where(is.numeric), ~dplyr::na_if(., Inf)), dplyr::across(where(is.numeric), ~dplyr::na_if(., -Inf))) %>% replace(is.na(.), 0)
  
  ## Alternative implementation that ignores given age numbers
  # TAB_EXPRESSION_DATA_AGES <- DataTableCombined %>% dplyr::select(AGE_ID, AGE) %>% distinct() %>% collect() %>%
  #   dplyr::mutate(AGE_MIN = 0) %>% dplyr::mutate(AGE_MAX = 0) %>% dplyr::mutate(AGE_MEAN = 0) %>%  dplyr::mutate(AGE_MEDIAN = 0)
  
  ### TAB_EXPRESSION_DATA_BASES #
  annotation_bgee <- getAnnotation(bgee)
  DataAnnotation <- as_tibble(unique(annotation_bgee[[2]][, c("Experiment.ID","Data.source.URL")]))
  DataAnnotation <- DataAnnotation %>% 
    dplyr::rename(DATA_BASE = Experiment.ID) %>% 
    dplyr::rename(URL = Data.source.URL) #%>% rename(NAME = Experiment.name) 
  AllDB <- DataTableCombined %>% dplyr::select(DATA_BASE) %>% distinct()
  TAB_EXPRESSION_DATA_BASES <- DataTableCombined %>% dplyr::select(DATA_BASE) %>% distinct() 
  
  if (!COMPUTE_IN_RAM) {
    dbWriteTable(conn = db_bgee_conn, name = "DataAnnotation", value = DataAnnotation, overwrite = T)
    DataAnnotation <- tbl(db_bgee_conn, "DataAnnotation") 
  } 
  TAB_EXPRESSION_DATA_BASES <- dplyr::left_join(TAB_EXPRESSION_DATA_BASES, DataAnnotation)
  
  ### TAB_EXPRESSION_DATA_GENDER_PROPERTIES #
  TAB_EXPRESSION_DATA_GENDER_PROPERTIES <- DataTableCombined %>% dplyr::select(GENDER,AGE,TISSUE) %>% distinct() %>% 
    dplyr::mutate(DEVELOPMENTAL_STAGE_SOURCE = "-") %>% 
    pivot_longer(!GENDER , names_to = "PROPERTY", values_to = "PROPERTY_VALUE") %>% distinct()  
  
  ### TAB_EXPRESSION_DATA_GENDERS #
  TAB_EXPRESSION_DATA_GENDERS <- DataTableCombined %>% dplyr::select(GENDER) %>% distinct() %>%
    dplyr::mutate(INFORMATION = paste0("A ",GENDER, " individual or population")) 
  
  ### TAB_EXPRESSION_DATA_HEALTH_STATE #
  TAB_EXPRESSION_DATA_HEALTH_STATE <- 
    DataTableCombined %>% dplyr::select(HEALTH_STATE) %>% distinct() %>%
    dplyr::mutate(INFORMATION = paste("This refers to from an ", HEALTH_STATE," individual.", sep = "")) 
  
  ### TAB_EXPRESSION_DATA_HEALTH_STATE_PROPERTIES #
  TAB_EXPRESSION_DATA_HEALTH_STATE_PROPERTIES <- 
    DataTableCombined %>% dplyr::select(HEALTH_STATE) %>% distinct() %>% 
    dplyr::mutate(PROPERTY = "HEALTH_STATE") %>% mutate(PROPERTY_VALUE = HEALTH_STATE) 
  
  ### TAB_EXPRESSION_DATA_RECORDS #
  KEYS <-  c("DATA_SOURCE_ID","DATA_BASE_REC_ID","TISSUE","HEALTH_STATE","GENDER","AGE_ID") 
  assign("TODAY", as.character(Sys.Date()) )
  TAB_EXPRESSION_DATA_RECORDS <- DataTableCombined %>% dplyr::select(KEYS) %>% distinct() %>% 
    dplyr::mutate(SAMPLE_SOURCE = "TISSUE") %>%
    dplyr::mutate(DATA_BASE = "RNAseq") %>%  
    dplyr::mutate(LAST_REFRESH_DATE =  TODAY) %>% 
    #mutate(LAST_REFRESH_DATE = as.character(as.POSIXlt(strptime(Sys.Date(), format = "%Y-%m-%d")))) %>% 
    #collect() %>% 
    dplyr::relocate(AGE_ID , .after = last_col())
  
  ### TAB_EXPRESSION_DATA_SAMPLE_SOURCE_PROPERTIES #
  TAB_EXPRESSION_DATA_SAMPLE_SOURCE_PROPERTIES <- 
    DataTableCombined %>% dplyr::select(TISSUE) %>% distinct() %>% 
    dplyr::rename(PROPERTY_VALUE = TISSUE) %>% dplyr::mutate(PROPERTY = "TISSUE_SOURCE") %>% 
    dplyr::mutate(SAMPLE_SOURCE = "TISSUE") %>% dplyr::arrange(PROPERTY_VALUE)
  
  ### TAB_EXPRESSION_DATA_SAMPLE_SOURCES #
  TAB_EXPRESSION_DATA_SAMPLE_SOURCES <- 
    tibble(SAMPLE_SOURCE = c("CELL LINE","PRIMARY CULTURE","TISSUE","UNSPECIFIED"),
           INFORMATION = c("The sample source is a cell line.",
                           "The sample source is a tissue culture started from cells, tissues, or organs taken directly from the organism.",
                           "The sample source is tissue.",
                           "The sample source is unknown."))
  
  ### TAB_EXPRESSION_DATA_TISSUE_PROPERTIES #
  TAB_EXPRESSION_DATA_TISSUE_PROPERTIES <- 
    DataTableCombined %>% dplyr::select(TISSUE) %>% distinct() %>% 
    dplyr::mutate(PROPERTY_VALUE = TISSUE) %>% 
    dplyr::mutate(PROPERTY = "TISSUE") %>% dplyr::arrange(PROPERTY_VALUE)
  
  ### TAB_EXPRESSION_DATA_TISSUES #
  TAB_EXPRESSION_DATA_TISSUES <-  
    DataTableCombined %>% dplyr::select(TISSUE) %>% distinct() %>% 
    dplyr::mutate(INFORMATION = "Organ") %>% dplyr::arrange(TISSUE)
  
  ### TAB_EXPRESSION_DATA_UNITS #
  TAB_EXPRESSION_DATA_UNITS <- 
    DataTableCombined %>% dplyr::select(UNIT) %>% distinct() %>% 
    dplyr::mutate(INFORMATION = "NGS data")
  
  ### TAB_GENE_NAMES #
  # linking the unique variant ID to the different identifier
  KEYS <- c("VARIANT_ID","VARIANT_NAME","SYMBOL","OFFICIAL_FULL_NAME","ENTREZID","SYNONYM","PROTEIN_ID","PREFFERED_NAME","OTHER_NAME")
  if (SPECIE != "Human" & SPECIE != "Sheep") {
    KEYS <- c(KEYS, "HOMOLOG", "HOMOLOG_SYMBOL")
  }
  TAB_GENE_NAMES <- DataTableCombined %>% dplyr::select(any_of(KEYS)) %>% distinct() %>%
    dplyr::rename(GENE_ID = ENTREZID) %>% mutate(GENE_ID = as.character(GENE_ID)) %>% distinct() %>% arrange(VARIANT_ID )
  
  TAB_GENE_NAMES <- TAB_GENE_NAMES %>% 
    pivot_longer(!VARIANT_ID, names_to = "NAME_TYPE", values_to = "GENE_NAME") %>% 
    dplyr::filter(GENE_NAME != "") %>% dplyr::rename(GENE_ID = VARIANT_ID) %>% distinct()
  
  ### TAB_GENE_NAME_TYPES #
  TAB_GENE_NAME_TYPES <- TAB_GENE_NAMES %>% dplyr::select(NAME_TYPE) %>% distinct() %>% 
    dplyr::mutate(INFORMATION = paste0("The identifier '", NAME_TYPE,"' is based on the Cran R biomaRt package.")) %>% compute()
  
  ### TAB_GENE_VARIANTS #
  TAB_GENE_VARIANTS <- DataTableCombined %>% dplyr::select(VARIANT_ID,VARIANT_NAME) %>% 
    distinct() %>% dplyr::mutate(GENE_ID = VARIANT_ID) 
  
  ### TAB_GENES #
  TAB_GENES <-  DataTableCombined %>% dplyr::select(VARIANT_ID) %>% distinct() %>% 
    dplyr::rename(GENE_ID = VARIANT_ID) 
  
  ### TAB_GLOBAL_STATISTICS #
  # function in original Access DB == AVG: 10^mean(Logarithmus([SAMPLE_COUNT]/[TOTAL_COUNT])/Logarithmus(10))
  TAB_GLOBAL_STATISTICS <- DataTableCombined %>% dplyr::select(UNIT, SAMPLE_COUNT,TOTAL_COUNT) %>% 
    group_by(UNIT) %>% dplyr::mutate(AVG = if_else(condition = UNIT %in% c("RPKM","FPKM"), true = 10^mean(log(SAMPLE_COUNT/TOTAL_COUNT)/log(10),na.rm = T), 
                                                   false = 10^mean(log(SAMPLE_COUNT)/log(10)))) %>% dplyr::select(UNIT, AVG) %>% distinct() 
  
  ##### Write data into  PK-Sim expression database #####
  assign(x = "TAB_CONTAINER_TISSUE", value = tibble(read.table("Code/TAB_CONTAINER_TISSUE.txt", header = 1, sep = "\t")))
  assign(x = "TAB_DTS_PROPERTIES", value = tibble(read.table("Code/TAB_DTS_PROPERTIES.txt", header = 1, sep = "\t")))
  Files <- ls(pattern = "^TAB_.{4}")
  
  #Connect to SQLite db
  db_PKsim_conn <- dbConnect(RSQLite::SQLite(), paste0(PATH2PKsim, DB_PKsim), synchronous = NULL)
  DB_Tables <- dbListTables(db_PKsim_conn)
  
  print("Write data to SQLite data base.")
  for (TAB in DB_Tables) {
    ifelse(isEmpty(grep("QRY", TAB)), dbRemoveTable(db_PKsim_conn, TAB), dbExecute(conn = db_PKsim_conn, paste0("DROP VIEW ", TAB,";")) )
  }
  
  # Load modified tab files and write to db
  for (i in 1:length(Files)) {
    TAB <- Files[i]
    # info output to screen
    print(paste0("Write \'", tolower(TAB), "\' to \'", SPECIE, "\' PK-Sim expression database"))
    
    tmp_table_names <- TYPE[[which(TAB_s %in% tolower(TAB))]]
    names(tmp_table_names) <- colnames(get(TAB))
    
    dbWriteTable(conn = db_PKsim_conn, name = tolower(TAB), value = get(TAB) %>% collect() , row.names = FALSE, 
                 overwrite = TRUE, field.types = tmp_table_names)
  }
  
  for (i in 1:length(VIEW)) { # if error occurs generate DB and show warning
    dbExecute(db_PKsim_conn, VIEW[i])
  } 
  
  for (i in 1:length(INDIZES)) { # if error occurs generate DB and show warning
    dbExecute(db_PKsim_conn, INDIZES[i])
  }  
  
  dbExecute(conn = db_PKsim_conn, statement = "VACUUM") # clears pre-allocated disc space for db
  dbDisconnect(db_PKsim_conn)
  dbDisconnect(db_bgee_conn)
  gc() # clears used memory
}
