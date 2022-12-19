PrepareBioMarts <- function(SPECIE = "Rat") {
  ### Function is designed to download species gene annotation and store the information (as RDS & DB files)
  ### The information is needed to:
  ###             1. Estimate RPKM from other RNAseq data (gene length needed; https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/)
  ###             2. Annotate a gene to enable search over multiple identifiers 
  ###                (e.g. NCBI gene id: 1543, ensembl id: ENSG00000140465, symbol: CYP1A1, etc.)
  ###             3. Map ortholog genes (homologs between species)
  ###
  ### Data storage ideally in an data base to enable lazy tibble operations
  ### List of human to ortholog animal ensembl id mapping (downloaded from: https://www.ensembl.org/biomart/martview)
  ### Length of a gene is estimated form start and stop region on chromosome
  ### A subset of genes is prepared including only AMDE relevant genes (reduce data load)
  ###
  ### ADME gene subset is defined via gene symbols compiled from different sources:
  ###     1. infos from: https://www.physiome.org/Course/2015winBioen498e/wk4/KellDobson08.pdf
  ###     2. pharmaadme.org
  ###     3. http://www.membranetransport.org/
  ###     4. http://www.tcdb.org/
  ##
  #
  ##
  ### List of available species:
  ## ensembl <- useMart("ensembl")
  ## ListofEnsemblSpecies <- listDatasets(ensembl)
  #
  # Notably Sheep homolog genes are not available from oaries_gene_ensembl
  # Alternative homology mappings via: OMA # https://bioconductor.org/packages/release/bioc/html/OmaDB.html
  
  #### Dependencies ####
  require("biomaRt")
  require("dplyr")
  #require("data.table")
  require("RSQLite")
  require("stringr")
  
  #### Preparations ####
  # adding the various identifier and information is extremely memory intensive and might only work with a 64-bit R version! 
  ALL_SPECIE <- c("Human","Monkey","Minipig","Dog","Mouse","Rat","Rabbit",
                  "Zebrafish","Cattle","Horse","Cat","GuineaPig","Chicken",
                  "Goat","Sheep","Turkey","Monkey_CrabEatingMacaque")
  if (!SPECIE %in% ALL_SPECIE) {
    stop(paste0("Given SPECIE: '", SPECIE, "' not in list. \n Must be one of: ", paste(ALL_SPECIE,collapse = " , ") ))
  }
  
  # Create annotation database, if exist only connection is made
  dir.create("./BioMarts/", showWarnings = F)
  Conn <- dbConnect(drv = RSQLite::SQLite(), paste0("./BioMarts/All_Species_BioMarts.DB"), synchronous = NULL)
  if (!dbExistsTable(conn = Conn, name = "Human_ADME") & SPECIE != "Human") {
    print(paste0("Human annotation do not exist in 'All_Species_BioMarts.DB' \n Human data is loaded first.") )
    PrepareMarts("Human")
  } 
  
  #### Load biomart data and store ####
  print(paste0("Loading ", SPECIE," gene inforamtion from biomart mirror."))
  
  switch(SPECIE,
         # Updates in Ensembl reference geneomes can cause  missmatches between organism Bgee Ensembl ID and reference IDs Ensembl IDs from biomaRt, fix is for not witchting to previous ensembl verions for specific species
         Human      = {ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")#, ensemblRedirec = FALSE)   
         }, Monkey  = {ensembl = useMart("ensembl", dataset = "mmulatta_gene_ensembl")   
         }, Minipig = {ensembl = useMart("ensembl", dataset = "sscrofa_gene_ensembl")    
         }, Dog     = {ensembl = useMart("ensembl", dataset = "clfamiliaris_gene_ensembl", host = "https://may2021.archive.ensembl.org") # naming error in previous biomart version can cause error --> fix is to use previous assambyl version
         }, Rat     = {ensembl = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
         }, Mouse   = {ensembl = useMart("ensembl", dataset = "mmusculus_gene_ensembl")  
         }, Rabbit  = {ensembl = useMart("ensembl", dataset = "ocuniculus_gene_ensembl") 
         }, Zebrafish  = {ensembl = useMart("ensembl", dataset = "drerio_gene_ensembl") 
         }, Cattle  = {ensembl = useMart("ensembl", dataset = "btaurus_gene_ensembl", host = "https://apr2022.archive.ensembl.org") 
         }, Horse  = {ensembl = useMart("ensembl", dataset = "ecaballus_gene_ensembl") 
         }, Cat  = {ensembl = useMart("ensembl", dataset = "fcatus_gene_ensembl") 
         }, GuineaPig  = {ensembl = useMart("ensembl", dataset = "cporcellus_gene_ensembl") 
         }, Chicken  = {ensembl = useMart("ensembl", dataset = "ggallus_gene_ensembl", host = "https://apr2022.archive.ensembl.org") # GRCg6a	GCA_000002315.5
         }, Turkey  = {ensembl = useMart("ensembl", dataset = "mgallopavo_gene_ensembl", host = "https://apr2022.archive.ensembl.org") 
         }, Goat  = {ensembl = useMart("ensembl", dataset = "chircus_gene_ensembl", host = "https://apr2022.archive.ensembl.org") 
         }, Sheep  = {ensembl = useMart("ensembl", dataset = "oaries_gene_ensembl", host = "https://apr2022.archive.ensembl.org") # oarambouillet_gene_ensembl has human homology mappings
         }, Monkey_CrabEatingMacaque  = {ensembl = useMart("ensembl", dataset = "mfascicularis_gene_ensembl") 
         }
  )
  
  #### Select annotation information and identifiers; mapping ortholog genes  ####
  print(paste0("Fetiching ", SPECIE," gene identifiers and annotations."))
  
  # Get SPECIE specific gene annotations
  SPECIE_ANNOTATION <- getBM(attributes = c('ensembl_gene_id',"description","entrezgene_id",
                                            "external_synonym","external_gene_name", "wikigene_name",
                                            "uniprot_gn_symbol"), mart = ensembl) # hgnc_symbol
  # Add SPECIE specific gene lengths and map human ortholog onto animal genes
  SPECIE_ANNOTATION <- left_join(SPECIE_ANNOTATION, getBM(attributes = c("ensembl_gene_id","start_position","end_position"), 
                                                          mart = ensembl), 
                                 by =  "ensembl_gene_id")
  SPECIE_ANNOTATION <- SPECIE_ANNOTATION %>% rename(., 
                                                    VARIANT_NAME = ensembl_gene_id, 
                                                    ENTREZID = entrezgene_id, 
                                                    PREFFERED_NAME = external_gene_name, 
                                                    SYMBOL = uniprot_gn_symbol, 
                                                    OFFICIAL_FULL_NAME = description, 
                                                    SYNONYM = external_synonym,
                                                    OTHER_NAME = wikigene_name)
  
  if (SPECIE != "Human" & SPECIE != "Sheep") { 
    print(paste0("Link ", SPECIE," gene identifiers to human orthologs based on ensembl IDs."))
    SPECIE_ANNOTATION <- left_join(SPECIE_ANNOTATION, 
                                   getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_associated_gene_name",
                                                        "hsapiens_homolog_ensembl_gene","hsapiens_homolog_perc_id_r1"),
                                         mart = ensembl), by =  c("VARIANT_NAME" = "ensembl_gene_id") )# "hsapiens_homolog_orthology_type"
    
    SPECIE_ANNOTATION <- SPECIE_ANNOTATION %>% 
      rename(., HOMOLOG_SYMBOL = hsapiens_homolog_associated_gene_name, 
             HOMOLOG = hsapiens_homolog_ensembl_gene, 
             HOMOLOGY_PERCENT = hsapiens_homolog_perc_id_r1)
  }
  # String cleaning 
  SPECIE_ANNOTATION <- tibble(SPECIE_ANNOTATION) %>% arrange(VARIANT_NAME)
  SPECIE_ANNOTATION <- SPECIE_ANNOTATION %>% 
    mutate(OFFICIAL_FULL_NAME = str_remove(string = OFFICIAL_FULL_NAME, pattern = "\\[Source:.*"))
  
  # Write annotation table to database
  dbWriteTable(conn = Conn, name = paste0(SPECIE,"_Annotations"), value = as_tibble(SPECIE_ANNOTATION) , overwrite = T)
  
  #### Only ADME genes ####
  # List of ADME gene abbreviations
  ADMEgene <- c("^ABC","^ADH","^AHR","^ALD","^ALDH","^ANXA","^AOX","^ARN","^ARNT","^ARS","^ATP",
                "^CACN","^CAT","^CBR","^CDA","^CES","^CFT","^CHS","^CHST","^CYB","^CYP","^DDO",
                "^DHR","^DHRS","^DPE","^DPY","^DPYD","^EPH","^EPHX","^FMO","^GPX","^GSR","^GSS",
                "^GST","^HAG","^HNF","^HNM","^HNMT","^HSD","^IAP","^KCN","^LOC","^MAT","^MET",
                "^MGS","^MGST","^MPO","^NAT","^NNM","^NNMT","^NOS","^NR1","^NUDT","^PDE","^PIAS","^PLG","^PNM",
                "^PNMT","^PON","^POR","^PPA","^PPAR","^PPARA","^RXR","^SER","^SLC","^SCN","^SOD",
                "^SUL","^SULF","^SULT","^SQST","^TAP","^TPM","^TPMT","^UGT","^UROC","^VKOR","^XDH","^XRC")
  
  # Select only ADME genes
  if (SPECIE != "Human" & SPECIE != "Sheep") {
    SPECIE_ANNOTATION <- SPECIE_ANNOTATION %>% filter(grepl(paste(ADMEgene, collapse = "|"), SYNONYM ) |
                                                        grepl(paste(ADMEgene, collapse = "|"), SYMBOL ) |  
                                                        grepl(paste(ADMEgene, collapse = "|"), HOMOLOG_SYMBOL ) |
                                                        grepl(paste(ADMEgene, collapse = "|"), PREFFERED_NAME ) |
                                                        grepl(paste(ADMEgene, collapse = "|"), OTHER_NAME ) 
                                                      
    ) %>% distinct() %>% arrange(VARIANT_NAME)
  } else {
    SPECIE_ANNOTATION <- SPECIE_ANNOTATION %>% filter(grepl(paste(ADMEgene, collapse = "|"), SYNONYM ) |
                                                        grepl(paste(ADMEgene, collapse = "|"), SYMBOL ) |  
                                                        grepl(paste(ADMEgene, collapse = "|"), PREFFERED_NAME ) |
                                                        grepl(paste(ADMEgene, collapse = "|"), OTHER_NAME ) 
    ) %>% distinct() %>% arrange(VARIANT_NAME)
  }
  # Save to database, compress database and close connection
  dbWriteTable(conn = Conn, name = paste0(SPECIE,"_ADME"), value = as_tibble(SPECIE_ANNOTATION), overwrite = T)
  dbExecute(conn = Conn, statement = "VACUUM") # reduce storage space
  dbDisconnect(Conn)
}
