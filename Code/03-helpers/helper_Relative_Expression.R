# Utility: extract expression values from a PK-Sim expression DB by Ensembl ID
# Relates to: https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/issues/6
#
# Usage (source this file after opening a DB connection):
#   conn <- DBI::dbConnect(RSQLite::SQLite(), "PK-Sim DBs/Human/GENEDB_human_ADME_ONLY_BgeeRelease_15_2.expressionDB")
#   result <- get_expression_by_ensembl_id("ENSG00000139618", conn)
#   DBI::dbDisconnect(conn)

#' Extract relative expression values for a gene using Ensembl ID
#' @param ensembl_id Character. The Ensembl gene ID to query (e.g., "ENSG00000139618")
#' @param db_conn SQLite connection to PKSim expression DB
#' @return Data frame with expression values for the gene
get_expression_by_ensembl_id <- function(ensembl_id, db_conn) {
  # Query tab_gene_names to get variant_id for the Ensembl ID
  # name_type 'ENSEMBL' stores Ensembl gene identifiers
  query_variant <- "SELECT gv.variant_id FROM tab_gene_names gn
     JOIN tab_gene_variants gv ON gn.gene_id = gv.gene_id
     WHERE gn.gene_name = ? AND gn.name_type = 'ENSEMBL'"
  variant_ids <- DBI::dbGetQuery(db_conn, query_variant, params = list(ensembl_id))
  if (nrow(variant_ids) == 0) {
    stop(sprintf("No expression data found for Ensembl ID: %s", ensembl_id))
  }
  # Query tab_expression_data_values for expression values using parameterized IN clause
  placeholders <- paste(rep("?", nrow(variant_ids)), collapse = ",")
  query_expr <- paste0(
    "SELECT * FROM tab_expression_data_values WHERE variant_id IN (", placeholders, ")"
  )
  expr_values <- DBI::dbGetQuery(db_conn, query_expr, params = as.list(variant_ids$variant_id))
  return(expr_values)
}

