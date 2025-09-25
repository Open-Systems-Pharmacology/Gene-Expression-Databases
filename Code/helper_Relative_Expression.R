# make a function that extract the realtive expression from PKSim expression DB given a gene name as input
# check issue: https://github.com/Open-Systems-Pharmacology/Gene-Expression-Databases/issues/6

# access PKsim DB
db_PKsim_conn <- DBI::dbConnect(RSQLite::SQLite(),
  paste0(PATH2PKsim, DB_PKsim),
  synchronous = "off",
  cache_size = -1000
)

get_expression_by_ensembl_id(ENSG00000139618, )
# build query to extract expression values for a gene using ensemble ID
#' Extract relative expression values for a gene using Ensembl ID
#' @param ensembl_id Character. The Ensembl gene ID to query (e.g., "ENSG00000139618")
#' @param db_conn SQLite connection to PKSim expression DB
#' @return Data frame with expression values for the gene
get_expression_by_ensembl_id <- function(ensembl_id, db_conn) {
  # Query tab_gene_names to get variant_id for the Ensembl ID
  query_variant <- sprintf(
    "SELECT gene_id FROM tab_gene_names WHERE gene_name = '%s' AND name_type = 'variant_name'",
    ensembl_id
  )
  variant_ids <- DBI::dbGetQuery(db_conn, query_variant)
  if (nrow(variant_ids) == 0) {
    stop(sprintf("No variant_id found for Ensembl ID: %s", ensembl_id))
  }
  # Query tab_expression_data_values for expression values
  query_expr <- sprintf(
    "SELECT * FROM tab_expression_data_values WHERE variant_id IN (%s)",
    paste(variant_ids$gene_id, collapse = ",")
  )
  expr_values <- DBI::dbGetQuery(db_conn, query_expr)
  return(expr_values)
}

# Example usage:
# result <- get_expression_by_ensembl_id("ENSG00000139618", db_PKsim_conn)
