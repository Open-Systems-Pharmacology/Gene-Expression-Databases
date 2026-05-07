# Helper: get proteins by name
get_proteins_by_name <- function(name, conn) {
  tab_gene_names <- RSQLite::dbReadTable(conn, "tab_gene_names")
  tab_gene_variants <- RSQLite::dbReadTable(conn, "tab_gene_variants")
  tab_expression_data_values <- RSQLite::dbReadTable(conn, "tab_expression_data_values")

  proteinsbyname <- tab_gene_names |>
    # dplyr::filter(grepl(name, gene_name)) |>
    dplyr::filter(stringr::str_detect(gene_name, paste(stringr::str_escape(name), collapse = "|"))) |>
    dplyr::group_by(gene_id) |>
    dplyr::summarise(
      gene_name = min(gene_name),
      name_type = min(name_type),
      symbol = dplyr::coalesce(
        dplyr::first(gene_name[name_type == "SYMBOL"], default = NA_character_),
        NA_character_
      ),
      gene_id_str = dplyr::coalesce(
        dplyr::first(gene_name[name_type == "GENE_ID"], default = NA_character_),
        NA_character_
      ),
      official_full_name = dplyr::coalesce(
        dplyr::first(gene_name[name_type == "OFFICIAL_FULL_NAME"], default = NA_character_),
        NA_character_
      )
    ) |>
    dplyr::ungroup()

  has_data_ids <- tab_gene_variants |>
    dplyr::filter(variant_id %in% tab_expression_data_values$variant_id) |>
    dplyr::pull(gene_id) |>
    unique()

  proteinsbyname |>
    dplyr::mutate(
      has_data = as.integer(gene_id %in% has_data_ids)
    )
}

# Helper: get expression data by gene id
get_expression_data_by_gene_id <- function(
    P_ID, conn,
    records_filter = NULL,
    ages_min = NULL,
    ages_max = NULL,
    unit_filter = NULL) {
  tab_genes <- dplyr::tbl(conn, "tab_genes")
  tab_gene_variants <- dplyr::tbl(conn, "tab_gene_variants")
  tab_expression_data_records <- dplyr::tbl(conn, "tab_expression_data_records")
  tab_expression_data_ages <- dplyr::tbl(conn, "tab_expression_data_ages")
  tab_expression_data_values <- dplyr::tbl(conn, "tab_expression_data_values")
  tab_global_statistics <- dplyr::tbl(conn, "tab_global_statistics")

  # Apply filters if provided
  if (!is.null(records_filter)) {
    tab_expression_data_records <- tab_expression_data_records |>
      # dplyr::filter(!stringr::str_detect(data_base_rec_id, stringr::fixed(records_filter)))
      # dplyr::filter(!stringr::str_detect(data_base_rec_id, records_filter))
      dplyr::filter(!stringr::str_detect(data_base_rec_id, stringr::fixed(records_filter)))
  }
  if (!is.null(ages_min) & !is.null(ages_max)) {
    tab_expression_data_ages <- tab_expression_data_ages |>
      dplyr::filter(age_min == 0 | age_min >= ages_min,
                    age_max == 0 | age_max <= ages_max)
  }
  if (!is.null(unit_filter)) {
    tab_expression_data_values <- tab_expression_data_values |> dplyr::filter(unit == unit_filter)
  }

  # Query 1: with age info
  query1 <- tab_genes |>
    dplyr::filter(gene_id %in% !!P_ID) |>
    dplyr::inner_join(tab_gene_variants, by = "gene_id") |>
    dplyr::inner_join(tab_expression_data_values, by = "variant_id") |>
    dplyr::inner_join(tab_expression_data_records, by = "data_source_id") |>
    dplyr::inner_join(tab_expression_data_ages, by = "age_id") |>
    dplyr::inner_join(tab_global_statistics, by = "unit") |>
    dplyr::filter(!is.na(tissue)) |>
    dplyr::mutate(
      ratio = sample_count / total_count,
      norm_value = (sample_count / total_count) / avg
    ) |>
    dplyr::select(
      gene_id, variant_name, data_base, data_base_rec_id, gender, tissue, health_state, sample_source,
      age_min, age_max, sample_count, total_count, ratio, norm_value, unit
    )

  # Query 2: without age info
  query2 <- tab_genes |>
    dplyr::filter(gene_id %in% !!P_ID) |>
    dplyr::inner_join(tab_gene_variants, by = "gene_id") |>
    dplyr::inner_join(tab_expression_data_values, by = "variant_id") |>
    dplyr::inner_join(tab_expression_data_records |> dplyr::filter(is.na(age_id)), by = "data_source_id") |>
    dplyr::inner_join(tab_global_statistics, by = "unit") |>
    dplyr::filter(!is.na(tissue)) |>
    dplyr::mutate(
      age_min = NA_real_,
      age_max = NA_real_,
      ratio = sample_count / total_count,
      norm_value = (sample_count / total_count) / avg
    ) |>
    dplyr::select(
      gene_id, variant_name, data_base, data_base_rec_id, gender, tissue, health_state, sample_source,
      age_min, age_max, sample_count, total_count, ratio, norm_value, unit
    )

  dplyr::bind_rows(
    query1 |> dplyr::collect(),
    query2 |> dplyr::collect()
  )
}

# Helper: get container tissue mapping
get_container_tissue_mapping <- function(conn) {
  tab_container_tissue <- RSQLite::dbReadTable(conn, "tab_container_tissue")
  tab_container_tissue |> dplyr::select(container, tissue)
}

# Helper: get hint information
get_hint <- function(conn, table_name, column, value) {
  table <- RSQLite::dbReadTable(conn, table_name)
  info <- table |> dplyr::filter(.data[[column]] == value)
  if (nrow(info) == 0) {
    return("")
  }
  info$INFORMATION[1]
}
# Example usage:
# gender_hint <- get_hint(tab_gender_hint, "gender", "male")
# tissue_hint <- get_hint(tab_tissue_hint, "tissue", "liver")
# health_state_hint <- get_hint(tab_health_state_hint, "health_state", "healthy")
# sample_source_hint <- get_hint(tab_sample_source_hint, "sample_source", "blood")
# unit_hint <- get_hint(tab_unit_hint, "unit", "TPM")
# name_type_hint <- get_hint(tab_name_type_hint, "name_type", "symbol")

# Helper: get database record properties
get_database_rec_properties <- function(conn, database, rec_id) {
  tab_database_rec_properties <- RSQLite::dbReadTable(conn, "tab_database_rec_properties")
  tab_database_rec_properties |>
    dplyr::filter(data_base == database, data_base_rec_id == rec_id) |>
    dplyr::mutate(property_string = paste(PROPERTY, PROPERTY_VALUE, sep = ": ")) |>
    dplyr::pull(property_string)
}

# Helper: get database record infos
get_database_rec_infos <- function(conn, database, rec_id) {
  tab_database_rec_info <- RSQLite::dbReadTable(conn, "tab_database_rec_info")
  tab_database_rec_info |>
    dplyr::filter(data_base == database, data_base_rec_id == rec_id) |>
    dplyr::select(-data_base, -data_base_rec_id) |>
    tidyr::pivot_longer(everything(), names_to = "column", values_to = "value") |>
    dplyr::mutate(info_string = paste(column, value, sep = ": ")) |>
    dplyr::pull(info_string)
}

# Helper: validate query columns
validate_query <- function(conn, table_name, columns) {
  table <- RSQLite::dbReadTable(conn, table_name)
  missing <- setdiff(columns, colnames(table))
  if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse = ", ")))
}

# Helper: clear cache (if using environments or lists for caching)
clear_cache <- function(cache_env) {
  rm(list = ls(envir = cache_env), envir = cache_env)
}
