# Comparison between old and new OSP-Expression database for humans #####
# Goal: Compare gene expression between new and old expression DB
# Task: Download explicit study and extract same explicit study from new Db

# free workspace
rm(list = ls())
setwd(here::here())

source("./Code/helper_SQL_Queries.R")
# 2.0 Load Old QSP Expression DB ####
# DB connection
genedb_human_old <- RSQLite::dbConnect(
  drv = RSQLite::SQLite(),
  "GENEDB_human.expressionDb"
)

enzymes <- xlsx::read.xlsx(
  file = "Code/Qualification/Proteins4Validation.xlsx",
  sheetIndex = 1,
  header = FALSE
) |>
  dplyr::rename("SYMBOL" = "X1")
enzymes <- enzymes |>
  dplyr::mutate(SYMBOL = trimws(SYMBOL, "both")) |>
  unlist() |>
  unname()

x_old <- get_proteins_by_name(name = enzymes, conn = genedb_human_old)

expression_values_old <- get_expression_data_by_gene_id(
  P_ID = x_old |>
    dplyr::filter(has_data == 1) |>
    dplyr::distinct() |>
    dplyr::pull(gene_id),
  conn = genedb_human_old,
  records_filter = "Fetal",
  ages_min = NULL,
  ages_max = NULL,
  unit_filter = "RT-PCR"
)
# Inspect expression profile with: View(expression_values_old)

tab_container_tissue_old <- RSQLite::dbReadTable(
  conn = genedb_human_old,
  name = "tab_container_tissue"
)
expression_profile_old <- dplyr::left_join(
  expression_values_old,
  tab_container_tissue_old
) |>
  dplyr::filter(!is.na(container)) |>
  dplyr::group_by(variant_name, container, unit) |>
  dplyr::mutate(norm_value_var = var(norm_value, y = NULL, na.rm = TRUE)) |>
  dplyr::mutate(norm_value_sd = sd(norm_value, na.rm = TRUE)) |>
  dplyr::mutate(norm_value = mean(norm_value, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::group_by(variant_name, unit) |>
  dplyr::mutate(Rel_Exp = norm_value / max(norm_value, na.rm = TRUE)) |>
  dplyr::filter(unit == "RT-PCR")
# Inspect expression profile with: View(expression_profile_old)

# 2.0 Load Old QSP Expression DB ####
# DB connection
genedb_human_new <- RSQLite::dbConnect(
  drv = RSQLite::SQLite(),
  "./PK-Sim DBs/Human/GENEDB_human_ADME_ONLY_BgeeRelease_15_2.expressionDB"
)

x_new <- get_proteins_by_name(name = enzymes, conn = genedb_human_new)
x_new <- x_new |> dplyr::filter(gene_name %in% enzymes)
# Inspect expression profile with: View(x_new)

expression_values_new <- get_expression_data_by_gene_id(
  P_ID = x_new |>
    dplyr::filter(has_data == 1) |>
    dplyr::distinct() |>
    dplyr::pull(gene_id),
  conn = genedb_human_new,
  records_filter = NULL,
  ages_min = min(expression_values_old$age_min),
  ages_max = max(expression_values_old$age_max),
  unit_filter = NULL
)

expression_values_new <- dplyr::left_join(
  expression_values_new,
  x_new |>
    dplyr::select(gene_id, gene_name) |>
    dplyr::distinct()
)

# Inspect expression profile with: View(expression_values_new)

tab_container_tissue_new <- RSQLite::dbReadTable(
  conn = genedb_human_new,
  name = "tab_container_tissue"
)

expression_profile_new <- dplyr::left_join(
  expression_values_new,
  tab_container_tissue_new |> dplyr::mutate(tissue = tissue)
) |>
  dplyr::filter(!grepl("-", container)) |>
  dplyr::group_by(variant_name, container, unit) |>
  dplyr::mutate(norm_value_var = var(norm_value, y = NULL, na.rm = TRUE)) |>
  dplyr::mutate(norm_value_sd = sd(norm_value, na.rm = TRUE)) |>
  dplyr::mutate(norm_value_max = max(norm_value, na.rm = TRUE)) |>
  dplyr::mutate(norm_value_min = min(norm_value, na.rm = TRUE)) |>
  dplyr::mutate(
    norm_value_geomean =
      #PKNCA::geomean(norm_value, na.rm = TRUE)
      10^mean(log10(norm_value), na.rm = TRUE)
  ) |>
  dplyr::mutate(norm_value_mean = mean(norm_value, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::group_by(variant_name, unit) |>
  dplyr::mutate(Rel_Exp = norm_value_mean / max(norm_value_mean, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::filter(unit == "TPM") |>
  dplyr::distinct()

# Inspect expression profile with: View(expression_profile_new)

# Compare gene expression in a plot genes in QSP-Model-Library
# Create Qualification folder if it doesn't exist
if (!dir.exists("Code/Qualification")) {
  dir.create("Code/Qualification")
}

for (enzyme in enzymes) {
  # Filter for enzyme
  old_data <- expression_profile_old |>
    dplyr::filter(variant_name == enzyme)
  new_data <- expression_profile_new |>
    dplyr::rename(variant_name = gene_name, gene_name = variant_name) |>
    dplyr::filter(variant_name == enzyme)
  # either old or new data is empty skip the itteration
  if (dim(old_data)[1] == 0 || dim(new_data)[1] == 0) {
    next
  }

  # Add source column
  old_data$DB <- "Old"
  new_data$DB <- "New"

  # Combine data
  plot_data <- dplyr::bind_rows(
    old_data |> dplyr::select(container, Rel_Exp, DB) |> dplyr::distinct(),
    new_data |> dplyr::select(container, Rel_Exp, DB) |> dplyr::distinct()
  )

  # Plot
  p1 <- ggplot2::ggplot(
    data = plot_data,
    mapping = ggplot2::aes(x = container, y = Rel_Exp, color = DB)
  ) +
    ggplot2::geom_point(
      size = 3,
      position = ggplot2::position_jitter(width = 0.2, height = 0)
    ) +
    ggplot2::labs(
      title = paste("Relative Expression of", enzyme, ": Old vs New DB"),
      x = "Tissue Container",
      y = "Relative Expression"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    # Prepare plot_data for points (Old) and violin (New)
    plot_data_points <- old_data |> dplyr::select(container, Rel_Exp, DB) |> dplyr::distinct()
    plot_data_violin <- new_data |>
      dplyr::mutate(Rel_Exp = norm_value / max(norm_value_mean)) |>
      #dplyr::mutate(Rel_Exp = ratio / max(norm_value)) |>
      dplyr::select(container, Rel_Exp, DB)

    # Plot: points for Old, violin for New
    p2 <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = plot_data,
        mapping = ggplot2::aes(x = container, y = Rel_Exp, color = DB),
        size = 3,
        position = ggplot2::position_jitter(width = 0.2, height = 0)
      ) +
      ggplot2::geom_violin(
        data = plot_data_violin,
        mapping = ggplot2::aes(x = container, y = Rel_Exp, fill = DB),
        alpha = 0.5,
        width = 0.8
      ) +
      ggplot2::labs(
        title = paste("Relative Expression of", enzyme, ": Old vs New DB"),
        x = "Tissue Container",
        y = "Relative Expression"
      ) +
      ggplot2::scale_y_log10() +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  # Save plot
  ggplot2::ggsave(
    filename = file.path(
      "Code/Qualification",
      paste0("RelExp_", enzyme, "_Old_vs_New.png")
    ),
    plot = p1,
    width = 10,
    height = 6,
    dpi = 300
  )
  ggplot2::ggsave(
    filename = file.path(
      "Code/Qualification",
      paste0("RelExp_", enzyme, "_Old_vs_New_violin.png")
    ),
    plot = p2,
    width = 12,
    height = 6,
    dpi = 300
  )
}
