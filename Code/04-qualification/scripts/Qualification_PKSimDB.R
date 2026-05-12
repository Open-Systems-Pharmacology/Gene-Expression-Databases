# Comparative Qualification: Old vs New OSP Expression Database
# Purpose: Compare gene expression between old (RT-PCR) and new (Bgee-based) DB
# Output: Violin plots by ADME family + validation summary statistics
#
# ENVIRONMENT REQUIREMENTS:
# R version: 4.4.1 (2024-06-14)
# Required packages (with tested versions):
#   - DBI (1.2.3)
#   - RSQLite (2.4.3)
#   - dplyr (1.2.1)
#   - readr (2.1.5)
#   - ggplot2 (4.0.0)
#   - tidyr (1.3.1)
#   - scales (1.4.0)
#   - here (1.0.1)

utils::globalVariables(c("variant_name", "container", "Rel_Exp", "DB", "n", ".x"))

rm(list = ls())
setwd(here::here())

# ============================================================================
# LIBRARY & CONFIG
# ============================================================================
suppressPackageStartupMessages({
  library(here)
  library(DBI)
  library(RSQLite)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

PATH <- here::here()
RELEASE <- "15_2"

# Load helper functions
source(paste0(PATH, "/Code/03-helpers/helper_SQL_Queries.R"))

# Configuration paths
config_dir <- paste0(PATH, "/Code/04-qualification/config")
level_output_dir <- paste0(PATH, "/Code/04-qualification/results/level2")
data_dir <- file.path(level_output_dir, "02_data")
plots_dir <- file.path(level_output_dir, "03_plots")

# Create output directories if they don't exist
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Read protein validation list
proteins_file <- paste0(config_dir, "/proteins-validation.txt")
enzymes <- readLines(proteins_file) |> trimws()

# Read container mapping
mapping_file <- paste0(config_dir, "/container-mapping.txt")
container_mapping_raw <- readLines(mapping_file)
container_mapping <- as.data.frame(
  do.call(rbind, strsplit(container_mapping_raw, "\t")),
  stringsAsFactors = FALSE
)
colnames(container_mapping) <- c("container", "tissue")

# ============================================================================
# LOAD DATABASES
# ============================================================================
message("Loading database connections...")

# Old DB connection - check multiple paths
genedb_old_path <- file.path(PATH, "GENEDB_human.expressionDb")
genedb_old_path_network <- "n:/FraProject/trial/cordesh/Developments/Gene-Expression-Databases/GENEDB_human.expressionDb"

genedb_old <- NULL

if (file.exists(genedb_old_path)) {
  message("  Found Old DB at: ", genedb_old_path)
  genedb_old <- RSQLite::dbConnect(RSQLite::SQLite(), genedb_old_path)
} else if (file.exists(genedb_old_path_network)) {
  message("  Found Old DB at network path: ", genedb_old_path_network)
  genedb_old <- RSQLite::dbConnect(RSQLite::SQLite(), genedb_old_path_network)
} else {
  message("  Old DB not found (checked: ", genedb_old_path, " and network path)")
}

# New DB connection
genedb_new_path <- file.path(PATH, "PK-Sim DBs/Human/GENEDB_human_ADME_ONLY_BgeeRelease_15_2.expressionDB")
if (!file.exists(genedb_new_path)) {
  warning("New DB not found at: ", genedb_new_path)
  genedb_new <- NULL
} else {
  genedb_new <- RSQLite::dbConnect(RSQLite::SQLite(), genedb_new_path)
}

if (is.null(genedb_new)) {
  stop("New DB file is required. Cannot proceed with qualification.")
}

# Determine scenario
SCENARIO_WITH_OLD_DB <- !is.null(genedb_old)

if (!SCENARIO_WITH_OLD_DB) {
  message("WARNING: Old DB not available. Using New DB only for qualification.")
}

# ============================================================================
# LOAD & PROCESS OLD DB
# ============================================================================
if (SCENARIO_WITH_OLD_DB) {
  message("Processing Old DB (fetal, RT-PCR)...")

  x_old <- get_proteins_by_name(name = enzymes, conn = genedb_old)

expression_values_old <- get_expression_data_by_gene_id(
  P_ID = x_old |>
    dplyr::filter(has_data == 1) |>
    dplyr::distinct() |>
    dplyr::pull(gene_id),
  conn = genedb_old,
  records_filter = "Fetal",
  ages_min = NULL,
  ages_max = NULL,
  unit_filter = "RT-PCR"
)

tab_container_tissue_old <- RSQLite::dbReadTable(genedb_old, "tab_container_tissue")

expression_profile_old <- dplyr::left_join(
  expression_values_old,
  tab_container_tissue_old
) |>
  dplyr::filter(!is.na(container)) |>
  dplyr::group_by(variant_name, container, unit) |>
  dplyr::mutate(
    norm_value_var = var(norm_value, na.rm = TRUE),
    norm_value_sd = sd(norm_value, na.rm = TRUE),
    norm_value = mean(norm_value, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |>
  dplyr::group_by(variant_name, unit) |>
  dplyr::mutate(Rel_Exp = norm_value / max(norm_value, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::filter(unit == "RT-PCR")

} else {
  # Old DB unavailable - create empty dataframe
  expression_profile_old <- tibble::tibble(
    variant_name = character(),
    container = character(),
    unit = character(),
    norm_value = numeric(),
    Rel_Exp = numeric()
  )
}

# ============================================================================
# LOAD & PROCESS NEW DB
# ============================================================================
message("Processing New DB (Bgee TPM)...")

x_new <- get_proteins_by_name(name = enzymes, conn = genedb_new)
x_new <- x_new |> dplyr::filter(gene_name %in% enzymes)

# Determine age range for New DB query
if (SCENARIO_WITH_OLD_DB && exists("expression_values_old") && nrow(expression_values_old) > 0) {
  ages_min_new <- min(expression_values_old$age_min, na.rm = TRUE)
  ages_max_new <- max(expression_values_old$age_max, na.rm = TRUE)
} else {
  # Use default ages if old DB not available
  ages_min_new <- NULL
  ages_max_new <- NULL
}

expression_values_new <- get_expression_data_by_gene_id(
  P_ID = x_new |>
    dplyr::filter(has_data == 1) |>
    dplyr::distinct() |>
    dplyr::pull(gene_id),
  conn = genedb_new,
  records_filter = NULL,
  ages_min = ages_min_new,
  ages_max = ages_max_new,
  unit_filter = NULL
)

expression_values_new <- dplyr::left_join(
  expression_values_new,
  x_new |> dplyr::select(gene_id, gene_name) |> dplyr::distinct()
)

tab_container_tissue_new <- RSQLite::dbReadTable(genedb_new, "tab_container_tissue")

expression_profile_new <- dplyr::left_join(
  expression_values_new,
  tab_container_tissue_new |> dplyr::mutate(tissue = tissue)
) |>
  dplyr::filter(!grepl("-", container)) |>
  dplyr::mutate(
    variant_name = dplyr::if_else(!is.na(gene_name) & nzchar(gene_name), gene_name, variant_name)
  ) |>
  dplyr::group_by(variant_name, unit) |>
  dplyr::mutate(Rel_Exp = norm_value / max(norm_value, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::filter(unit == "TPM") |>
  dplyr::filter(is.finite(norm_value), norm_value > 0, is.finite(Rel_Exp), Rel_Exp > 0) |>
  dplyr::distinct()

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Classify gene family based on name
classify_family <- function(gene_name) {
  if (grepl("^CYP", gene_name)) return("CYP")
  if (grepl("^ABC", gene_name)) return("ABC")
  if (grepl("^UGT", gene_name)) return("UGT")
  if (grepl("^SULT", gene_name)) return("SULT")
  if (grepl("^SLC|^OAT|^OATP|^SLCO", gene_name)) return("SLC/Transporter")
  if (grepl("^CES", gene_name)) return("CES")
  return("Other")
}

#' OSP-themed ggplot2 theme
theme_osp <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 13, hjust = 0.5, margin = ggplot2::margin(b = 10)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10, hjust = 0.5, color = "gray40", margin = ggplot2::margin(b = 8)
      ),
      axis.title = ggplot2::element_text(face = "bold", size = 10),
      axis.text = ggplot2::element_text(size = 9, color = "gray20"),
      axis.line = ggplot2::element_line(color = "gray50", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(color = "gray50", linewidth = 0.3),
      strip.text = ggplot2::element_text(
        face = "bold", size = 9, color = "white",
        margin = ggplot2::margin(t = 5, b = 5)
      ),
      strip.background = ggplot2::element_rect(
        fill = "#2E7D9A", color = NA
      ),
      panel.border = ggplot2::element_rect(color = "gray80", fill = NA, linewidth = 0.3),
      panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.2),
      panel.spacing = ggplot2::unit(0.8, "lines"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 9),
      legend.text = ggplot2::element_text(size = 8),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
    )
}

#' Create composite violin plot for a gene family
plot_family_violins <- function(old_data, new_data, family_name, output_path) {
  # Prepare data for plotting
  plot_data_old <- old_data |>
    dplyr::select(variant_name, container, Rel_Exp) |>
    dplyr::distinct() |>
    dplyr::mutate(DB = "Old (RT-PCR)")

  plot_data_new <- new_data |>
    dplyr::select(variant_name, container, Rel_Exp) |>
    dplyr::mutate(DB = "New (Bgee TPM)")

  plot_data <- dplyr::bind_rows(plot_data_old, plot_data_new)

  # Reorder containers logically (if possible)
  container_order <- c(
    "Liver", "Kidney", "Intestine", "Heart", "Brain", "Lung",
    "Muscle", "Fat", "Bone", "Blood", "Skin"
  )
  plot_data <- plot_data |>
    dplyr::mutate(
      container = factor(container, levels = unique(c(
        container_order[container_order %in% unique(plot_data$container)],
        setdiff(unique(plot_data$container), container_order)
      )))
    )

  plot_data_old <- plot_data_old |>
    dplyr::mutate(container = factor(container, levels = levels(plot_data$container)))

  plot_data_new <- plot_data_new |>
    dplyr::mutate(container = factor(container, levels = levels(plot_data$container)))

  violin_data_new <- plot_data_new |>
    dplyr::group_by(variant_name, container) |>
    dplyr::filter(dplyr::n() >= 2) |>
    dplyr::ungroup()

  # Determine number of facets needed
  n_genes <- length(unique(plot_data$variant_name))
  n_facet_cols <- min(3, ceiling(sqrt(n_genes)))
  n_facet_rows <- ceiling(n_genes / n_facet_cols)

  # Create plot: new DB as violin distribution, old DB as comparison points
  p <- ggplot2::ggplot(plot_data_new, ggplot2::aes(x = container, y = Rel_Exp, fill = DB, color = DB)) +
    ggplot2::geom_violin(
      data = violin_data_new,
      position = ggplot2::position_dodge(width = 0.85),
      alpha = 0.15,
      show.legend = TRUE,
      trim = FALSE,
      scale = "width",
      linewidth = 0.2
    ) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.85),
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.4,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = plot_data_new,
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.85),
      size = 1.2,
      alpha = 0.35,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = plot_data_old,
      position = ggplot2::position_jitterdodge(jitter.width = 0.06, dodge.width = 0.85),
      size = 2.6,
      alpha = 0.85,
      shape = 18
    ) +
    ggplot2::scale_y_log10(
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::scale_fill_manual(
      values = c("New (Bgee TPM)" = "#4A90E2", "Old (RT-PCR)" = "#666666"),
      drop = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("Old (RT-PCR)" = "#666666", "New (Bgee TPM)" = "#4A90E2"),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = paste("Human Gene Expression Comparison:", family_name),
      subtitle = "New human data shown as violin distribution, Old human data as points",
      x = "Tissue/Container",
      y = "Relative Expression (log10)",
      fill = "Database",
      color = "Database"
    ) +
    ggplot2::facet_wrap(~variant_name, nrow = n_facet_rows, ncol = n_facet_cols) +
    theme_osp() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom"
    )

  # Save plot
  ggplot2::ggsave(
    plot = p,
    filename = output_path,
    width = 12,
    height = 3 + 3 * n_facet_rows,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  message("Saved: ", output_path)
}

#' Calculate validation statistics
compute_validation_stats <- function(old_data, new_data, gene_name) {
  old_filtered <- as.numeric(old_data |>
    dplyr::filter(variant_name == gene_name) |>
    dplyr::pull(Rel_Exp))

  new_filtered <- as.numeric(new_data |>
    dplyr::filter(variant_name == gene_name) |>
    dplyr::pull(Rel_Exp))

  n_tissues_old <- length(unique(old_data |>
    dplyr::filter(variant_name == gene_name) |>
    dplyr::pull(container)))
  n_tissues_new <- length(unique(new_data |>
    dplyr::filter(variant_name == gene_name) |>
    dplyr::pull(container)))

  # Correlation needs paired values; align by shared containers first
  old_by_container <- old_data |>
    dplyr::filter(variant_name == gene_name) |>
    dplyr::group_by(container) |>
    dplyr::summarise(Rel_Exp_old = mean(Rel_Exp, na.rm = TRUE), .groups = "drop")

  new_by_container <- new_data |>
    dplyr::filter(variant_name == gene_name) |>
    dplyr::group_by(container) |>
    dplyr::summarise(Rel_Exp_new = mean(Rel_Exp, na.rm = TRUE), .groups = "drop")

  paired_values <- dplyr::inner_join(old_by_container, new_by_container, by = "container")

  if (nrow(paired_values) >= 2) {
    cor_val <- stats::cor(
      paired_values$Rel_Exp_old,
      paired_values$Rel_Exp_new,
      use = "pairwise.complete.obs"
    )
  } else {
    cor_val <- NA_real_
  }

  # Threshold 0.7 for container-level paired means correlation:
  # Container aggregation (via old_by_container/new_by_container) computes means per tissue.
  # This container-level approach reduces noise from individual measurements and is more
  # appropriate for gene family validation. The 0.7 threshold reflects high concordance
  # between old and new databases at the container (tissue-level) granularity.
  status <- if (!is.na(cor_val) && cor_val > 0.7) "PASS" else if (is.na(cor_val)) "MISSING" else "WARN"

  tibble::tibble(
    gene_name = gene_name,
    family = classify_family(gene_name),
    n_tissues_old = n_tissues_old,
    n_tissues_new = n_tissues_new,
    mean_rel_exp_old = mean(old_filtered, na.rm = TRUE),
    sd_rel_exp_old = sd(old_filtered, na.rm = TRUE),
    mean_rel_exp_new = mean(new_filtered, na.rm = TRUE),
    sd_rel_exp_new = sd(new_filtered, na.rm = TRUE),
    correlation_old_new = cor_val,
    n_data_points_old = length(old_filtered),
    n_data_points_new = length(new_filtered),
    validation_status = status
  )
}

# ============================================================================
# GENERATE PLOTS BY FAMILY
# ============================================================================
message("Generating plots by gene family...")

families <- c("CYP", "ABC", "UGT", "SULT", "SLC/Transporter", "CES", "Other")

for (family in families) {
  # Filter proteins for this family
  family_genes <- sapply(enzymes, classify_family)
  family_proteins <- enzymes[family_genes == family]

  if (length(family_proteins) == 0) {
    message("  Skipping ", family, " (no genes)")
    next
  }

  message("  Processing ", family, " (", length(family_proteins), " genes)")

  # Filter expression data
  old_family <- expression_profile_old |>
    dplyr::filter(variant_name %in% family_proteins)

  new_family <- expression_profile_new |>
    dplyr::filter(variant_name %in% family_proteins)

  if (nrow(old_family) == 0 && nrow(new_family) == 0) {
    message("    No data found for ", family)
    next
  }

  # Create plot filename
  family_filename <- tolower(gsub("/", "_", family))
  plot_filename <- file.path(plots_dir, paste0("human_old_vs_new_", family_filename, ".png"))

  # Generate plot
  plot_family_violins(old_family, new_family, family, plot_filename)
}

# ============================================================================
# GENERATE VALIDATION SUMMARY
# ============================================================================
message("Computing validation statistics...")

summary_stats <- tibble::tibble()

for (enzyme in enzymes) {
  old_enzyme <- expression_profile_old |> dplyr::filter(variant_name == enzyme)
  new_enzyme <- expression_profile_new |>
    dplyr::filter(variant_name == enzyme)

  if (nrow(old_enzyme) > 0 || nrow(new_enzyme) > 0) {
    stats <- compute_validation_stats(expression_profile_old, expression_profile_new, enzyme)
    summary_stats <- dplyr::bind_rows(summary_stats, stats)
  }
}

# Export summary
summary_file <- file.path(data_dir, "validation_summary.csv")
readr::write_csv(summary_stats, summary_file)
message("Exported validation summary: ", summary_file)

# Print summary and ensure guaranteed cleanup via tryCatch/finally
tryCatch(
  {
    message("\n=== VALIDATION SUMMARY ===")
    print(summary_stats |
      dplyr::select(gene_name, family, validation_status, correlation_old_new))
  },
  finally = {
    # Guaranteed cleanup: always runs regardless of errors
    if (!is.null(genedb_old) && DBI::dbIsValid(genedb_old)) {
      DBI::dbDisconnect(genedb_old)
    }
    if (!is.null(genedb_new) && DBI::dbIsValid(genedb_new)) {
      DBI::dbDisconnect(genedb_new)
    }
    message("\nQualification complete!")
  }
)
