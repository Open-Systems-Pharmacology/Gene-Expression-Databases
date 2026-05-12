# Cross-Species Qualification: 40-gene ADME panel across 5 species
# Purpose: Compare gene expression patterns across human, macaque, dog, rat, mouse
# Output: Individual gene plots (40 total) showing tissue distributions per species
#
# ENVIRONMENT REQUIREMENTS:
# R version: 4.4.1 (2024-06-14)
# Required packages (with tested versions):
#   - DBI (1.2.3)
#   - RSQLite (2.4.3)
#   - dplyr (1.2.1)
#   - readr (2.1.5)
#   - ggplot2 (4.0.0)
#   - scales (1.4.0)
#   - here (1.0.1)

invisible(utils::globalVariables(c(
  "Rel_Exp", "container", "family", "gene_id", "gene_name", "has_data",
  "norm_value", "species_coverage", "species_label", "tissue"
)))

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
  library(scales)
})

PATH <- here::here()
RELEASE <- "15_2"

# Load helper functions
helper_env <- new.env(parent = globalenv())
sys.source(paste0(PATH, "/Code/03-helpers/helper_SQL_Queries.R"), envir = helper_env)
get_proteins_by_name <- helper_env$get_proteins_by_name
get_expression_data_by_gene_id <- helper_env$get_expression_data_by_gene_id

# Configuration paths
config_dir <- paste0(PATH, "/Code/04-qualification/config")
level_output_dir <- paste0(PATH, "/Code/04-qualification/results/level3")
data_dir <- file.path(level_output_dir, "02_data")
cross_species_dir <- file.path(level_output_dir, "03_plots")

# Create output directories
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cross_species_dir, showWarnings = FALSE, recursive = TRUE)

# Read selected genes
genes_file <- paste0(config_dir, "/cross-species-genes.txt")
gene_targets <- readr::read_tsv(
  genes_file,
  col_names = c("gene_name", "family"),
  show_col_types = FALSE,
  progress = FALSE
) |>
  dplyr::filter(!is.na(.data$gene_name), nzchar(.data$gene_name)) |>
  dplyr::mutate(
    gene_name = trimws(.data$gene_name),
    family = trimws(.data$family)
  )

selected_genes <- gene_targets$gene_name

species_config <- tibble::tribble(
  ~species, ~species_label, ~db_path,
  "Human", "Human", file.path(PATH, "PK-Sim DBs/Human", paste0("GENEDB_human_ADME_ONLY_BgeeRelease_", RELEASE, ".expressionDB")),
  "Monkey_mulatta", "Macaque", file.path(PATH, "PK-Sim DBs/Monkey_mulatta", paste0("GENEDB_monkey_mulatta_ADME_ONLY_BgeeRelease_", RELEASE, ".expressionDB")),
  "Dog", "Dog", file.path(PATH, "PK-Sim DBs/Dog", paste0("GENEDB_dog_ADME_ONLY_BgeeRelease_", RELEASE, ".expressionDB")),
  "Rat", "Rat", file.path(PATH, "PK-Sim DBs/Rat", paste0("GENEDB_rat_ADME_ONLY_BgeeRelease_", RELEASE, ".expressionDB")),
  "Mouse", "Mouse", file.path(PATH, "PK-Sim DBs/Mouse", paste0("GENEDB_mouse_ADME_ONLY_BgeeRelease_", RELEASE, ".expressionDB"))
)

species_colors <- c(
  "Human" = "#1F77B4",
  "Macaque" = "#D62728",
  "Dog" = "#2CA02C",
  "Rat" = "#FF7F0E",
  "Mouse" = "#9467BD"
)

preferred_containers <- c(
  "Liver", "Kidney", "SmallIntestine", "LargeIntestine",
  "Lung", "Heart", "Brain", "Spleen"
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

classify_family <- function(gene_name) {
  if (grepl("^CYP", gene_name)) return("CYP")
  if (grepl("^GST", gene_name)) return("GST")
  if (grepl("^UGT", gene_name)) return("UGT")
  if (grepl("^SULT", gene_name)) return("SULT")
  if (grepl("^ABC", gene_name)) return("ABC")
  if (grepl("^SLCO|^OATP", gene_name)) return("SLCO")
  if (grepl("^OAT", gene_name) || grepl("^SLC22A(6|7|8|9|11)$", gene_name)) return("OAT")
  if (grepl("^SLC|^OAT", gene_name)) return("SLC")
  return("Other")
}

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
      panel.border = ggplot2::element_rect(color = "gray80", fill = NA, linewidth = 0.3),
      panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.2),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 9),
      legend.text = ggplot2::element_text(size = 8),
      legend.key = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
    )
}

extract_species_expression <- function(species_name, species_label, db_path, gene_targets) {
  if (!file.exists(db_path)) {
    warning("DB not found for ", species_label, ": ", db_path)
    return(tibble::tibble())
  }

  conn <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(try(RSQLite::dbDisconnect(conn), silent = TRUE), add = TRUE)

  gene_hits <- lapply(seq_len(nrow(gene_targets)), function(index) {
    target <- gene_targets[index, ]

    matches <- get_proteins_by_name(name = target$gene_name, conn = conn) |>
      dplyr::filter(.data$has_data == 1)

    if (nrow(matches) == 0) {
      return(NULL)
    }

    exact_symbol_match <- matches |>
      dplyr::filter(.data$symbol == target$gene_name)

    exact_name_match <- matches |>
      dplyr::filter(.data$gene_name == target$gene_name)

    best_match <- if (nrow(exact_symbol_match) > 0) {
      exact_symbol_match
    } else if (nrow(exact_name_match) > 0) {
      exact_name_match
    } else {
      matches
    }

    best_match |>
      dplyr::slice_head(n = 1) |>
      dplyr::transmute(
        gene_id = .data$gene_id,
        gene_name = target$gene_name,
        family = dplyr::if_else(
          is.na(target$family) | !nzchar(target$family),
          classify_family(target$gene_name),
          target$family
        )
      )
  }) |>
    dplyr::bind_rows() |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()

  if (nrow(gene_hits) == 0) {
    return(tibble::tibble())
  }

  expression_values <- get_expression_data_by_gene_id(
    P_ID = gene_hits$gene_id,
    conn = conn,
    records_filter = NULL,
    ages_min = NULL,
    ages_max = NULL,
    unit_filter = "TPM"
  )

  if (nrow(expression_values) == 0) {
    return(tibble::tibble())
  }

  tissue_map <- unique(RSQLite::dbReadTable(conn, "tab_container_tissue")[, c("tissue", "container")])

  expression_values |>
    dplyr::left_join(gene_hits, by = "gene_id") |>
    dplyr::left_join(tissue_map, by = "tissue", relationship = "many-to-many") |>
    dplyr::mutate(container = dplyr::if_else(is.na(.data$container), .data$tissue, .data$container)) |>
    dplyr::filter(!is.na(.data$container), !grepl("-", .data$container), is.finite(.data$norm_value), .data$norm_value > 0) |>
    dplyr::group_by(.data$gene_name) |>
    dplyr::mutate(Rel_Exp = .data$norm_value / max(.data$norm_value, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      species = species_name,
      species_label = species_label
    )
}

select_containers_for_gene <- function(gene_data) {
  coverage <- gene_data |>
    dplyr::distinct(.data$container, .data$species_label) |>
    dplyr::count(.data$container, name = "species_coverage") |>
    dplyr::arrange(dplyr::desc(.data$species_coverage), .data$container)

  preferred <- coverage |>
    dplyr::filter(.data$container %in% preferred_containers, .data$species_coverage >= 2) |>
    dplyr::pull(.data$container)

  if (length(preferred) >= 3) {
    return(preferred)
  }

  fallback <- coverage |>
    dplyr::filter(.data$species_coverage >= 3) |>
    dplyr::slice_head(n = 8) |>
    dplyr::pull(.data$container)

  if (length(fallback) > 0) {
    return(fallback)
  }

  coverage |>
    dplyr::slice_head(n = min(6, nrow(coverage))) |>
    dplyr::pull(.data$container)
}

plot_cross_species_gene <- function(selected_gene, gene_data, output_path) {
  plot_data <- gene_data |>
    dplyr::filter(.data$gene_name == selected_gene)

  if (nrow(plot_data) == 0) {
    return(invisible(NULL))
  }

  containers_to_plot <- select_containers_for_gene(plot_data)
  plot_data <- plot_data |>
    dplyr::filter(.data$container %in% containers_to_plot)

  if (nrow(plot_data) == 0) {
    return(invisible(NULL))
  }

  coverage_levels <- plot_data |>
    dplyr::distinct(.data$container, .data$species_label) |>
    dplyr::count(.data$container, name = "species_coverage") |>
    dplyr::arrange(dplyr::desc(.data$species_coverage), .data$container) |>
    dplyr::pull(.data$container)

  plot_data <- plot_data |>
    dplyr::mutate(
      container = factor(.data$container, levels = coverage_levels),
      species_label = factor(.data$species_label, levels = names(species_colors))
    )

  violin_data <- plot_data |>
    dplyr::group_by(.data$container, .data$species_label) |>
    dplyr::filter(dplyr::n() >= 2) |>
    dplyr::ungroup()

  gene_family <- unique(plot_data$family)[1]

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$container, y = .data$Rel_Exp, fill = .data$species_label, color = .data$species_label)
  ) +
    ggplot2::geom_violin(
      data = violin_data,
      position = ggplot2::position_dodge(width = 0.85),
      alpha = 0.15,
      trim = FALSE,
      scale = "width",
      linewidth = 0.2
    ) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.85),
      width = 0.16,
      outlier.shape = NA,
      alpha = 0.4,
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitterdodge(jitter.width = 0.12, dodge.width = 0.85),
      size = 1.2,
      alpha = 0.35
    ) +
    ggplot2::scale_y_log10(
      labels = function(value) parse(text = paste0("10^", value))
    ) +
    ggplot2::scale_fill_manual(values = species_colors, drop = FALSE) +
    ggplot2::scale_color_manual(values = species_colors, drop = FALSE) +
    ggplot2::labs(
      title = paste("Cross-Species Expression Comparison:", selected_gene),
      subtitle = paste(gene_family, "family | New Bgee-based PK-Sim DBs | Human, Macaque, Dog, Rat, Mouse"),
      x = "Tissue/Container",
      y = "Relative Expression (log10)",
      fill = "Species",
      color = "Species"
    ) +
    theme_osp() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom"
    )

  ggplot2::ggsave(
    filename = output_path,
    plot = p,
    width = 12,
    height = 7,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  message("Saved: ", output_path)
  invisible(output_path)
}

# ============================================================================
# LOAD & PROCESS SPECIES DATA
# ============================================================================
message("Loading cross-species expression data...")

species_expression <- lapply(seq_len(nrow(species_config)), function(index) {
  cfg <- species_config[index, ]
  message("  Processing ", cfg$species_label, "...")
  extract_species_expression(
    species_name = cfg$species,
    species_label = cfg$species_label,
    db_path = cfg$db_path,
    gene_targets = gene_targets
  )
})

cross_species_data <- dplyr::bind_rows(species_expression)

if (nrow(cross_species_data) == 0) {
  stop("No cross-species data available for the selected genes.")
}

cross_species_data <- cross_species_data |>
  dplyr::mutate(species_label = factor(species_label, levels = names(species_colors)))

# Export sample-level data
sample_level_file <- file.path(data_dir, "cross_species_selected_genes.csv")
readr::write_csv(cross_species_data, sample_level_file)
message("Exported sample-level cross-species data: ", sample_level_file)

# Export summary table
summary_file <- file.path(data_dir, "cross_species_summary.csv")
cross_species_summary <- cross_species_data |>
  dplyr::group_by(.data$gene_name, .data$family, .data$species_label, .data$container) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    mean_rel_exp = mean(.data$Rel_Exp, na.rm = TRUE),
    median_rel_exp = median(.data$Rel_Exp, na.rm = TRUE),
    max_rel_exp = max(.data$Rel_Exp, na.rm = TRUE),
    .groups = "drop"
  )
readr::write_csv(cross_species_summary, summary_file)
message("Exported cross-species summary: ", summary_file)

# ============================================================================
# GENERATE GENE-SPECIFIC PLOTS
# ============================================================================
message("Generating cross-species plots...")

for (gene in selected_genes) {
  gene_data <- cross_species_data |>
    dplyr::filter(.data$gene_name == gene)

  if (nrow(gene_data) == 0) {
    message("  Skipping ", gene, " (no data across selected species)")
    next
  }

  plot_path <- file.path(cross_species_dir, paste0("cross_species_", tolower(gene), ".png"))
  plot_cross_species_gene(gene, cross_species_data, plot_path)
}

message("\nCross-species qualification complete!")