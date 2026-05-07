# Technical Qualification: BgeeDB to PKSimDB Integration
# Purpose: Validate raw data integrity when integrating Bgee source → PKSimDB
# Output: Scatter plot with technical validation metrics by ADME family

rm(list = ls())
setwd(here::here())

# ============================================================================
# LIBRARY & CONFIG
# ============================================================================
suppressPackageStartupMessages({
  library(here)
  library(BgeeDB)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

PATH <- here::here()
RELEASE <- "15_2"

# Configuration paths
config_dir <- paste0(PATH, "/Code/04-qualification/Qualification/01_config")
data_dir <- paste0(PATH, "/Code/04-qualification/Qualification/02_data")
plots_dir <- paste0(PATH, "/Code/04-qualification/Qualification/03_plots")

# Create output directories
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Read protein validation list
proteins_file <- paste0(config_dir, "/proteins-validation.txt")
proteins4validation_list <- readLines(proteins_file) |> trimws()

# ============================================================================
# FETCH & CACHE BGEE DATA
# ============================================================================
message("Loading Bgee data (GSE30611, ERX011211)...")

bgee_cache_path <- file.path(data_dir, "DATA_HUMAN_GSE30611_ERX011211_BgeeDB.csv")

if (!file.exists(bgee_cache_path)) {
  message("  Fetching from Bgee API...")
  bgee <- BgeeDB::Bgee$new(
    species = "Homo_sapiens",
    dataType = "rna_seq",
    sendStats = FALSE,
    release = "15_2"
  )

  GSE30611_ERX011211_BgeeDB <- BgeeDB::getSampleProcessedData(
    myBegeeObject = bgee,
    experimentId = "GSE30611",
    sampleId = "ERX011211"
  )

  readr::write_csv(GSE30611_ERX011211_BgeeDB, bgee_cache_path)
  message("  Cached to: ", bgee_cache_path)
} else {
  message("  Loading from cache: ", bgee_cache_path)
  GSE30611_ERX011211_BgeeDB <- readr::read_csv(bgee_cache_path, show_col_types = FALSE)
}

# ============================================================================
# LOAD OSP DB DATA
# ============================================================================
message("Loading OSP DB data...")

osp_cache_path <- file.path(data_dir, "DATA_HUMAN_GSE30611_ERX011211_OSP_DB.csv")

if (!file.exists(osp_cache_path)) {
  message("  Querying OSP DB...")
  # Note: This assumes tab_expression_data_records and related tables are loaded
  # from a previously opened database connection (from GeneratePKsimDB or similar)
  message("  WARNING: OSP DB tables not available. Skipping OSP data extraction.")
  GSE30611_ERX011211_OSP_DB <- NULL
} else {
  message("  Loading from cache: ", osp_cache_path)
  GSE30611_ERX011211_OSP_DB <- readr::read_csv(osp_cache_path, show_col_types = FALSE)
}

# ============================================================================
# PREPARE DATA & CREATE X-Y PLOTS BY GENE FAMILY
# ============================================================================
message("Generating technical validation plots by gene family...")

# Classify genes and add family information to Bgee data
bgee_classified <- GSE30611_ERX011211_BgeeDB |>
  dplyr::filter(!is.na(TPM), TPM > 0) |>
  dplyr::mutate(
    Gene_Family = dplyr::case_when(
      grepl("^CYP", Gene.ID) ~ "CYP",
      grepl("^ABC", Gene.ID) ~ "ABC",
      grepl("^UGT", Gene.ID) ~ "UGT",
      grepl("^SULT", Gene.ID) ~ "SULT",
      grepl("^SLC|^OAT|^OATP|^SLCO", Gene.ID) ~ "SLC/Transporter",
      grepl("^CES", Gene.ID) ~ "CES",
      TRUE ~ "Other"
    )
  )

# Define families and colors for consistent plotting
families_to_plot <- c("CYP", "ABC", "UGT", "SULT", "SLC/Transporter", "CES", "Other")
family_colors <- c(
  "CYP" = "#E74C3C", "ABC" = "#3498DB", "UGT" = "#2ECC71",
  "SULT" = "#F39C12", "SLC/Transporter" = "#9B59B6", "CES" = "#1ABC9C", "Other" = "#95A5A6"
)

# ============================================================================
# SCENARIO 1: OSP DB Available - Create X-Y comparison plots
# ============================================================================
if (!is.null(GSE30611_ERX011211_OSP_DB)) {
  message("  OSP DB available: Creating Bgee → OSP comparison plots")

  # Join Bgee and OSP data
  joint_data <- dplyr::full_join(
    GSE30611_ERX011211_BgeeDB |>
      dplyr::filter(Detection.flag == "present", !is.na(TPM), TPM > 0) |>
      dplyr::select(Gene.ID, TPM),
    GSE30611_ERX011211_OSP_DB |>
      dplyr::select(variant_name, sample_count, gene_name),
    by = dplyr::join_by(Gene.ID == variant_name)
  ) |>
    dplyr::filter(!is.na(TPM), !is.na(sample_count)) |>
    dplyr::mutate(
      Gene_Family = dplyr::case_when(
        grepl("^CYP", gene_name) ~ "CYP",
        grepl("^ABC", gene_name) ~ "ABC",
        grepl("^UGT", gene_name) ~ "UGT",
        grepl("^SULT", gene_name) ~ "SULT",
        grepl("^SLC|^OAT|^OATP|^SLCO", gene_name) ~ "SLC/Transporter",
        grepl("^CES", gene_name) ~ "CES",
        TRUE ~ "Other"
      )
    )

  # Generate X-Y scatter plots for each family
  for (family in families_to_plot) {
    family_data <- joint_data |> dplyr::filter(Gene_Family == family)

    if (nrow(family_data) == 0) {
      message("    Skipping ", family, " (no data)")
      next
    }

    message("    Plotting ", family, " (", nrow(family_data), " genes)")

    p <- ggplot2::ggplot(family_data, ggplot2::aes(x = TPM, y = sample_count)) +
      ggplot2::geom_point(
        color = family_colors[family],
        size = 3.5,
        alpha = 0.6,
        stroke = 0.5,
        shape = 21,
        fill = family_colors[family]
      ) +
       ggrepel::geom_text_repel(
         ggplot2::aes(label = gene_name),
         size = 2.5,
         max.overlaps = 20,
         alpha = 0.7,
         box.padding = ggplot2::unit(0.3, "lines"),
         point.padding = ggplot2::unit(0.3, "lines")
       ) +
      ggplot2::geom_smooth(
        method = "lm",
        se = FALSE,
        color = "gray40",
        linetype = "dashed",
         linewidth = 0.6,
        alpha = 0.5
      ) +
      ggplot2::scale_x_log10(
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      ggplot2::scale_y_log10(
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      ggplot2::labs(
        title = paste("Technical Validation:", family),
        subtitle = "Bgee TPM vs OSP DB Sample Count | GSE30611 ERX011211",
        x = "Bgee TPM (log10)",
        y = "OSP DB Sample Count (log10)"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, color = "gray40", hjust = 0.5),
        axis.title = ggplot2::element_text(face = "bold", size = 10),
        axis.text = ggplot2::element_text(size = 9),
        panel.grid.major = ggplot2::element_line(color = "gray90", size = 0.2),
        panel.grid.minor = ggplot2::element_line(color = "gray95", size = 0.1),
        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
      )

    family_filename <- tolower(gsub("/", "_", family))
    plot_path <- file.path(plots_dir, paste0("technical_validation_xy_", family_filename, ".png"))
    ggplot2::ggsave(plot_path, p, width = 8, height = 7, dpi = 300, bg = "white")
    message("      Saved: ", plot_path)
  }

  # Export comparison data
  comparison_file <- file.path(data_dir, "technical_validation_data.csv")
  readr::write_csv(joint_data, comparison_file)
  message("Exported comparison data: ", comparison_file)

} else {
  # ============================================================================
  # SCENARIO 2: OSP DB Unavailable - Create Bgee overview plots with all samples
  # ============================================================================
  message("  OSP DB unavailable: Creating Bgee-focused overview plots")

  # Generate distribution plots by family
  for (family in families_to_plot) {
    family_data <- bgee_classified |> dplyr::filter(Gene_Family == family)

    if (nrow(family_data) == 0) {
      message("    Skipping ", family, " (no data)")
      next
    }

    message("    Plotting ", family, " (", length(unique(family_data$Gene.ID)), " genes)")

    # Create scatter plot: TPM vs Detection Flag colored by gene
    p <- ggplot2::ggplot(family_data, ggplot2::aes(x = TPM, y = Detection.flag)) +
      ggplot2::geom_jitter(
        height = 0.15,
        width = 0,
        size = 2.5,
        alpha = 0.5,
        color = family_colors[family],
        shape = 21,
        fill = family_colors[family]
      ) +
      ggplot2::scale_x_log10(
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      ggplot2::labs(
        title = paste("Technical Overview:", family),
        subtitle = "Bgee GSE30611 ERX011211 | All Detection Flags",
        x = "TPM (log10)",
        y = "Detection Flag"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, color = "gray40", hjust = 0.5),
        axis.title = ggplot2::element_text(face = "bold", size = 10),
        axis.text = ggplot2::element_text(size = 9),
        panel.grid.major = ggplot2::element_line(color = "gray90", size = 0.2),
        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
      )

    family_filename <- tolower(gsub("/", "_", family))
    plot_path <- file.path(plots_dir, paste0("technical_bgee_overview_", family_filename, ".png"))
    ggplot2::ggsave(plot_path, p, width = 8, height = 6, dpi = 300, bg = "white")
    message("      Saved: ", plot_path)
  }

  # Export Bgee overview data
  overview_file <- file.path(data_dir, "technical_bgee_overview.csv")
  readr::write_csv(bgee_classified, overview_file)
  message("Exported Bgee overview: ", overview_file)
}

message("\nTechnical qualification complete!")
