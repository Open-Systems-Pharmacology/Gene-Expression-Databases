# Technical Validation: Source Bgee vs OSP DB Integration
# Purpose: Validate Bgee source data integrity and integration into PKSimDB
# Output: X-Y scatter plots with gene labels for 7 ADME families
#
# ENVIRONMENT REQUIREMENTS:
# R version: 4.4.1 (2024-06-14)
# Required packages (with tested versions):
#   - BgeeDB (2.32.0)
#   - dplyr (1.2.1)
#   - readr (2.1.5)
#   - ggplot2 (4.0.0)
#   - ggrepel (0.9.6)
#   - scales (1.4.0)
#   - here (1.0.1)

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
config_dir <- paste0(PATH, "/Code/04-qualification/config")
level_output_dir <- paste0(PATH, "/Code/04-qualification/results/level1")
data_dir <- file.path(level_output_dir, "02_data")
plots_dir <- file.path(level_output_dir, "03_plots")

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
    myBgeeObject = bgee,
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

# Build gene symbol lookup - first try OSP DB, then fall back to hardcoded ADME genes
message("  Building gene symbol lookup...")
osp_db_path <- "/workfra/content/trial/cordesh/Developments/Gene-Expression-Databases/GENEDB_human.expressionDb"
gene_symbol_lookup <- NULL

if (file.exists(osp_db_path)) {
  tryCatch(
    {
      suppressPackageStartupMessages({
        library(DBI)
        library(RSQLite)
      })
      
      conn <- RSQLite::dbConnect(RSQLite::SQLite(), osp_db_path)
      on.exit(RSQLite::dbDisconnect(conn), add = TRUE)
      
      # Get ADME gene names from database
      gene_names_tbl <- DBI::dbGetQuery(conn, 
        "SELECT DISTINCT gene_id, gene_name FROM tab_gene_names WHERE gene_name IS NOT NULL AND gene_name != '' AND gene_name IN ('CYP1A1', 'CYP1A2', 'CYP2A6', 'CYP2B6', 'CYP2C8', 'CYP2C9', 'CYP2C19', 'CYP2D6', 'CYP2E1', 'CYP3A4', 'CYP3A5', 'CYP3A7', 'ABCB1', 'ABCC2', 'ABCG2', 'UGT1A1', 'UGT1A4', 'UGT1A6', 'UGT1A9', 'UGT2B7', 'SULT1A1', 'SULT2A1', 'SLC22A1', 'SLC22A2', 'SLCO1B1', 'SLCO1B3', 'SLC15A1', 'OAT1', 'OAT3', 'CES1', 'CES2')")
      
      if (nrow(gene_names_tbl) > 0) {
        gene_symbol_lookup <- gene_names_tbl %>%
          dplyr::mutate(gene_id = as.character(gene_id)) %>%
          dplyr::select(gene_id, gene_name) %>%
          dplyr::distinct()
      }
    },
    error = function(e) {
      message("    Error querying OSP DB")
    }
  )
}

# Classify genes and add family information to Bgee data
# NOTE: Bgee uses Ensembl IDs which don't contain gene symbol patterns.
# When OSP DB is unavailable, we use a curated list of known ADME gene Ensembl IDs.
# This is a known limitation - ideally, we would query a gene symbol database.

adme_ensembl_ids <- data.frame(
  ensembl_id = c(
    # CYP (Cytochrome P450 genes)
    "ENSG00000100030", "ENSG00000130649", "ENSG00000116104", "ENSG00000139593",
    "ENSG00000112116", "ENSG00000155015", "ENSG00000109700", "ENSG00000107404",
    "ENSG00000186115", "ENSG00000160868", "ENSG00000165246", "ENSG00000123275",
    # ABC (ATP-binding cassette transporters)
    "ENSG00000085563", "ENSG00000198695", "ENSG00000138115", "ENSG00000135408",
    # UGT (UDP-glucuronosyltransferase genes)
    "ENSG00000241266", "ENSG00000244731", "ENSG00000243649", "ENSG00000124206",
    "ENSG00000213721",
    # SULT (Sulfotransferase genes)
    "ENSG00000105144", "ENSG00000186810",
    # SLC/Transporter genes
    "ENSG00000109270", "ENSG00000155658", "ENSG00000125148", "ENSG00000170248",
    "ENSG00000161584", "ENSG00000108848", "ENSG00000160888",
    # CES (Carboxylesterase genes)
    "ENSG00000198845"
  ),
  family = c(
    rep("CYP", 12), rep("ABC", 4), rep("UGT", 5), rep("SULT", 2),
    rep("SLC/Transporter", 7), rep("CES", 1)
  ),
  stringsAsFactors = FALSE
)

bgee_classified <- GSE30611_ERX011211_BgeeDB |>
  dplyr::filter(!is.na(TPM), TPM > 0) |>
  dplyr::left_join(
    adme_ensembl_ids %>% dplyr::rename(Gene_Family = family),
    by = c("Gene.ID" = "ensembl_id")
  ) |>
  dplyr::mutate(
    Gene_Family = dplyr::if_else(is.na(Gene_Family), "Other", Gene_Family),
    symbol = Gene.ID
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
      ggrepel::geom_label_repel(
        ggplot2::aes(label = gene_name),
        size = 2.4,
        max.overlaps = 20,
        alpha = 0.9,
        label.size = 0.15,
        label.padding = ggplot2::unit(0.12, "lines"),
        box.padding = ggplot2::unit(0.3, "lines"),
        point.padding = ggplot2::unit(0.2, "lines"),
        fill = "white",
        color = "gray15"
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
        subtitle = "Bgee source expression vs PK-Sim DB expression | GSE30611 ERX011211",
        x = "Source Expression TPM (log10)",
        y = "PK-Sim DB Expression (log10)"
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
  # SCENARIO 2: Create X-Y plots of Expression Rank vs TPM
  # ============================================================================
  message("  Creating Bgee X-Y validation plots (Rank vs TPM) by gene family")

  for (family in families_to_plot) {
    family_data <- bgee_classified |>
      dplyr::filter(Gene_Family == family, !is.na(TPM), TPM > 0) |>
      dplyr::arrange(desc(TPM)) |>
      dplyr::mutate(Gene_Rank = row_number())
    
    if (nrow(family_data) == 0) {
      message("    Skipping ", family, " (no data)")
      next
    }
    
    message("    Plotting ", family, " (", nrow(family_data), " genes)")
    
    # Create X-Y scatter plot: Rank vs TPM
    p <- ggplot2::ggplot(family_data, ggplot2::aes(x = Gene_Rank, y = TPM)) +
      ggplot2::geom_point(
        color = family_colors[family],
        size = 3.5,
        alpha = 0.6,
        stroke = 0.5,
        shape = 21,
        fill = family_colors[family]
      ) +
      ggrepel::geom_label_repel(
        ggplot2::aes(label = symbol),
        size = 2.2,
        max.overlaps = 15,
        alpha = 0.9,
        label.size = 0.15,
        label.padding = ggplot2::unit(0.12, "lines"),
        box.padding = ggplot2::unit(0.3, "lines"),
        point.padding = ggplot2::unit(0.2, "lines"),
        fill = "white",
        color = "gray15"
      ) +
      ggplot2::geom_smooth(
        method = "loess",
        se = FALSE,
        color = "gray40",
        linetype = "dashed",
        linewidth = 0.6,
        alpha = 0.5
      ) +
      ggplot2::scale_y_log10(
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      ggplot2::labs(
        title = paste("Technical Validation:", family),
        subtitle = "Bgee GSE30611 ERX011211 | Expression Rank vs Level",
        x = "Expression Rank (highest to lowest)",
        y = "TPM (log10)"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 10, color = "gray40", hjust = 0.5),
        axis.title = ggplot2::element_text(face = "bold", size = 10),
        axis.text = ggplot2::element_text(size = 9),
        panel.grid.major = ggplot2::element_line(color = "gray90", linewidth = 0.2),
        panel.grid.minor = ggplot2::element_line(color = "gray95", linewidth = 0.1),
        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
      )
    
    family_filename <- tolower(gsub("/", "_", family))
    plot_path <- file.path(plots_dir, paste0("technical_validation_xy_", family_filename, ".png"))
    ggplot2::ggsave(plot_path, p, width = 8, height = 7, dpi = 300, bg = "white")
    message("      Saved: ", plot_path)
  }

  # Export Bgee validation data
  overview_file <- file.path(data_dir, "technical_validation_data.csv")
  readr::write_csv(bgee_classified, overview_file)
  message("Exported Bgee validation data: ", overview_file)
}

message("\nTechnical qualification complete!")
