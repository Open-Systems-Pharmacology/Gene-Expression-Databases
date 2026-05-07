# Technical qualification BgeeDB 2 OSP_DB #####
# Goal: Ensure that data in BgeeDB and OSP_DB are identical
# Task: Choose example study e.g. ERP003613 and plot data from both sources

# free workspace
rm(list = ls())

# set working and export directory
setwd(here::here())

proteins4validation <- xlsx::read.xlsx(
  file = "Code/Qualification/Proteins4Validation.xlsx",
  sheetIndex = 1, header = FALSE
) |>
  dplyr::rename("SYMBOL" = "X1")
proteins4validation <- proteins4validation |>
  dplyr::mutate(SYMBOL = trimws(SYMBOL, "both"))

# 1.0 Load original Bgee data ####
# List available data: BgeeDB::listBgeeRelease()
# List available species: BgeeDB::listBgeeSpecies(release = "15")

# Fetch data form DB
# if file exist skip extraction just load
if (!file.exists("Code/Qualification/DATA_HUMAN_GSE30611_ERX011211_BgeeDB.csv")) {
  bgee <- BgeeDB::Bgee$new(
    species = "Homo_sapiens",
    dataType = "rna_seq",
    sendStats = FALSE,
    release = "15_2"
  )

  GSE30611_ERX011211_BgeeDB <-
    BgeeDB::getSampleProcessedData(bgee,
      myBegeeObject = bgee,
      experimentId = "GSE30611",
      sampleId = "ERX011211"
    )
  # Export
  readr::write_csv(
    x = GSE30611_ERX011211_BgeeDB,
    file = "Code/Qualification/DATA_HUMAN_GSE30611_ERX011211_BgeeDB.csv"
  )
} else {
  GSE30611_ERX011211_BgeeDB <-
    readr::read_csv(
      file = "Code/Qualification/DATA_HUMAN_GSE30611_ERX011211_BgeeDB.csv"
    )
}

# 1.1 Load OSP expression DB ####

# if DATA_HUMAN_GSE30611_ERX011211_OSP_DB.csv file exist skip extraction
if (!file.exists("Code/Qualification/DATA_HUMAN_GSE30611_ERX011211_OSP_DB.csv")) {
  GSE30611_ERX011211_OSP_DB <- tab_expression_data_records |>
    dplyr::left_join(., tab_expression_data_values) |>
    dplyr::left_join(., tab_expression_data_properties |>
      dplyr::filter(property == "data_base") |>
      dplyr::select(-property)) |>
    dplyr::filter(unit == "TPM") |>
    dplyr::left_join(., tab_gene_variants) |>
    dplyr::left_join(., tab_gene_names |>
      dplyr::filter(name_type %in% c("SYMBOL")) |>
      dplyr::select(-name_type)) |>
    dplyr::collect()

  readr::write_csv(
    x = GSE30611_ERX011211_OSP_DB,
    file = "Code/Qualification/DATA_HUMAN_GSE30611_ERX011211_OSP_DB.csv"
  )
} else {
  GSE30611_ERX011211_OSP_DB <-
    readr::read_csv(
      file = "Code/Qualification/DATA_HUMAN_GSE30611_ERX011211_OSP_DB.csv"
    )
}

# Define plot theme
# Unified plot function for all validation plots
make_validation_plot <- function(data, title = NULL) {
  ggplot2::ggplot(data = data) +
    ggplot2::geom_smooth(ggplot2::aes(x = TPM, y = sample_count),
      method = "lm", se = FALSE, color = "black"
    ) +
    ggplot2::geom_point(ggplot2::aes(x = TPM, y = sample_count)) +
    ggplot2::scale_x_log10(
      limits = c(0.01, 10000),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::scale_y_log10(
      limits = c(0.01, 10000),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::xlab("BgeeDB Source TPM") +
    ggplot2::ylab("OSP suite DB sample count (TPM)") +
    ggrepel::geom_label_repel(
      ggplot2::aes(x = TPM, y = sample_count, label = gene_name),
      box.padding = attr(theme_quali_plots, "label_repel")$box.padding,
      point.padding = attr(theme_quali_plots, "label_repel")$point.padding,
      force = attr(theme_quali_plots, "label_repel")$force,
      segment.color = attr(theme_quali_plots, "label_repel")$segment.color,
      segment.size = attr(theme_quali_plots, "label_repel")$segment.size,
      max.overlaps = attr(theme_quali_plots, "label_repel")$max.overlaps,
      size = attr(theme_quali_plots, "label_repel")$size
    ) +
    (if (!is.null(title)) ggplot2::ggtitle(label = title) else NULL) +
    theme_quali_plots
}
# Square plot theme: width = height = 7.4 cm (half A5 width)
theme_quali_plots <- ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    # Axis properties
    axis.title = ggplot2::element_text(face = "bold", size = 8),
    axis.text = ggplot2::element_text(color = "black", size = 7),
    axis.line = ggplot2::element_line(size = 0.5),
    axis.ticks = ggplot2::element_line(size = 0.5),
    # Legend properties
    legend.position = "bottom",
    legend.title = ggplot2::element_text(face = "bold", size = 8),
    legend.text = ggplot2::element_text(size = 7),
    legend.key.size = grid::unit(0.7, "lines"),
    # Title properties
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 9),
    plot.subtitle = ggplot2::element_text(size = 8),
    plot.caption = ggplot2::element_text(size = 7),
    # Facet/strip properties
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 7),
    # Panel/border properties
    panel.border =
      ggplot2::element_rect(color = "black", fill = NA, size = 0.7),
    panel.grid.major =
      ggplot2::element_line(size = 0.3, linetype = "dotted", color = "grey90"),
    panel.grid.minor =
      ggplot2::element_line(size = 0.2, linetype = "dotted", color = "grey95"),
    # ggrepel label text size
    label.text = ggplot2::element_text(size = 3)
  )

attr(theme_quali_plots, "label_repel") <- list(
  box.padding = 0.15,
  point.padding = 0.15,
  force = 30,
  segment.color = "grey50",
  segment.size = 0.15,
  max.overlaps = Inf,
  size = 2
)
# set size to half A4 page
attr(theme_quali_plots, "ggsave_size") <- list(
  width = 10.5,
  height = 10.5,
  units = "cm"
)
# Define reference genes to validate (ADME focus)
save_validation_plot <- function(plot, filename) {
  ggplot2::ggsave(
    plot = plot,
    filename = filename,
    device = "png",
    width = attr(theme_quali_plots, "ggsave_size")$width,
    height = attr(theme_quali_plots, "ggsave_size")$height,
    units = attr(theme_quali_plots, "ggsave_size")$units
  )
}

#### 1.2 compare Bgee DB based and new OSP expression DB values ####
joint_bgeedb_osp_db <- dplyr::full_join(
  GSE30611_ERX011211_BgeeDB,
  GSE30611_ERX011211_OSP_DB |>
    dplyr::select(
      property_value, data_base_rec_id,
      variant_name, tissue, gender,
      sample_count, gene_name
    ) |>
    dplyr::distinct(),
  by = dplyr::join_by(
    Experiment.ID == property_value,
    Library.ID == data_base_rec_id,
    Gene.ID == variant_name,
    Anatomical.entity.name == tissue,
    Sex == gender
  ),
  keep = FALSE
)

joint_bgeedb_osp_db <- dplyr::full_join(
  GSE30611_ERX011211_OSP_DB,
  GSE30611_ERX011211_BgeeDB |>
    dplyr::filter(Detection.flag == "present") |>
    dplyr::select(Gene.ID, TPM),
  by = dplyr::join_by(variant_name == Gene.ID)
)

## All data ####
data_val <- joint_bgeedb_osp_db |>
  dplyr::filter(gene_name %in% Proteins4Validation$SYMBOL)

## All CYP data ####
cyp_data <- joint_bgeedb_osp_db |>
  dplyr::filter(grepl("^CYP", gene_name))

## ABC transporter data ####
abc_data <- joint_bgeedb_osp_db |>
  dplyr::filter(grepl("^ABC", gene_name))

## OAT transporter data ####
ces_data <- joint_bgeedb_osp_db |>
  dplyr::filter(grepl("^CES", gene_name))

## SLC transporter data ####
slc_data <- joint_bgeedb_osp_db |>
  dplyr::filter(grepl("^SLC", gene_name))

## UGT data ####
ugt_data <- joint_bgeedb_osp_db |>
  dplyr::filter(grepl("^UGT", gene_name))

## SULT data ####
sult_data <- joint_bgeedb_osp_db |>
  dplyr::filter(grepl("^SULT", gene_name))

## plot all data
plot_list <- list(
  list(
    data = data_val,
    filename = "./Code/Qualification/Technical_OSP_Library_Processes.png",
    title = "HUMAN GSE30611 ERX011211"
  ),
  list(
    data = cyp_data,
    filename = "./Code/Qualification/Technical_validation_CYP.png",
    title = "HUMAN GSE30611 ERX011211"
  ),
  list(
    data = abc_data,
    filename = "./Code/Qualification/Technical_validation_ABC.png",
    title = "HUMAN GSE30611 ERX011211"
  ),
  list(
    data = ces_data,
    filename = "./Code/Qualification/Technical_validation_CES.png",
    title = "HUMAN GSE30611 ERX011211"
  ),
  list(
    data = slc_data,
    filename = "./Code/Qualification/Technical_validation_SLC.png",
    title = "HUMAN GSE30611 ERX011211"
  ),
  list(
    data = ugt_data,
    filename = "./Code/Qualification/Technical_validation_UGT.png",
    title = "HUMAN GSE30611 ERX011211"
  ),
  list(
    data = sult_data,
    filename = "./Code/Qualification/Technical_validation_SULT.png",
    title = "HUMAN GSE30611 ERX011211"
  )
)

# save plots
for (plt in plot_list) {
  p <- make_validation_plot(plt$data, title = plt$title)
  save_validation_plot(p, plt$filename)
}
