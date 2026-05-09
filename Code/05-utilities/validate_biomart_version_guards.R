#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(biomaRt))

entries <- list(
  list(species = "Cat", dataset = "fcatus_gene_ensembl", hosts = c("https://may2021.archive.ensembl.org", "https://www.ensembl.org"), expected = "^Felis_catus_9\\.0$"),
  list(species = "Cattle", dataset = "btaurus_gene_ensembl", hosts = c("https://may2021.archive.ensembl.org", "https://www.ensembl.org"), expected = "^ARS-UCD1\\.2$"),
  list(species = "Chicken", dataset = "ggallus_gene_ensembl", hosts = c("https://apr2022.archive.ensembl.org", "https://may2021.archive.ensembl.org"), expected = "^GRCg6a$"),
  list(species = "Dog", dataset = "clfamiliaris_gene_ensembl", hosts = c("https://may2021.archive.ensembl.org", "https://www.ensembl.org"), expected = "^CanFam3\\.1$"),
  list(species = "Goat", dataset = "chircus_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^ARS1$"),
  list(species = "Guineapig", dataset = "cporcellus_gene_ensembl", hosts = c("https://may2025.archive.ensembl.org", "https://www.ensembl.org"), expected = "^Cavpor3\\.0$"),
  list(species = "Horse", dataset = "ecaballus_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^EquCab3\\.0$"),
  list(species = "Human", dataset = "hsapiens_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^GRCh38\\.p14$"),
  list(species = "Minipig", dataset = "sscrofa_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^Sscrofa11\\.1$"),
  list(species = "Monkey_fascicularis", dataset = "mfascicularis_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^Macaca_fascicularis_6\\.0$"),
  list(species = "Monkey_mulatta", dataset = "mmulatta_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^Mmul_10$"),
  list(species = "Monkey_PigTailed", dataset = "mnemestrina_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^Mnem_1\\.0$"),
  list(species = "Mouse", dataset = "mmusculus_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^GRCm39$"),
  list(species = "Rabbit", dataset = "ocuniculus_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^OryCun2\\.0$"),
  list(species = "Rat", dataset = "rnorvegicus_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^GRCr8$"),
  list(species = "Sheep", dataset = "oaries_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^ARS-UI_Ramb_v3\\.0$"),
  list(species = "Turkey", dataset = "mgallopavo_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^Turkey_5\\.1$"),
  list(species = "Zebrafish", dataset = "drerio_gene_ensembl", hosts = c("https://oct2024.archive.ensembl.org", "https://www.ensembl.org"), expected = "^GRCz11$")
)

check_entry <- function(e) {
  attempts <- c()
  for (h in e$hosts) {
    host_label <- ifelse(is.na(h), "<default>", h)
    res <- tryCatch({
      mart <- if (is.na(h)) {
        biomaRt::useMart("ensembl", dataset = e$dataset)
      } else {
        biomaRt::useMart("ensembl", dataset = e$dataset, host = h)
      }
      meta <- if (is.na(h)) {
        biomaRt::useMart("ensembl")
      } else {
        biomaRt::useMart("ensembl", host = h)
      }
      ds <- biomaRt::listDatasets(meta)
      row <- ds[ds$dataset == e$dataset, , drop = FALSE]
      if (nrow(row) == 0) {
        list(ok = FALSE, host = host_label, version = "", msg = "dataset_not_listed")
      } else {
        v <- as.character(row$version[1])
        if (grepl(e$expected, v, ignore.case = TRUE)) {
          list(ok = TRUE, host = host_label, version = v, msg = "version_match")
        } else {
          list(ok = FALSE, host = host_label, version = v, msg = paste0("version_mismatch expected:", e$expected))
        }
      }
    }, error = function(err) {
      list(ok = FALSE, host = host_label, version = "", msg = conditionMessage(err))
    })

    attempts <- c(attempts, paste0("[", res$host, "] ", res$msg, if (nzchar(res$version)) paste0(" (", res$version, ")") else ""))
    if (isTRUE(res$ok)) {
      return(list(
        species = e$species,
        dataset = e$dataset,
        status = "OK",
        host = res$host,
        version = res$version,
        expected = e$expected,
        details = paste(attempts, collapse = " || ")
      ))
    }
  }

  list(
    species = e$species,
    dataset = e$dataset,
    status = "FAIL",
    host = "",
    version = "",
    expected = e$expected,
    details = paste(attempts, collapse = " || ")
  )
}

results <- lapply(entries, check_entry)

cat("species\tdataset\tstatus\thost\tversion\texpected\tdetails\n")
for (r in results) {
  cat(paste(
    r$species,
    r$dataset,
    r$status,
    r$host,
    r$version,
    r$expected,
    gsub("\t", " ", r$details),
    sep = "\t"
  ), "\n", sep = "")
}

n_fail <- sum(vapply(results, function(x) identical(x$status, "FAIL"), logical(1)))
cat("\nSUMMARY: ", length(results) - n_fail, " OK / ", n_fail, " FAIL\n", sep = "")

if (n_fail > 0) {
  quit(status = 1)
}
