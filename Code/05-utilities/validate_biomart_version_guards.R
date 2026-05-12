#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(biomaRt))

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_file_arg) == 1) {
  script_path <- normalizePath(sub("^--file=", "", script_file_arg), winslash = "/")
} else {
  script_path <- normalizePath("Code/05-utilities/validate_biomart_version_guards.R", winslash = "/")
}

repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/")
source(file.path(repo_root, "Code/03-helpers/helper_BioMart_Pins.R"))
source(file.path(repo_root, "Code/03-helpers/helper_Species.R"))

entries <- get_biomart_entries()
entry_species <- sort(vapply(entries, function(e) e$species, character(1)))
if (!identical(unname(entry_species), sort(ALL_SPECIE))) {
  stop("Shared BioMart pin definitions are out of sync with helper_Species.R")
}

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
