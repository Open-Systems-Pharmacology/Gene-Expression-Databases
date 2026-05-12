# Shared BioMart dataset/host/version pins for reproducibility-critical mapping.
# This is the single source of truth used by both PrepareBioMarts and
# validation utilities.
get_biomart_entries <- function() {
  list(
    Cat = list(species = "Cat", dataset = "fcatus_gene_ensembl", hosts = c("https://may2021.archive.ensembl.org", "https://www.ensembl.org"), expected = "^Felis_catus_9\\.0$"),
    Cattle = list(species = "Cattle", dataset = "btaurus_gene_ensembl", hosts = c("https://may2021.archive.ensembl.org", "https://www.ensembl.org"), expected = "^ARS-UCD1\\.2$"),
    Chicken = list(species = "Chicken", dataset = "ggallus_gene_ensembl", hosts = c("https://apr2022.archive.ensembl.org", "https://may2021.archive.ensembl.org"), expected = "^GRCg6a$"),
    Dog = list(species = "Dog", dataset = "clfamiliaris_gene_ensembl", hosts = c("https://may2021.archive.ensembl.org", "https://www.ensembl.org"), expected = "^CanFam3\\.1$"),
    Goat = list(species = "Goat", dataset = "chircus_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^ARS1$"),
    Guineapig = list(species = "Guineapig", dataset = "cporcellus_gene_ensembl", hosts = c("https://may2025.archive.ensembl.org", "https://www.ensembl.org"), expected = "^Cavpor3\\.0$"),
    Horse = list(species = "Horse", dataset = "ecaballus_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^EquCab3\\.0$"),
    Human = list(species = "Human", dataset = "hsapiens_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^GRCh38\\.p14$"),
    Minipig = list(species = "Minipig", dataset = "sscrofa_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^Sscrofa11\\.1$"),
    Monkey_fascicularis = list(species = "Monkey_fascicularis", dataset = "mfascicularis_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^Macaca_fascicularis_6\\.0$"),
    Monkey_mulatta = list(species = "Monkey_mulatta", dataset = "mmulatta_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^Mmul_10$"),
    Monkey_PigTailed = list(species = "Monkey_PigTailed", dataset = "mnemestrina_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^Mnem_1\\.0$"),
    Mouse = list(species = "Mouse", dataset = "mmusculus_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^GRCm39$"),
    Rabbit = list(species = "Rabbit", dataset = "ocuniculus_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^OryCun2\\.0$"),
    Rat = list(species = "Rat", dataset = "rnorvegicus_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^GRCr8$"),
    Sheep = list(species = "Sheep", dataset = "oaries_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^ARS-UI_Ramb_v3\\.0$"),
    Turkey = list(species = "Turkey", dataset = "mgallopavo_gene_ensembl", hosts = c(NA_character_, "https://www.ensembl.org", "https://useast.ensembl.org", "https://uswest.ensembl.org", "https://asia.ensembl.org"), expected = "^Turkey_5\\.1$"),
    Zebrafish = list(species = "Zebrafish", dataset = "drerio_gene_ensembl", hosts = c("https://oct2024.archive.ensembl.org", "https://www.ensembl.org"), expected = "^GRCz11$")
  )
}
