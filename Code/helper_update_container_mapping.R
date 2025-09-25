# Update container to tissue mapping in PK-Sim expression databases
# Mapping of tissue samples to PK-Sim containers is a rather biased
# process based on individual desssion making, thus alignment in the
# comunity and enabling a quick remapping might be beneficial

# set working and export directory
setwd(here::here())
PATH <- getwd()

source(paste0(PATH, "/Code/helper_SQL_Commands.R"))

# Load tab_container_tissue.txt
assign(
  x = "tab_container_tissue",
  value = tibble::tibble(read.table(paste0(PATH, "/Code/tab_container_tissue.txt"),
    header = 1, sep = "\t", quote = '"'
  ))
)

# Get all .expressionDB files in PK-Sim DBs folder
db_files <- list.files(
  path = "PK-Sim DBs",
  pattern = "\\.expressionDB$",
  full.names = TRUE,
  recursive = TRUE
)

# Iterate through each DB and update tab_container_tissue
for (db_path in db_files) {
  db_PKsim_conn <- DBI::dbConnect(RSQLite::SQLite(),
    db_path,
    synchronous = "off",
    cache_size = -1000
  )

  # Write updated table to PK-Sim DB
  DBI::dbWriteTable(
    conn = db_PKsim_conn,
    name = "tab_container_tissue",
    value = tab_container_tissue,
    row.names = FALSE,
    overwrite = TRUE,
    append = FALSE
  )
  # Remove old tab_container_tissue table if exists
  # if ("tab_container_tissue" %in% DBI::dbListTables(db_PKsim_conn)) {
  #  DBI::dbRemoveTable(db_PKsim_conn, "tab_container_tissue")
  # }
  #
  # Add/update the view for tab_container_tissue
  # DBI::dbExecute(
  #  conn = db_PKsim_conn,
  #  statement = "DROP VIEW IF EXISTS QRY_CONTAINER_TISSUE"
  # )
  #
  # Add/update the view for tab_container_tissue
  # DBI::dbExecute(
  #  conn = db_PKsim_conn,
  #  statement = "DROP TABLE IF EXISTS TAB_CONTAINER_TISSUE;"
  # )
  #
  # DBI::dbExecute(
  #  conn = db_PKsim_conn,
  #  statement = "CREATE VIEW QRY_CONTAINER_TISSUE AS SELECT TAB_CONTAINER_TISSUE.CONTAINER as CONTAINER, TAB_CONTAINER_TISSUE.TISSUE as TISSUE FROM TAB_CONTAINER_TISSUE"
  #  # "CREATE VIEW QRY_CONTAINER_TISSUE AS SELECT container, tissue FROM tab_container_tissue"
  # )
  DBI::dbDisconnect(db_PKsim_conn)
}
