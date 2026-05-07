# Add indizes to all PK-Sim expression databases

setwd(here::here())
PATH <- getwd()

source(paste0(PATH, "/Code/helper_SQL_Commands.R")) # loads INDIZES

db_files <- list.files(
  path = "PK-Sim DBs",
  pattern = "\\.expressionDB$",
  full.names = TRUE,
  recursive = TRUE
)

for (db_path in db_files) {
  db_PKsim_conn <- DBI::dbConnect(RSQLite::SQLite(),
    db_path,
    synchronous = "off",
    cache_size = -1000
  )

  print(paste0("Remove all indexes from '", db_path, "' PK-Sim expression database"))
  # Get all indexes in the database
  idx_info <- DBI::dbGetQuery(db_PKsim_conn, "SELECT name FROM sqlite_master WHERE type = 'index';")
  for (idx in idx_info$name) {
    tryCatch(
      DBI::dbExecute(db_PKsim_conn, paste0("DROP INDEX IF EXISTS ", idx)),
      error = function(e) warning(paste("Index removal failed for", db_path, ":", idx, "\n", e$message))
    )
  }

  print(paste0("Set table indizes for '", db_path, "' PK-Sim expression database"))
  for (i in seq_along(INDIZES)) {
    tryCatch(
      DBI::dbExecute(db_PKsim_conn, INDIZES[i]),
      error = function(e) warning(paste("Index creation failed for", db_path, ":", INDIZES[i], "\n", e$message))
    )
  }

  DBI::dbDisconnect(db_PKsim_conn)
}
