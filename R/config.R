# Project configuration
config <- list(
  project_directory = here::here(), # Use here() for portability
  version = "20251103_r2r",

  # Paths
  paths = list(
    data = "data",
    raw = "data/raw",
    processed = "data/processed",
    results = "results",
    tables = "tables"
  ),

  # Data.table options
  dt_options = list(
    datatable.print.class = TRUE,
    datatable.print.trunc.cols = TRUE,
    datatable.print.nrows = 5
  )
)

# Apply configuration
do.call(options, config$dt_options)
