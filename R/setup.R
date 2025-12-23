# Common library loader with conditional loading
load_libraries <- function(packages) {
  lapply(
    packages,
    \(pkg)
      suppressPackageStartupMessages({
        library(pkg, character.only = TRUE)
      })
  ) |>
    invisible()
}

# Define package groups
pkg_groups <- list(
  core = c("arrow", "cli", "glue", "data.table", "dplyr", "here"),
  phewas = "PheWAS",
  bigquery = "bigrquery",
  analysis = c(
    "comorbidity",
    "MatchIt",
    "survival",
    "survRM2",
    "caret",
    "survival",
    "survminer",
    "coxphf"
  ),
  summarize = c("gtsummary", "gt", "labelled"),
  plotting = c("ggplot2", "patchwork"),
  utils = c("purrr", "janitor", "optparse", "furrr", "future")
)
