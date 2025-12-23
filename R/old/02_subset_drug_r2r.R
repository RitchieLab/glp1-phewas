# libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
  library(arrow)
  library(cli)
  library(glue)
  library(data.table)
  library(PheWAS)
  library(dplyr)
  library(purrr)
})

options(datatable.print.class = TRUE)
options(datatable.print.trunc.cols = TRUE)

project_directory <- "/home/jupyter/workspaces/duplicateofglp1rapharmacogenomics"
data_directory <- file.path(project_directory, "data")
raw_directory <- file.path(data_directory, "raw")
processed_directory <- file.path(data_directory, "processed")
version <- "20251103_r2r"

dir.create(processed_directory, showWarnings = FALSE)

drug <- open_dataset(glue("{raw_directory}/drug/")) |>
  collect()
setDT(drug)
drug[, `:=`(
  person_id = as.character(person_id),
  date = as.Date(drug_exposure_start_datetime)
)]
setkey(drug, person_id)

double_check_codes <- fread("tables/t2d_glp1_drug.codes")

write_dataset(
  drug[drug_concept_id %in% double_check_codes[, concept_id], ],
  glue("{processed_directory}/drug_{version}"),
  format = "parquet",
  compression = "zstd",
  compression_level = 12
)
