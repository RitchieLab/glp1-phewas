# libraries --------------------------------------------------------------------
source("R/config.R")
source("R/setup.R")
source("R/fn/helper-functions.R")
load_libraries(c(pkg_groups$core, pkg_groups$phewas, pkg_groups$utils))

project_directory <- config$project_directory
data_directory <- file.path(project_directory, config$paths$data)
raw_directory <- file.path(project_directory, config$paths$raw)
processed_directory <- file.path(project_directory, config$paths$processed)
tables_directory <- file.path(project_directory, config$paths$tables)
version <- config$version

dir.create(processed_directory, showWarnings = FALSE)

# data -------------------------------------------------------------------------
demo <- open_dataset(glue("{raw_directory}/demo_{version}/")) |> collect()
setDT(demo)
demo[, person_id := as.character(person_id)]
setkey(demo, person_id)

icd <- open_dataset(glue("{raw_directory}/icd_{version}/")) |>
  select(
    person_id,
    vocabulary_id = source_vocabulary,
    code = source_concept_code,
    date = condition_start_datetime
  ) |>
  collect()
setDT(icd)
icd[, `:=`(person_id = as.character(person_id), date = as.Date(date))]
setkey(icd, person_id)

# load phecode mapping tables
## local files downloaded from https://github.com/PheWAS/PhecodeX on December 9, 2025
phecodeX_labels <- fread(glue("{tables_directory}/phecodeX_labels.csv"))
phecodeX_rollup_map <- fread(glue("{tables_directory}/phecodeX_rollup_map.csv"))
phecodeX_map <- fread(glue("{tables_directory}/phecodeX_map.csv"))
phecodeX_sex <- fread(glue("{tables_directory}/phecodeX_sex.csv"))

# map ICD codes to phecodes ----------------------------------------------------
phecodes <- mapCodesToPhecodes(
  input = icd,
  vocabulary.map = phecodeX_map,
  rollup.map = phecodeX_rollup_map,
  make.distinct = TRUE
)

# remove phecode sex-subject sex discordant observations
id_sex <- demo[, .(
  person_id,
  sex = fcase(
    gender == "Female",
    "F",
    gender == "Male",
    "M"
  )
)]

phe_res <- restrictPhecodesBySex(phecodes, id_sex, phecodeX_sex)

# identify first occurrence of each phecode
setkey(phe_res, person_id, phecode, date)
phe_res_min <- phe_res[phe_res[, .I[1], .(person_id, phecode)]$V1]

# generate dsb-version of phenome
phe_dsb <- merge(
  phe_res_min,
  demo[, .(person_id, dob = as.Date(date_of_birth))],
  by = "person_id",
  all.x = TRUE
)[, dsb := as.numeric(date - dob)][, .(person_id, phecode, dsb)]

# set keys
setkey(phe_res, person_id)
setkey(phe_res_min, person_id)
setkey(phe_dsb, person_id)

# save
psave(
  phe_res,
  glue("{processed_directory}/phenome_{version}.parquet")
)
psave(
  phe_res_min,
  glue("{processed_directory}/phenome_first_{version}.parquet")
)
psave(
  phe_dsb,
  glue("{processed_directory}/phenome_first_dsb_{version}.parquet")
)

cli_alert_success("phenome processing complete 🥳")
