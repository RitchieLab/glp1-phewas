# libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
    library(arrow)
    library(cli)
    library(glue)
    library(data.table)
    # library(CCS)
    library(dplyr)
    library(purrr)
    library(janitor)
})

options(datatable.print.class = TRUE)
options(datatable.print.trunc.cols = TRUE)

project_directory <- "/home/jupyter/workspaces/duplicateofglp1rapharmacogenomics"
data_directory <- file.path(project_directory, "data")
raw_directory <- file.path(data_directory, "raw")
processed_directory <- file.path(data_directory, "processed")
version <- "20250129"

dir.create(processed_directory, showWarnings = FALSE)

# functions --------------------------------------------------------------------
# simplifies saving later
psave <- function(obj, file) {
    write_parquet(
        x = obj,
        sink = file,
        compression = "zstd",
        compression_level = 12
    )
}

remove_single_quotes <- function(dt) {
    stopifnot(is.data.table(dt)) # Ensure input is a data.table

    # Loop through character columns and remove single quotes
    char_cols <- names(dt)[sapply(dt, is.character)]
    dt[,
        (char_cols) := lapply(.SD, function(x) gsub("'", "", x)),
        .SDcols = char_cols
    ]

    return(dt[])
}

# data -------------------------------------------------------------------------
demo <- open_dataset(glue("{raw_directory}/demo/")) |> collect()
setDT(demo)
demo[, person_id := as.character(person_id)]
setkey(demo, person_id)

cohort <- fread(file.path(processed_directory, "cohort.csv"))

icd <- open_dataset(glue("{raw_directory}/icd/")) |>
    filter(
        source_vocabulary == "ICD10CM" & person_id %in% cohort[, person_id]
    ) |>
    select(
        person_id,
        vocabulary_id = source_vocabulary,
        code = source_concept_code,
        date = condition_start_datetime
    ) |>
    collect()
setDT(icd)
icd[, `:=`(person_id = as.character(person_id), date = as.Date(date))]
setkey(icd, code)
icd[, code := gsub("\\.", "", code)]


ccsr_map <- fread(file.path(raw_directory, "DXCCSR_v2025-1.csv")) |>
    clean_names() |>
    remove_single_quotes()

# map ICD codes to CCSR --------------------------------------------------------
ccsr_res <- merge(
    icd,
    ccsr_map[, .(
        code = icd_10_cm_code,
        ccsr = ccsr_category_1
    )],
    by = "code",
    all = FALSE
)


# save
psave(
    ccsr_res,
    glue("{processed_directory}/CCSR_{version}.parquet")
)


cli_alert_success("ccsr processing complete 🥳")
