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

project_directory   <- "/home/jupyter/workspaces/duplicateofglp1rapharmacogenomics"
data_directory      <- file.path(project_directory, "data")
raw_directory       <- file.path(data_directory, "raw")
processed_directory <- file.path(data_directory, "processed")
version             <- "20250129"

dir.create(processed_directory, showWarnings = FALSE)

# functions --------------------------------------------------------------------
# simplifies saving later
psave <- function(obj, file) {
    write_parquet(x = obj, sink = file, compression = "zstd", compression_level = 12)
}

# remove phecode sex-subject sex discordant observations
restrictPhecodesBySex <- function(phecodes, id_sex, phe_sex) {
    # merge in subject sex
    data <- merge.data.table(phecodes, id_sex, by = "person_id", all.x = TRUE)
    # merge in phecode sex
    phe_sex[, phe_sex := fcase(
        male_only, "M",
        female_only, "F",
        default = "B"
    )]
    data <- merge.data.table(data, phe_sex[, .(phecode, phe_sex)], all.x = TRUE, by = "phecode")
    # drop discordant observations
    data[
        phe_sex == "B" |
        (sex == "M" & phe_sex == "M") |
        (sex == "F" & phe_sex == "F"),
    ][, !c("sex", "phe_sex")]
}

# turn phenome into phecode indicator matrix
generate_pim <- function(
    phe_dat,
    id_var     = "person_id",
    chunk_size = 10000,
    id_sex     = NULL,
    phe_sex    = NULL,
    progress   = TRUE,
    phe_var    = "phecode",
    n_var      = "N",
    sex_var   = "sex"
  ) {
    # generate ID chunks
    ids       <- unique(phe_dat[[id_var]])
    id_chunks <- base::split(ids, ceiling(seq_along(ids) / chunk_size))
        
    # iterate over the chunks
    res_list <- lapply(
        cli_progress_along(seq_along(id_chunks)), \(i) {
            data_sub <- phe_dat[phe_dat[[id_var]] %in% id_chunks[[i]], ]
            dcast(
                data_sub,
                formula   = get(id_var) ~ get(phe_var),
                value.var = n_var
            )
        }
    )
    cli_progress_done()
    
    res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
    res[is.na(res)] <- 0
    
    if (!is.null(id_sex) && !is.null(phe_sex)) {
        male_codes   <- phe_sex[male_only == TRUE, phecode]
        male_codes   <- male_codes[male_codes %in% names(res)]
        female_codes <- phe_sex[female_only == TRUE, phecode]
        female_codes <- female_codes[female_codes %in% names(res)]
        
        male_ids   <- id_sex[tolower(substr(id_sex[[sex_var]], 1, 1)) == "m", ][[id_var]]
        female_ids <- id_sex[tolower(substr(id_sex[[sex_var]], 1, 1)) == "f", ][[id_var]]
        
        for (var in male_codes) {
        set(res, i = which(res[[id_var]] %in% female_ids), j = var, value = NA)
        }
        for (var in female_codes) {
        set(res, i = which(res[[id_var]] %in% male_ids), j = var, value = NA)
        }
    }
    
    return(res)
    
}

# summarize phecode indicator matrix
quick_pim_summary <- function(pim, names = phecodeX_labels[, phenotype]) {
    result <- rbindlist(
        lapply(
            names(pim)[names(pim) %in% names],
            \(x) {
                data.table(
                    phecode = x,
                    cases    = sum(pim[[x]] == 1, na.rm = TRUE),
                    controls = sum(pim[[x]] == 0, na.rm = TRUE),
                    missing  = sum(is.na(pim[[x]]), na.rm = TRUE)
                )
            }
        ), use.names = TRUE, fill = TRUE
    )
    result[, total := cases + controls + missing][]
}

# data -------------------------------------------------------------------------
demo <- open_dataset(glue("{raw_directory}/demo/")) |> collect()
setDT(demo)
demo[, person_id := as.character(person_id)]
setkey(demo, person_id)

icd <- open_dataset(glue("{raw_directory}/icd/")) |> 
    select(
        person_id,
        vocabulary_id = source_vocabulary,
        code          = source_concept_code,
        date          = condition_start_datetime
    ) |>
    collect()
setDT(icd)
icd[, `:=` (person_id = as.character(person_id), date = as.Date(date))]
setkey(icd, person_id)

# load phecode mapping tables
phecodeX_labels     <- fread("https://github.com/PheWAS/PhecodeX/raw/main/phecodeX_R_labels.csv", showProgress = FALSE)
phecodeX_rollup_map <- fread("https://github.com/PheWAS/PhecodeX/raw/main/phecodeX_R_rollup_map.csv", showProgress = FALSE)
phecodeX_map        <- fread("https://github.com/PheWAS/PhecodeX/raw/main/phecodeX_R_map.csv", showProgress = FALSE) ## if you are using ICD-10 (not CM), load phecodeX_R_map_ICD_10_WHO.csv instead
phecodeX_sex        <- fread("https://github.com/PheWAS/PhecodeX/raw/main/phecodeX_R_sex.csv", showProgress = FALSE)

# map ICD codes to phecodes ----------------------------------------------------
phecodes <- mapCodesToPhecodes(
    input          = icd,
    vocabulary.map = phecodeX_map,
    rollup.map     = phecodeX_rollup_map,
    make.distinct  = TRUE
)

# remove phecode sex-subject sex discordant observations
id_sex <- demo[, .(
    person_id,
    sex = fcase(
        gender == "Female", "F",
        gender == "Male", "M"
        )
    )]

phe_res <- restrictPhecodesBySex(phecodes, id_sex, phecodeX_sex)

# identify first occurrence of each phecode
setkey(phe_res, person_id, phecode, date)
phe_res_min <- phe_res[phe_res[, .I[1], .(person_id, phecode)]$V1]

# count occurrences of each phecode
phe_res_n   <- phe_res[, .N, .(person_id, phecode)]

# generate dsb-version of phenome
phe_dsb <- merge(
    phe_res_min,
    demo[, .(person_id, dob = as.Date(date_of_birth))],
    by    = "person_id",
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

# #### REMOVE R OBJECTS BEFORE PROCEEDING ####
# # otherwise you will run out of memory
# rm(phe_res, phe_res_min, phe_dsb)
# gc()

# # generate phecode indicator matrix
# # single occurrence counts
# pim1 <- generate_pim(
#     phe_dat    = phe_res_n[, .(person_id, phecode, N = 1)],
#     phe_sex    = phecodeX_sex,
#     id_sex     = id_sex,
#     chunk_size = 10000
# )
# pim1_sum <- quick_pim_summary(pim1)
# psave(
#     pim1,
#     glue("{processed_directory}/pim1_{version}.parquet")
# )
# rm(pim1)
# gc()

# # two occurrences count
# pim2 <- generate_pim(
#     phe_dat    = phe_res_n[N >= 2][, .(person_id, phecode, N = 1)],
#     phe_sex    = phecodeX_sex,
#     id_sex     = id_sex,
#     chunk_size = 10000
# )
# pim2_sum <- quick_pim_summary(pim2)
# psave(
#     pim2,
#     glue("{processed_directory}/pim2_{version}.parquet")
# )
# rm(pim2)
# gc()

# phe_sum <- merge(
#     pim1_sum[, .(phecode, n1_cases = cases, n1_controls = controls, n1_miss = missing)],
#     pim2_sum[, .(phecode, n2_cases = cases, n2_controls = controls, n2_miss = missing)]
# )[, `:=` (
#     n1_prop = n1_cases / (n1_cases + n1_controls),
#     n2_prop = n2_cases / (n2_cases + n2_controls)
# )][]

# table_directory <- file.path(project_directory, "tables")
# dir.create(table_directory, showWarnings = FALSE)

# fwrite(phe_sum, file.path(table_directory, "aou_phe_sum.txt"))

cli_alert_success("phenome processing complete 🥳")
