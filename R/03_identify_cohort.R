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
version <- "20250129"

dir.create(processed_directory, showWarnings = FALSE)

# demographics -----------------------------------------------------------------
demo <- open_dataset(file.path(raw_directory, "demo")) |>
    collect()
setDT(demo)
demo[, `:=`(
    person_id = as.character(person_id),
    dob = as.Date(date_of_birth),
    sex = fcase(
        sex_at_birth == "Male",
        "Male",
        sex_at_birth == "Female",
        "Female"
    )
)]
setkey(demo, person_id)
print(
    paste0(
        "There are ",
        demo[, uniqueN(person_id)],
        " individuals in the demographics dataset."
    )
)

demo <- demo[!is.na(sex), ]
print(
    paste0(
        "There are ",
        demo[, uniqueN(person_id)],
        " individuals in the demographics dataset with non-missing Male/Female sex assigned at birth."
    )
)

# identify subset --------------------------------------------------------------
phe_res <- open_dataset(glue(
    "{processed_directory}/phenome_{version}.parquet"
)) |>
    collect()
setDT(phe_res)
print(
    paste0(
        "There are ",
        phe_res[, uniqueN(person_id)],
        " individuals in the phenome dataset."
    )
)

# individuals with t2d
phe_res_t2d <- phe_res[
    phecode == "EM_202.2",
    .(t2d_date = min(date)),
    person_id
]
print(
    paste0(
        "There are ",
        phe_res_t2d[, uniqueN(person_id)],
        " individuals with T2D."
    )
)

drug <- open_dataset(glue("{processed_directory}/drug/")) |>
    collect()
setDT(drug)
drug[, date := as.Date(drug_exposure_start_datetime)]

drug_sub <- drug[, .(drug_date = max(date)), person_id]
print(
    paste0(
        "There are ",
        drug_sub[, uniqueN(person_id)],
        " individuals with drug data."
    )
)

# individuals with demographic, t2d, and drug data
merged <- merge(
    demo,
    phe_res_t2d,
    by = "person_id",
    all = FALSE
) |>
    merge(
        unique(drug_sub),
        by = "person_id",
        all = FALSE
    )
print(
    paste0(
        "There are ",
        merged[, uniqueN(person_id)],
        " individuals with demographic, T2D, and drug data."
    )
)

# individuals with qualifying drug after t2d diagnosis
merged <- merged[drug_date > t2d_date, ]
print(
    paste0(
        "There are ",
        merged[, uniqueN(person_id)],
        " individuals with qualifying drug after T2D diagnosis."
    )
)

# individuals with qualifying drug on or after January 1, 2018
merged <- merged[drug_date >= as.Date("2018-01-01"), ]
print(
    paste0(
        "There are ",
        merged[, uniqueN(person_id)],
        " individuals with qualifying drug on or after January 1, 2018."
    )
)

# one-year follow-up subset ----------------------------------------------------
icd <- open_dataset(file.path(raw_directory, "icd")) |>
    mutate(person_id = as.character(person_id)) |>
    filter(person_id %in% merged[, unique(person_id)]) |>
    select(person_id, date = condition_start_datetime) |>
    collect()
setDT(icd)
icd[, `:=`(
    person_id = as.character(person_id),
    date = as.Date(date)
)]
setkey(icd, person_id)
icd <- merge(icd, merged[, .(person_id, drug_date)], by = "person_id")
icd <- icd[date <= drug_date, ]
icd_n <- unique(icd[, .(person_id, date)])[, .(encounters = .N), person_id]
icd_min <- icd[icd[, .I[which.min(date)], person_id][, V1]][, .(
    person_id,
    icd_date = date,
    drug_date
)]
icd_min <- merge(
    icd_min,
    icd_n,
    by = "person_id"
)
icd_min[, followup := round(as.numeric(drug_date - icd_date) / 365.25, 1)]
icd_min <- icd_min[followup >= 1, ]
merged <- merge(
    merged,
    icd_min[, !c("drug_date")],
    by = "person_id",
    all = FALSE
)
print(
    paste0(
        "There are ",
        merged[, uniqueN(person_id)],
        " individuals with one-year follow-up."
    )
)

# drug subset ------------------------------------------------------------------
drug_2018 <- drug[
    date >= as.Date("2018-01-01") & person_id %in% merged[, unique(person_id)],
]
drug_2017 <- drug[
    date >= as.Date("2017-01-01") & person_id %in% merged[, unique(person_id)],
]

drug_min_2018 <- drug_2018[, .(min_date_2018 = min(date)), person_id]
drug_min_2017 <- drug_2017[, .(min_date_2017 = min(date)), person_id]

drug_min <- merge(drug_min_2018, drug_min_2017, by = "person_id", all = TRUE)

# remove individuals with qualifying drug within 120 days of 2018-based T0 prescription
rm_bc_hx_rx <- drug_min[
    min_date_2018 != min_date_2017 &
        (as.numeric(min_date_2018 - min_date_2017) < 120),
    person_id
]

merged <- merged[!person_id %in% rm_bc_hx_rx, ]
print(
    paste0(
        "There are ",
        merged[, uniqueN(person_id)],
        " individuals without disqualifying prescription within 120 days."
    )
)

mappings <- fread(
    file.path(project_directory, "tables", "omop_mappings.csv"),
    na.strings = c("", "NA")
)[!is.na(drug_class), ]

# remove individuals with without qualifying prescription on or after January 1, 2018
sub <- merge(
    drug[
        person_id %in% merged[, person_id] & date >= as.Date("2018-01-01"),
        .(person_id, drug_concept_id, date)
    ],
    mappings[, .(drug_concept_id = concept_id, drug_class, search_term)],
    by = "drug_concept_id",
    all.x = TRUE
)
sub <- sub[sub[, .I[date == min(date)], .(person_id, drug_class)][, V1], ]
sub <- sub[sub[, .I[date == min(date)], .(person_id)][, V1], ]
sub <- unique(sub[, .(person_id, date, drug_class, search_term)])
sub[, uniqueN(person_id)]

sub2 <- sub[!person_id %in% sub[, .N, person_id][N > 1, person_id], ]
sub2[, uniqueN(person_id)]
print(
    paste0(
        "There are ",
        sub2[, uniqueN(person_id)],
        " individuals with a unique qualifying prescription on or after January 1, 2018."
    )
)


# remove individuals with contraindications
phe_sub <- merge(
    phe_res[
        person_id %in% sub2[, person_id],
        .(person_id, phe_date = date, phecode)
    ],
    sub2,
    by = "person_id"
)[data.table::between(phe_date, date - (5 * 365.25), date), ]

drop_phe <- phe_sub[
    phecode %in%
        c(
            # MEN type II
            "GE_961.22",
            "GE_961.21",
            # gastroparesis
            "GI_516.4",
            # dialysis or kidney transplant
            "GU_582.3",
            "SS_847.2"
        ),
][, .(phe_date = max(phe_date)), person_id]

icd <- open_dataset(file.path(raw_directory, "icd")) |>
    mutate(person_id = as.character(person_id)) |>
    filter(
        person_id %in%
            sub2[, person_id] &
            source_concept_code %in%
                c(
                    # medullary thyroid cancer
                    "C73",
                    "Z85.850",
                    # Hypoglycemia with coma
                    "E15",
                    "E13.641",
                    "E11.641"
                )
    ) |>
    select(
        person_id,
        code = source_concept_code,
        vocabulary_id = source_vocabulary,
        date = condition_start_datetime,
        description = source_concept_name
    ) |>
    collect()
setDT(icd)
icd[, date := as.Date(date)]
drop_icd <- merge(
    icd[, .(person_id, icd_date = date, code)],
    sub2,
    by = "person_id"
)[data.table::between(icd_date, date - (5 * 365.25), date), ][,
    .(icd_date = max(icd_date)),
    person_id
]
rm(icd)

lab <- open_dataset(file.path(raw_directory, "lab")) |>
    mutate(person_id = as.character(person_id)) |>
    filter(
        person_id %in% sub2[, person_id] & standard_concept_code %in% "77147-7"
    ) |>
    collect()
setDT(lab)

lab[, `:=`(
    value = as.numeric(value_as_number),
    date = as.Date(measurement_datetime)
)]

drop_lab <- merge(
    lab[!is.na(value), .(person_id, value, lab_date = date)],
    sub2,
    by = "person_id"
)[data.table::between(lab_date, date - (5 * 365.25), date), ][
    value < 30,
    .(person_id, lab_date)
]

contra_tab <- merge(
    drop_phe,
    drop_icd,
    by = "person_id",
    all = TRUE
) |>
    merge(
        drop_lab,
        by = "person_id",
        all = TRUE
    )

drop_lab <- lab[value < 30, .SD[which.min(date)], person_id][, .(
    id = person_id,
    date,
    value
)]

drop_contra <- rbindlist(
    list(drop_phe, drop_icd, drop_lab),
    fill = TRUE,
    use.names = TRUE
)[, .SD[which.min(date)], id][, .(id, contra_date = date)]

# individuals without contraindications
sub3 <- sub2[!person_id %in% drop_contra[, id], ]
print(
    paste0(
        "There are ",
        sub3[, uniqueN(person_id)],
        " individuals without contraindications."
    )
)
setnames(sub3, "date", "drug_date")

# add some additional information
sub3 <- merge(
    sub3,
    merged[, .(
        person_id,
        dob,
        t2d_date,
        t2d_age = round(
            as.numeric(t2d_date - as.Date(date_of_birth)) / 365.25,
            1
        ),
        encounters,
        followup
    )],
    by = "person_id"
)

death <- open_dataset(file.path(raw_directory, "death")) |>
    mutate(person_id = as.character(person_id)) |>
    filter(person_id %in% sub3[, person_id]) |>
    select(person_id, death_date) |>
    collect()
setDT(death)

drug_tmp <- merge(
    drug[
        person_id %in% sub3[, person_id] & date >= as.Date("2018-01-01"),
        .(person_id, date, concept_id = drug_concept_id)
    ],
    mappings[, .(concept_id, drug_class, search_term)],
    by = "concept_id"
)[,
    drug_class := fcase(
        drug_class == "GLP-1 RA",
        "glp1",
        drug_class == "DPP4i",
        "dpp4",
        drug_class == "SGLT2i",
        "sglt2"
    )
]

sema_min_date <- dcast(
    drug_tmp[
        search_term == "semaglutide",
        .(min_date = min(date)),
        .(person_id, search_term)
    ],
    person_id ~ paste0(search_term, "_min_date"),
    value.var = "min_date"
)
sema_max_date <- dcast(
    drug_tmp[
        search_term == "semaglutide",
        .(max_date = max(date)),
        .(person_id, search_term)
    ],
    person_id ~ paste0(search_term, "_max_date"),
    value.var = "max_date"
)

other_min_dates <- dcast(
    drug_tmp[, .(min_date = min(date)), .(person_id, drug_class)],
    person_id ~ paste0(drug_class, "_min_date"),
    value.var = "min_date"
)
other_max_dates <- dcast(
    drug_tmp[, .(max_date = max(date)), .(person_id, drug_class)],
    person_id ~ paste0(drug_class, "_max_date"),
    value.var = "max_date"
)

date_dt <- Reduce(
    \(x, y) merge(x, y, by = "person_id", all = TRUE),
    list(
        sema_min_date,
        sema_max_date,
        other_min_dates,
        other_max_dates,
        death,
        sub3[, .(person_id, drug_class, search_term)]
    )
)[, cutoff_date := as.Date("2023-10-01")]


# set semaglutide_date to missing for individuals with alternate drug within 120 days
date_dt[
    semaglutide_min_date != glp1_min_date &
        (semaglutide_min_date -
            pmax(glp1_min_date, dpp4_min_date, sglt2_min_date, na.rm = TRUE)) <
            120,
    semaglutide_min_date := NA
]

# max dates for per-protocol analyses
# GLP1 RA vs SGLT2i
date_dt[
    drug_class == "GLP-1 RA",
    max_glp1_sglt2_date := pmin(
        glp1_max_date + 90,
        sglt2_min_date,
        death_date,
        cutoff_date,
        na.rm = TRUE
    )
]
date_dt[
    drug_class == "SGLT2i",
    max_glp1_sglt2_date := pmin(
        glp1_min_date,
        sglt2_max_date + 90,
        death_date,
        cutoff_date,
        na.rm = TRUE
    )
]

# GLP1 RA vs DPP4i
date_dt[
    drug_class == "GLP-1 RA",
    max_glp1_dpp4_date := pmin(
        glp1_max_date + 90,
        dpp4_min_date,
        death_date,
        cutoff_date,
        na.rm = TRUE
    )
]
date_dt[
    drug_class == "DPP4i",
    max_glp1_dpp4_date := pmin(
        glp1_min_date,
        dpp4_max_date + 90,
        death_date,
        cutoff_date,
        na.rm = TRUE
    )
]

# Semaglutide vs SGLT2i
date_dt[
    search_term == "semaglutide",
    max_sema_sglt2_date := pmin(
        semaglutide_max_date + 90,
        sglt2_min_date,
        death_date,
        cutoff_date,
        na.rm = TRUE
    )
]
date_dt[
    drug_class == "SGLT2i",
    max_sema_sglt2_date := pmin(
        semaglutide_min_date,
        sglt2_max_date + 90,
        death_date,
        cutoff_date,
        na.rm = TRUE
    )
]

# Semaglutide vs DPP4i
date_dt[
    search_term == "semaglutide",
    max_sema_dpp4_date := pmin(
        semaglutide_max_date + 90,
        dpp4_min_date,
        death_date,
        cutoff_date,
        na.rm = TRUE
    )
]
date_dt[
    drug_class == "DPP4i",
    max_sema_dpp4_date := pmin(
        semaglutide_min_date,
        dpp4_max_date + 90,
        death_date,
        cutoff_date,
        na.rm = TRUE
    )
]

sub3 <- merge(
    sub3,
    date_dt[, !c("drug_class", "search_term", "cutoff_date")],
    by = "person_id"
)

fwrite(
    sub3,
    file.path(processed_directory, "cohort.csv"),
    row.names = FALSE
)

cli_alert_success("Cohort identified and saved. 🥳")
