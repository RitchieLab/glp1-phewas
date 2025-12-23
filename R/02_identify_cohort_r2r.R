# libraries --------------------------------------------------------------------
source("R/config.R")
source("R/setup.R")
source("R/fn/helper-functions.R")
load_libraries(c(pkg_groups$core, pkg_groups$phewas, pkg_groups$utils))

project_directory <- config$project_directory
data_directory <- file.path(project_directory, config$paths$data)
raw_directory <- file.path(project_directory, config$paths$raw)
processed_directory <- file.path(project_directory, config$paths$processed)
version <- config$version

dir.create(processed_directory, showWarnings = FALSE)

# demographics -----------------------------------------------------------------
demo <- open_dataset(file.path(raw_directory, glue("demo_{version}"))) |>
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

double_check_codes <- fread(
  file.path(
    project_directory,
    "data",
    "raw",
    glue("drug_concept_ids_{config$version}.tsv")
  )
)

drug <- open_dataset(glue("{raw_directory}/drug_{version}/")) |>
  filter(drug_concept_id %in% double_check_codes[, descendant_id]) |>
  collect()
setDT(drug)
drug[, `:=`(
  person_id = as.character(person_id),
  date = as.Date(drug_exposure_start_datetime)
)]

drug_sub <- drug[
  date >= as.Date("2018-01-01"),
  .(drug_date_min = min(date), drug_date_max = max(date)),
  person_id
]
print(
  paste0(
    "There are ",
    drug[, uniqueN(person_id)],
    " individuals with drug data."
  )
)

print(
  paste0(
    "There are ",
    drug_sub[, uniqueN(person_id)],
    " individuals with a qualifying drug on or after January 1, 2018."
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
merged <- merged[drug_date_max > t2d_date, ]
print(
  paste0(
    "There are ",
    merged[, uniqueN(person_id)],
    " individuals with qualifying drug after T2D diagnosis."
  )
)

# individuals with qualifying drug on or after January 1, 2018
merged <- merged[drug_date_min >= as.Date("2018-01-01"), ]
print(
  paste0(
    "There are ",
    merged[, uniqueN(person_id)],
    " individuals with qualifying drug on or after January 1, 2018."
  )
)

# one-year follow-up subset ----------------------------------------------------
icd <- open_dataset(file.path(raw_directory, glue("icd_{version}"))) |>
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
icd <- merge(icd, merged[, .(person_id, drug_date_min)], by = "person_id")
icd <- icd[date <= drug_date_min, ]

icd_n <- unique(icd[, .(person_id, date)])[,
  .(n_unique_dx_days = .N),
  person_id
]
icd_min <- icd[icd[, .I[which.min(date)], person_id][, V1]][, .(
  person_id,
  icd_date = date,
  drug_date_min
)]
icd_min <- merge(
  icd_min,
  icd_n,
  by = "person_id"
)
icd_min[, followup := round(as.numeric(drug_date_min - icd_date) / 365.25, 1)]
icd_min <- icd_min[followup >= 1, ]
merged <- merge(
  merged,
  icd_min[, !c("drug_date_min")],
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

# add visits -------------------------------------------------------------------
visits <- open_dataset(file.path(raw_directory, glue("visit_{version}"))) |>
  rename(person_id = PERSON_ID, date = visit_start_datetime) |>
  mutate(person_id = as.character(person_id)) |>
  filter(person_id %in% merged[, unique(person_id)]) |>
  unique() |>
  collect()
setDT(visits)
visits[, `:=`(
  person_id = as.character(person_id),
  date = as.Date(date)
)]
setkey(visits, person_id)

visits <- merge(visits, merged[, .(person_id, drug_date_min)], by = "person_id")
visits <- visits[date <= drug_date_min, ]

visits <- unique(visits[, .(person_id, date)])[, .(n_visits = .N), person_id]

merged <- merge.data.table(merged, visits, by = "person_id", all.x = TRUE)
merged[is.na(n_visits), n_visits := 0]

# drug subset ------------------------------------------------------------------
drug_2018 <- drug[
  date >= as.Date("2018-01-01") & person_id %in% merged[, unique(person_id)],
]
drug_2017 <- drug[
  date < as.Date("2018-01-01") &
    date >= as.Date("2017-01-01") &
    person_id %in% merged[, unique(person_id)],
]

drug_min_2018 <- drug_2018[, .(min_date_2018 = min(date)), person_id]
drug_max_2017 <- drug_2017[, .(max_date_2017 = max(date)), person_id]

drug_min <- merge(drug_min_2018, drug_max_2017, by = "person_id", all = TRUE)

# remove individuals with qualifying drug within 120 days of 2018-based T0 prescription
rm_bc_hx_rx <- drug_min[
  !is.na(max_date_2017) &
    (as.numeric(min_date_2018 - max_date_2017) < 120),
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

all_descendants <- fread(
  file.path(
    project_directory,
    "data",
    "raw",
    glue("drug_concept_ids_{version}.tsv")
  )
)

# remove individuals with without qualifying prescription on or after January 1, 2018
sub <- merge(
  drug[
    person_id %in% merged[, person_id] & date >= as.Date("2018-01-01"),
    .(person_id, drug_concept_id, date)
  ],
  all_descendants[, .(
    drug_concept_id = descendant_id,
    drug_ingredient_name = tolower(name),
    drug_class
  )],
  by = "drug_concept_id",
  allow.cartesian = TRUE,
  all.x = TRUE
)

sub <- sub[sub[, .I[date == min(date)], .(person_id)][, V1], ]
sub <- unique(sub[, .(person_id, date, drug_class, drug_ingredient_name)])
sub[, uniqueN(person_id)]

# if multiple ids for same drug ingredient - keep and assign drug ingredient
# if ids for different drugs of same class - remove
# if ids for different drug classes - remove
sub[, n_ingredients := uniqueN(drug_ingredient_name), by = person_id]
sub[, n_classes := uniqueN(drug_class), by = person_id]

sub <- unique(sub[
  n_ingredients == 1 & # all same drug_name (ingredient) - KEEP
    n_classes == 1 # all same class - KEEP
])

print(
  paste0(
    "There are ",
    sub[, uniqueN(person_id)],
    " individuals with a unique qualifying drug prescription on or after January 1, 2018"
  )
)

sub2 <- unique(sub[, .(person_id, date, drug_ingredient_name, drug_class)])


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

icd <- open_dataset(file.path(raw_directory, glue("icd_{version}"))) |>
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

lab <- open_dataset(file.path(raw_directory, glue("lab_{version}"))) |>
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
    n_visits,
    n_unique_dx_days,
    followup
  )],
  by = "person_id"
)

death <- open_dataset(file.path(raw_directory, glue("death_{version}"))) |>
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
  all_descendants[, .(
    concept_id = descendant_id,
    drug_class,
    drug_ingredient_name = tolower(name)
  )],
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
    drug_ingredient_name == "semaglutide",
    .(min_date = min(date)),
    .(person_id, drug_ingredient_name)
  ],
  person_id ~ paste0(drug_ingredient_name, "_min_date"),
  value.var = "min_date"
)
sema_max_date <- dcast(
  drug_tmp[
    drug_ingredient_name == "semaglutide",
    .(max_date = max(date)),
    .(person_id, drug_ingredient_name)
  ],
  person_id ~ paste0(drug_ingredient_name, "_max_date"),
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
    sub3[, .(person_id, drug_class, drug_ingredient_name)]
  )
)[, cutoff_date := as.Date("2023-10-01")]


# set semaglutide_date to missing for individuals in GLP-1 RA group whose first drug wasn't semaglutide
date_dt[
  drug_ingredient_name != "semaglutide",
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
  drug_ingredient_name == "semaglutide",
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
  drug_ingredient_name == "semaglutide",
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
  date_dt[, !c("drug_class", "drug_ingredient_name", "cutoff_date")],
  by = "person_id"
)

# metformin and insulin data
mni <- open_dataset(file.path(
  raw_directory,
  glue("metformin_and_insulins_{version}")
)) |>
  filter(person_id %in% sub3[, unique(person_id)]) |>
  select(
    person_id,
    standard_concept_name,
    drug_concept_id,
    date = drug_exposure_start_datetime
  ) |>
  collect()
setDT(mni)
mni[, `:=`(
  person_id = as.character(person_id),
  date = as.Date(date)
)]

mni_ids <- fread(file.path(
  project_directory,
  "data",
  "raw",
  glue("metformin_and_insulins_ids_{version}.tsv")
))

mni <- merge(
  mni,
  mni_ids[, .(drug_concept_id = descendant_id, drug_class)],
  allow.cartesian = TRUE,
  by = "drug_concept_id"
)

mni <- merge(
  mni,
  sub3[, .(person_id, index_date = drug_date)],
  by = "person_id"
)

ins <- mni[
  date <= index_date &
    date >= (index_date - 365 * 5) &
    drug_class == "Insulins",
  {
    min_date = min(date)
    max_date = max(date)
    duration = as.integer(max_date - min_date)
    list(
      ins_min_date = min_date,
      ins_max_date = max_date,
      ins_dur = duration,
      insulin = 1L
    )
  },
  person_id
]
print(glue("There are {ins[, uniqueN(person_id)]} people with insulin data"))

mets <- mni[
  date <= index_date &
    date >= (index_date - 365 * 5) &
    drug_class == "Metformin",
  {
    min_date = min(date)
    max_date = max(date)
    duration = as.integer(max_date - min_date)
    list(
      met_min_date = min_date,
      met_max_date = max_date,
      met_duration = duration,
      metformin = 1L
    )
  },
  person_id
]
print(glue(
  "There are {mets[, uniqueN(person_id)]} people with metformin data."
))

sub3 <- merge.data.table(
  sub3,
  ins,
  by = "person_id",
  all.x = TRUE
)
sub3 <- merge.data.table(
  sub3,
  mets,
  by = "person_id",
  all.x = TRUE
)
sub3[is.na(insulin), insulin := 0L]
sub3[is.na(metformin), metformin := 0L]

fwrite(
  sub3,
  file.path(processed_directory, glue("cohort_{version}.csv"))
)

cli_alert_success("Cohort identified and saved. 🥳")
