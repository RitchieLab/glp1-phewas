# libraries --------------------------------------------------------------------
source("R/config.R")
source("R/setup.R")
source("R/fn/helper-functions.R")
load_libraries(c(
  pkg_groups$core,
  pkg_groups$analysis,
  pkg_groups$summarize,
  pkg_groups$utils
))

project_directory <- config$project_directory
data_directory <- file.path(project_directory, config$paths$data)
raw_directory <- file.path(project_directory, config$paths$raw)
processed_directory <- file.path(project_directory, config$paths$processed)
tables_directory <- file.path(project_directory, config$paths$tables)
version <- config$version

dir.create(processed_directory, showWarnings = FALSE)

# cohort -----------------------------------------------------------------------
cohort <- fread(file.path(processed_directory, glue("cohort_{version}.csv")))
cohort[, person_id := as.character(person_id)]
cohort_vars <- names(cohort)
for (i in cohort_vars) {
  if ("IDate" %in% class(cohort[[i]])) {
    cohort[[i]] <- as.Date(cohort[[i]])
  }
}
setkey(cohort, person_id)

# demographics -----------------------------------------------------------------
demo <- open_dataset(file.path(raw_directory, glue("demo_{version}"))) |>
  mutate(person_id = as.character(person_id)) |>
  filter(person_id %in% cohort[, person_id]) |>
  select(
    person_id,
    gender,
    dob = date_of_birth,
    race,
    ethnicity,
    sex_at_birth
  ) |>
  collect()
setDT(demo)
setkey(demo, person_id)

demo[, `:=`(
  dob = as.Date(dob),
  race = fifelse(
    race %in%
      c(
        "None Indicated",
        "I prefer not to answer",
        "None of these",
        "PMI: Skip"
      ),
    "Unknown",
    race
  ),
  ethnicity = fifelse(
    ethnicity %in%
      c(
        "No matching concept",
        "PMI: Prefer Not To Answer",
        "What Race Ethnicity: Race Ethnicity None Of These",
        "PMI: Skip"
      ),
    "Unknown",
    ethnicity
  ),
  sex = fcase(
    sex_at_birth == "Female",
    "Female",
    sex_at_birth == "Male",
    "Male"
  )
)]

demo[, `:=`(
  race_eth = fcase(
    race == "White" & ethnicity == "Not Hispanic or Latino",
    "NHW",
    race == "Black or African American" &
      ethnicity == "Not Hispanic or Latino",
    "NHB",
    race == "Asian" & ethnicity == "Not Hispanic or Latino",
    "NHA",
    ethnicity == "Hispanic or Latino",
    "HISP",
    default = "Other/Unknown"
  )
)]

demo <- demo[, .(person_id, sex, race_eth)]

# bmi --------------------------------------------------------------------------
bmi <- open_dataset(file.path(raw_directory, glue("bmi_{version}"))) |>
  mutate(person_id = as.character(person_id)) |>
  filter(
    person_id %in%
      cohort[, person_id] &
      measurement_concept_id %in% "3038553"
  ) |>
  collect()
setDT(bmi)
setkey(bmi, person_id)

bmi[, `:=`(
  date = as.Date(measurement_datetime),
  bmi = value_as_number
)]
bmi <- unique(bmi[!is.na(bmi), .(person_id, bmi_last = bmi, date)])

bmi <- merge(
  bmi,
  cohort[, .(person_id, drug_date)],
  by = "person_id"
)

bmi <- bmi[date <= drug_date, .SD[which.max(date)], person_id]

bmi[, `:=`(
  bmi_last_cat = fcase(
    bmi_last < 18.5,
    "Underweight",
    bmi_last >= 18.5 & bmi_last < 25,
    "Normal",
    bmi_last >= 25 & bmi_last < 30,
    "Overweight",
    bmi_last >= 30,
    "Obese"
  )
)]
bmi <- bmi[, .(person_id, bmi_last, bmi_last_date = date, bmi_last_cat)]
bmi <- merge(bmi, cohort[, .(person_id)], by = "person_id", all = TRUE)
bmi[is.na(bmi_last_cat), bmi_last_cat := "No record"]

# ses --------------------------------------------------------------------------
ses <- open_dataset(file.path(raw_directory, glue("ses_{version}"))) |>
  mutate(person_id = as.character(person_id)) |>
  filter(person_id %in% cohort[, person_id]) |>
  rename(
    date = observation_datetime
  ) |>
  collect()
setDT(ses)
setkey(ses, person_id)

ses[, date := as.Date(date)]

ses <- ses[, .(
  person_id,
  high_school_education,
  no_health_insurance,
  median_income,
  deprivation_index
)]
ses <- merge(
  ses,
  cohort[, .(person_id)],
  by = "person_id",
  all = TRUE
)
ses_vars <- names(ses)[names(ses) != "person_id"]
for (i in ses_vars) categorize_quartiles(ses, i)

# cci score --------------------------------------------------------------------
icd <- open_dataset(file.path(raw_directory, glue("icd_{version}"))) |>
  mutate(person_id = as.character(person_id)) |>
  filter(person_id %in% cohort[, person_id]) |>
  select(
    person_id,
    code = source_concept_code,
    vocabulary_id = source_vocabulary,
    date = condition_start_datetime
  ) |>
  collect()
setDT(icd)
icd[, date := as.Date(date)]
setkey(icd, person_id)

icd <- merge(
  icd,
  cohort[, .(person_id, drug_date)],
  by = "person_id"
)
icd <- icd[date <= drug_date, ]

cci <- cci_score(icd)

cci <- merge(cci, cohort[, .(person_id)], by = "person_id", all = TRUE)
cci[is.na(cci)] <- 0
cci[, `:=`(
  cci_score_cat = fcase(
    cci_score == 0,
    "Very low",
    cci_score %in% 1:2,
    "Low",
    cci_score %in% 3:5,
    "Moderate",
    cci_score %in% 5:10,
    "High",
    cci_score > 10,
    "Very high"
  )
)]

# labs -------------------------------------------------------------------------
lab_maps <- fread(
  file.path(tables_directory, "omop_mappings.csv"),
  na.strings = c("", "NA")
)[group == "lab", ][, `:=`(
  abbrev = fcase(
    study_name == "Glucose",
    "glucose",
    study_name == "HbA1c",
    "hba1c",
    study_name == "Creatinine",
    "creatinine",
    study_name == "eGFR",
    "egfr",
    study_name == "Cholesterol",
    "cholesterol",
    study_name == "AST",
    "ast",
    study_name == "ALT",
    "alt",
    study_name == "C-reactive protein",
    "crp",
    study_name == "Amylase",
    "amylase",
    study_name == "LDL",
    "ldl",
    study_name == "HDL",
    "hdl",
    study_name == "Albumin/creatinine ratio",
    "acr",
    study_name == "Diastolic blood pressure",
    "dbp",
    study_name == "Systolic blood pressure",
    "sbp"
  )
)][]
lab <- open_dataset(file.path(raw_directory, glue("lab_{version}"))) |>
  mutate(person_id = as.character(person_id)) |>
  select(
    person_id,
    concept_id = measurement_concept_id,
    concept_name = standard_concept_name,
    date = measurement_datetime
  ) |>
  filter(
    person_id %in%
      cohort[, person_id] &
      concept_id %in% lab_maps[, concept_id]
  ) |>
  collect()
setDT(lab)
lab[, date := as.Date(date)]
setkey(lab, person_id)

lab <- merge(
  lab,
  cohort[, .(person_id, drug_date)],
  by = "person_id"
)

lab <- unique(lab[
  data.table::between(date, drug_date - (365.25 * 2), drug_date),
  .(person_id, concept_id, concept_name, date)
])
lab <- lab[, .N, .(person_id, concept_id)]

lab_count <- merge(
  lab,
  lab_maps[, .(concept_id, lab = abbrev)],
  by = "concept_id",
  all = TRUE
)
lab_count <- dcast(
  lab_count[, !c("concept_id")],
  person_id ~ lab,
  value.var = "N",
  fill = 0
)
lab_count <- merge(
  lab_count,
  cohort[, .(person_id)],
  by = "person_id",
  all = TRUE
)
lab_count[is.na(lab_count)] <- 0

lab_vars <- names(lab_count)[names(lab_count) != "person_id"]
for (i in lab_vars) {
  set(
    lab_count,
    j = paste0(i, "_cat"),
    value = factor(
      fcase(
        lab_count[[i]] == 0,
        "None",
        lab_count[[i]] %in% 1:3,
        "One to three",
        lab_count[[i]] > 3,
        "More than three"
      ),
      levels = c("None", "One to three", "More than three")
    )
  )
}

# merge ------------------------------------------------------------------------
analytic_dataset <- Reduce(
  \(x, y) merge(x, y, by = "person_id", all = TRUE),
  list(
    cohort,
    demo,
    bmi,
    ses,
    cci,
    lab_count
  )
)

analytic_dataset[,
  high_school_education_cat := low_mid_high(high_school_education_cat)
]
analytic_dataset[,
  no_health_insurance_cat := low_mid_high(no_health_insurance_cat)
]
analytic_dataset[, median_income_cat := low_mid_high(median_income_cat)]

analytic_dataset[, `:=`(
  cci_score_cat = factor(
    cci_score_cat,
    levels = c("Very low", "Low", "Moderate", "High", "Very high")
  ),
  bmi_last_cat = factor(
    bmi_last_cat,
    levels = c("Underweight", "Normal", "Overweight", "Obese", "No record")
  ),
  race_eth = factor(
    race_eth,
    levels = c("NHW", "NHB", "HISP", "NHA", "Other/Unknown")
  ),
  t2d_age_cat = factor(
    fcase(
      t2d_age < 45,
      "Under 45",
      t2d_age >= 45 & t2d_age < 55,
      "[45, 55)",
      t2d_age >= 55 & t2d_age < 65,
      "[55, 65)",
      t2d_age >= 65,
      "65 and older"
    ),
    levels = c("Under 45", "[45, 55)", "[55, 65)", "65 and older")
  ),
  followup_cat = factor(
    fcase(
      followup < 5,
      "Less than 5 years",
      followup >= 5 & followup < 10,
      "[5, 10)",
      followup >= 10 & followup < 20,
      "[10, 20)",
      followup >= 20,
      "20 years and more"
    ),
    levels = c(
      "Less than 5 years",
      "[5, 10)",
      "[10, 20)",
      "20 years and more"
    )
  ),
  visits_cat = factor(
    fcase(
      n_visits < 50,
      "Less than 50",
      n_visits >= 50 & n_visits < 100,
      "[50, 100)",
      n_visits >= 100 & n_visits < 200,
      "[100, 200)",
      n_visits >= 200,
      "200 and more"
    ),
    levels = c("Less than 50", "[50, 100)", "[100, 200)", "200 and more")
  ),
  dx_days_cat = factor(
    fcase(
      n_unique_dx_days < 50,
      "Less than 50",
      n_unique_dx_days >= 50 & n_unique_dx_days < 100,
      "[50, 100)",
      n_unique_dx_days >= 100 & n_unique_dx_days < 200,
      "[100, 200)",
      n_unique_dx_days >= 200,
      "200 and more"
    ),
    levels = c("Less than 50", "[50, 100)", "[100, 200)", "200 and more")
  )
)]

# label variables
var_label(analytic_dataset$person_id) <- "Person ID"
var_label(analytic_dataset$drug_date) <- "Date of drug initiation"
var_label(analytic_dataset$dob) <- "Date of birth"
var_label(analytic_dataset$drug_class) <- "Drug class"
var_label(analytic_dataset$t2d_age) <- "Age at T2D diagnosis"
var_label(analytic_dataset$sex) <- "Sex"
var_label(analytic_dataset$race_eth) <- "Race/Ethnicity"
var_label(analytic_dataset$bmi_last_cat) <- "BMI category"
var_label(analytic_dataset$followup) <- "Follow-up time (years)"
var_label(
  analytic_dataset$n_unique_dx_days
) <- "Number of unique days with ICD billing diagnosis"
var_label(
  analytic_dataset$dx_days_cat
) <- "Number of unique days with ICD billing diagnosis, categorized"
var_label(analytic_dataset$n_visits) <- "Number of unique days with a visit"
var_label(
  analytic_dataset$visits_cat
) <- "Number of unique days with a visit, categorized"
var_label(
  analytic_dataset$high_school_education_cat
) <- "High school education (neighborhood)"
var_label(
  analytic_dataset$no_health_insurance_cat
) <- "No health insurance (neighborhood)"
var_label(analytic_dataset$median_income_cat) <- "Median income (neighborhood)"
var_label(
  analytic_dataset$deprivation_index_cat
) <- "Deprivation index (neighborhood)"
var_label(analytic_dataset$mi) <- "Myocardial infarction"
var_label(analytic_dataset$chf) <- "Congestive heart failure"
var_label(analytic_dataset$pvd) <- "Peripheral vascular disorders"
var_label(analytic_dataset$cevd) <- "Cerebrovascular disease"
var_label(analytic_dataset$dementia) <- "Dementia"
var_label(analytic_dataset$cpd) <- "Chronic pulmonary disease"
var_label(analytic_dataset$rheumd) <- "Rheumatoid disease"
var_label(analytic_dataset$pud) <- "Peptic ulcer disease"
var_label(analytic_dataset$mld) <- "Mild liver disease"
var_label(analytic_dataset$diab) <- "Diabetes without complications"
var_label(analytic_dataset$diabwc) <- "Diabetes with complications"
var_label(analytic_dataset$hp) <- "Hemiplegia or paraplegia"
var_label(analytic_dataset$rend) <- "Renal disease"
var_label(analytic_dataset$canc) <- "Cancer (any malignancy)"
var_label(analytic_dataset$msld) <- "Moderate/severe liver disease"
var_label(analytic_dataset$metacanc) <- "Metastatic solid tumor"
var_label(analytic_dataset$aids) <- "HIV/AIDS"
var_label(analytic_dataset$cci_score) <- "Charlson comorbidity index [CCI]"
var_label(analytic_dataset$cci_score_cat) <- "CCI category"
var_label(analytic_dataset$cholesterol) <- "Cholesterol lab (count)"
var_label(analytic_dataset$hba1c) <- "HbA1c lab (count)"
var_label(analytic_dataset$glucose) <- "Glucose lab (count)"
var_label(analytic_dataset$dbp) <- "Diastolic blood pressure (count)"
var_label(analytic_dataset$egfr) <- "eGFR lab (count)"
var_label(analytic_dataset$acr) <- "ACR lab (count)"
var_label(analytic_dataset$alt) <- "ALT lab (count)"
var_label(analytic_dataset$amylase) <- "Amylase lab (count)"
var_label(analytic_dataset$ast) <- "AST lab (count)"
var_label(analytic_dataset$creatinine) <- "Creatinine lab (count)"
var_label(analytic_dataset$crp) <- "CRP lab (count)"
var_label(analytic_dataset$hdl) <- "HDL lab (count)"
var_label(analytic_dataset$ldl) <- "LDL lab (count)"
var_label(analytic_dataset$sbp) <- "Systolic blood pressure (count)"

# generate summary tables ------------------------------------------------------
summary_table <- analytic_dataset[, .(
  drug_class,
  # demo
  t2d_age,
  t2d_age_cat,
  sex,
  race_eth,
  bmi_last_cat,
  # healthcare utilization
  followup,
  followup_cat,
  n_visits,
  visits_cat,
  n_unique_dx_days,
  dx_days_cat,
  # neighborhood
  high_school_education_cat,
  no_health_insurance_cat,
  median_income_cat,
  # comorbidities
  mi,
  chf,
  pvd,
  cevd,
  dementia,
  cpd,
  rheumd,
  pud,
  mld,
  diab,
  diabwc,
  hp,
  rend,
  canc,
  msld,
  metacanc,
  aids,
  cci_score,
  cci_score_cat,
  # labs
  acr,
  alt,
  amylase,
  ast,
  cholesterol,
  creatinine,
  crp,
  dbp,
  egfr,
  glucose,
  hba1c,
  hdl,
  ldl,
  sbp
)] |>
  tab_sum(
    by_var = "drug_class",
    mod_span_head = c("stat_1", "stat_2", "stat_3") ~ "**Drug class**"
  )

demo_table <- analytic_dataset[, .(
  drug_class,
  # demo
  t2d_age,
  t2d_age_cat,
  sex,
  race_eth,
  bmi_last_cat,
  # healthcare utilization
  followup,
  followup_cat,
  n_visits,
  visits_cat,
  n_unique_dx_days,
  dx_days_cat,
  # neighborhood
  high_school_education_cat,
  no_health_insurance_cat,
  median_income_cat
)] |>
  tab_sum(
    by_var = "drug_class",
    mod_span_head = c("stat_1", "stat_2", "stat_3") ~ "**Drug class**"
  )

comorbid_table <- analytic_dataset[, .(
  drug_class,
  # comorbidities
  mi,
  chf,
  pvd,
  cevd,
  dementia,
  cpd,
  rheumd,
  pud,
  mld,
  diab,
  diabwc,
  hp,
  rend,
  canc,
  msld,
  metacanc,
  aids,
  cci_score,
  cci_score_cat
)] |>
  tab_sum(
    by_var = "drug_class",
    mod_span_head = c("stat_1", "stat_2", "stat_3") ~ "**Drug class**"
  )

lab_table <- analytic_dataset[, .(
  drug_class,
  # labs
  acr,
  alt,
  amylase,
  ast,
  cholesterol,
  creatinine,
  crp,
  dbp,
  egfr,
  glucose,
  hba1c,
  hdl,
  ldl,
  sbp,
  # labs - categories
  acr_cat,
  alt_cat,
  amylase_cat,
  ast_cat,
  cholesterol_cat,
  creatinine_cat,
  crp_cat,
  dbp_cat,
  egfr_cat,
  glucose_cat,
  hba1c_cat,
  hdl_cat,
  ldl_cat,
  sbp_cat
)] |>
  tab_sum(
    by_var = "drug_class",
    mod_span_head = c("stat_1", "stat_2", "stat_3") ~ "**Drug class**"
  )

analytic_dataset[, `:=`(
  sema_or_not = fcase(
    drug_ingredient_name == "semaglutide",
    "Semaglutide",
    drug_class == "GLP-1 RA",
    "Non-semaglutide GLP-1 RA",
    default = NA
  ),
  sema_vs_oth = fcase(
    drug_ingredient_name == "semaglutide",
    "Semaglutide",
    drug_class == "SGLT2i",
    "SGLT2i",
    drug_class == "DPP4i",
    "DPP4i",
    default = NA
  )
)]

sema_table <- analytic_dataset[
  !is.na(sema_or_not),
  .(
    sema_or_not,
    # demo
    t2d_age,
    t2d_age_cat,
    sex,
    race_eth,
    bmi_last_cat,
    # healthcare utilization
    followup,
    followup_cat,
    n_visits,
    visits_cat,
    n_unique_dx_days,
    dx_days_cat,
    # neighborhood
    high_school_education_cat,
    no_health_insurance_cat,
    median_income_cat,
    # comorbidities
    mi,
    chf,
    pvd,
    cevd,
    dementia,
    cpd,
    rheumd,
    pud,
    mld,
    diab,
    diabwc,
    hp,
    rend,
    canc,
    msld,
    metacanc,
    aids,
    cci_score,
    cci_score_cat,
    # labs
    acr,
    alt,
    amylase,
    ast,
    cholesterol,
    creatinine,
    crp,
    dbp,
    egfr,
    glucose,
    hba1c,
    hdl,
    ldl,
    sbp
  )
] |>
  tab_sum(
    by_var = "sema_or_not",
    mod_span_head = c("stat_1", "stat_2") ~ "**GLP-1 RA type**"
  )

sema_oth_table <- analytic_dataset[
  !is.na(sema_vs_oth),
  .(
    sema_vs_oth,
    # demo
    t2d_age,
    t2d_age_cat,
    sex,
    race_eth,
    bmi_last_cat,
    # healthcare utilization
    followup,
    followup_cat,
    n_visits,
    visits_cat,
    n_unique_dx_days,
    dx_days_cat,
    # neighborhood
    high_school_education_cat,
    no_health_insurance_cat,
    median_income_cat,
    # comorbidities
    mi,
    chf,
    pvd,
    cevd,
    dementia,
    cpd,
    rheumd,
    pud,
    mld,
    diab,
    diabwc,
    hp,
    rend,
    canc,
    msld,
    metacanc,
    aids,
    cci_score,
    cci_score_cat,
    # labs
    acr,
    alt,
    amylase,
    ast,
    cholesterol,
    creatinine,
    crp,
    dbp,
    egfr,
    glucose,
    hba1c,
    hdl,
    ldl,
    sbp
  )
] |>
  tab_sum(
    by_var = "sema_vs_oth",
    mod_span_head = c("stat_1", "stat_2", "stat_3") ~ "**Drug type**"
  )


gt::gtsave(
  data = as_gt(summary_table),
  filename = glue("summary_table1_{version}.html"),
  path = "tables/"
)

gt::gtsave(
  data = as_gt(demo_table),
  filename = glue("demo_table1_{version}.html"),
  path = "tables/"
)

gt::gtsave(
  data = as_gt(comorbid_table),
  filename = glue("comorbid_table1_{version}.html"),
  path = "tables/"
)

gt::gtsave(
  data = as_gt(lab_table),
  filename = glue("lab_table1_{version}.html"),
  path = "tables/"
)

gt::gtsave(
  data = as_gt(sema_table),
  filename = glue("sema_table1_{version}.html"),
  path = "tables/"
)

gt::gtsave(
  data = as_gt(sema_oth_table),
  filename = glue("sema_vs_oth_table1_{version}.html"),
  path = "tables/"
)

psave(
  analytic_dataset,
  file.path(
    processed_directory,
    glue("analytic_dataset_{version}.parquet")
  )
)

cli_alert_success("processed covariate data 🥳")
