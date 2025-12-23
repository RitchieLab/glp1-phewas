# libraries --------------------------------------------------------------------
source("R/config.R")
source("R/setup.R")
source("R/fn/helper-functions.R")
load_libraries(c(
  pkg_groups$core,
  pkg_groups$analysis,
  pkg_groups$utils
))

outcomes <- readLines("tables/outcomes_for_stratification")
set.seed(123)

# options ----------------------------------------------------------------------
option_list <- list(
  make_option(
    c("-p", "--project_directory"),
    type = "character",
    default = config$project_directory,
    help = "Path to project directory [default = %default]"
  ),
  make_option(
    c("-v", "--version"),
    type = "character",
    default = config$version,
    help = "Data version [default = %default]"
  ),
  make_option(
    c("-t", "--treatment"),
    type = "character",
    default = "GLP-1 RA",
    help = "Treatment drug, either 'GLP-1 RA' or 'Semaglutide' [default = %default]"
  ),
  make_option(
    c("-c", "--comparator"),
    type = "character",
    default = "SGLT2i",
    help = "Comparator drug, either 'SGLT2i' or 'DPP4i' [default = %default]"
  ),
  make_option(
    c("-l", "--match_cal_width"),
    type = "numeric",
    default = "0.1",
    help = "Caliper width for matching, usually 0.1 or 0.25 [default = %default]"
  ),
  make_option(
    c("-s", "--seed"),
    type = "numeric",
    default = "4571",
    help = "Seed for reproducibility [default = %default]"
  ),
  make_option(
    c("-d", "--distance_approach"),
    type = "character",
    default = "lasso",
    help = "Model for calculating propensity score [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

set.seed(opt$seed)
plan("multicore", workers = 8)
options(future.globals.maxSize = 500 * 1024^2)

project_directory <- opt$project_directory
data_directory <- file.path(project_directory, config$path$data)
raw_directory <- file.path(project_directory, config$paths$raw)
processed_directory <- file.path(project_directory, config$paths$processed)
results_directory <- file.path(project_directory, config$paths$results)
version <- opt$version
treatment <- opt$treatment
comparator <- opt$comparator
match_cal_width <- opt$match_cal_width
dist_approach <- opt$distance_approach

if (
  !(treatment %in% c("GLP-1 RA", "Semaglutide")) ||
    !(comparator %in% c("SGLT2i", "DPP4i"))
) {
  stop(
    "Invalid treatment ('GLP-1 RA' or 'Semaglutide') or comparator ('SGLT2i' or 'DPP4i')"
  )
}

max_date_var <- fcase(
  treatment == "GLP-1 RA" & comparator == "SGLT2i",
  "max_glp1_sglt2_date",
  treatment == "GLP-1 RA" & comparator == "DPP4i",
  "max_glp1_dpp4_date",
  treatment == "Semaglutide" & comparator == "SGLT2i",
  "max_sema_sglt2_date",
  treatment == "Semaglutide" & comparator == "DPP4i",
  "max_sema_dpp4_date"
)

dir.create(processed_directory, showWarnings = FALSE)
dir.create(results_directory, showWarnings = FALSE)

psm_covariates <- c(
  # sociodemographics
  "sex",
  "race_eth",
  "t2d_age",
  "drug_age",
  "bmi_last_cat",
  "drug_date",
  "metformin",
  "insulin",
  # neighborhood-level characteristics
  "high_school_education_cat",
  "no_health_insurance_cat",
  "median_income_cat",
  # health history (components of CCI)
  "mi",
  "chf",
  "pvd",
  "cevd",
  "dementia",
  "cpd",
  "rheumd",
  "pud",
  "mld",
  "hp",
  "rend",
  "canc",
  "msld",
  "metacanc",
  "aids",
  "diab",
  "diabwc",
  # healthcare utilization
  "followup",
  "n_visits",
  # lab counts
  "dbp",
  "glucose",
  "hba1c",
  "egfr",
  "crp",
  "ast"
)

# load data --------------------------------------------------------------------
dt <- read_parquet(
  file.path(
    processed_directory,
    paste0("analytic_dataset_", version, ".parquet")
  )
)
setDT(dt)

# initialize outcome variable
if (treatment == "GLP-1 RA") {
  yes <- dt[drug_class == "GLP-1 RA", unique(person_id)]
  no <- dt[drug_class == comparator, unique(person_id)]
} else if (treatment == "Semaglutide") {
  yes <- dt[drug_ingredient_name == "semaglutide", unique(person_id)]
  no <- dt[
    !is.na(get(max_date_var)) &
      !(person_id %in% yes) &
      drug_class == comparator,
    unique(person_id)
  ]
} else {
  stop("Invalid treatment - must be 'GLP-1 RA' or 'Semaglutide'")
}

# create treatment variable
dt[,
  trt_drug := fcase(
    person_id %in% yes,
    1,
    person_id %in% no,
    0,
    default = NA_real_
  )
]
dt <- dt[!is.na(trt_drug), ]

# load phenome data
phe <- open_dataset(
  glue("{processed_directory}/phenome_{version}.parquet")
) |>
  filter(person_id %in% dt[, person_id]) |>
  collect()
setDT(phe)

setkey(dt, NULL)
setkey(phe, NULL)

phe2 <- data.table::merge.data.table(
  phe,
  dt[, .(person_id, drug_date)],
  by = "person_id"
)

# split phenome data into pre- and post-treatment
phe_pre <- phe2[between(date, round(drug_date - (365.25 * 2)), drug_date), ][,
  drug_date := NULL
]
phe_post <- phe2[date > drug_date, ][, drug_date := NULL]

setkey(phe_pre, phecode)
setkey(phe_post, phecode)

# analysis ---------------------------------------------------------------------
## analysis_function: perform PheWAS analysis
analysis_function <- function(
  data,
  phe_pre,
  phe_post,
  psm_covs = psm_covariates,
  phewas_code = "EM_236", # outome of interest
  trt_var = "trt_drug",
  event_var = "outcome",
  phecode_var = "phecode",
  .max_date_var = "max_glp1_sglt2_date",
  cor_cut = 0.9,
  taus = 1:3,
  caliper_width = 0.1,
  filter_pred = FALSE,
  distance_approach = "lasso",
  mah_vars = NULL
) {
  # identify cases and non-eligible controls
  not_eligible <- phe_pre[.(phewas_code), unique(person_id)]
  cases <- phe_post[.(phewas_code), unique(person_id)]

  # restrict to cases and eliglible controls
  dt2 <- data[!person_id %in% not_eligible, ]

  # create outcome variable
  dt2[, outcome := as.numeric(person_id %in% cases)]

  # subset phe_post to only include the first occurrence of the phecode
  phe_post_sub <- phe_post[
    get(phecode_var) == phewas_code &
      person_id %in% setdiff(cases, not_eligible),
  ]
  phe_post_sub <- phe_post_sub[
    phe_post_sub[, .I[which.min(date)], by = person_id][, V1],
  ]

  # merge in outcome occurrence date to the dataset
  dt3 <- merge(
    dt2,
    phe_post_sub[, .(person_id, phe_date = date)],
    by = "person_id",
    all.x = TRUE
  )

  # create end of follow-up date; phecode occurrence date or 2023-10-01 for intention-to-treat
  # and minimum of phecode occurrence date, maximum GLP-1 date + 90, SGLT2 initiation date, or death date for per-protocol
  dt3[, `:=`(
    max_itt_date = pmin(phe_date, as.Date("2023-10-01"), na.rm = TRUE),
    max_pp_date = pmin(phe_date, get(.max_date_var), na.rm = TRUE)
  )]

  dt3[, `:=`(
    outcome_itt = as.numeric(person_id %in% cases),
    outcome_pp = as.numeric(person_id %in% cases & phe_date == max_pp_date)
  )]

  # create time to end of follow-up variable
  dt3[, `:=`(
    time_itt = as.numeric(max_itt_date - drug_date),
    time_pp = as.numeric(max_pp_date - drug_date)
  )]

  # create matching formula
  match_formula <- paste0(trt_var, " ~ ", paste0(psm_covs, collapse = " + "))

  # match on propensity score
  match_obj <- MatchIt::matchit(
    as.formula(match_formula),
    data = dt3,
    method = "nearest",
    distance = distance_approach,
    caliper = caliper_width,
    ratio = 1,
    replace = FALSE,
    mahvars = mah_vars
  )
  matched_data <- match.data(match_obj)

  if (
    matched_data[, sum(outcome_itt, na.rm = TRUE)] < 50 ||
      matched_data[, sum(outcome_pp, na.rm = TRUE)] < 50
  ) {
    return(data.table(note = "Insufficient cases after matching"))
  } else {
    # identify variables to adjust for in the Cox model based on residual standard mean differences
    adj_vars <- data.table::as.data.table(
      summary(match_obj)$sum.matched,
      keep.rownames = "term"
    )
    setnames(
      adj_vars,
      old = names(adj_vars),
      new = janitor::make_clean_names(names(adj_vars))
    )
    adj_vars <- adj_vars[abs(std_mean_diff) > 0.1, term]
    adj_vars <- matching_patterns(
      patterns = c("glp1", psm_covariates),
      character_vector = adj_vars
    )

    # filter predictors based on near-zero variance and high correlation
    if (filter_pred == TRUE) {
      matched_data <- filter_predictors(
        matched_data,
        corr_cutoff = cor_cut,
        ignore.vars = c(
          "outcome_itt",
          "outcome_pp",
          "time_itt",
          "time_pp",
          "trt_drug"
        )
      )
      adj_vars <- adj_vars[adj_vars %in% names(matched_data)]
    }

    # create Cox PH formulas
    if (length(adj_vars) > 0) {
      covar_string <- paste0(c(trt_var, adj_vars), collapse = " + ")
    } else {
      covar_string <- trt_var
    }

    itt_formula <- as.formula(paste0(
      "Surv(time_itt, outcome_itt) ~ ",
      covar_string
    ))
    pp_formula <- as.formula(paste0(
      "Surv(time_pp, outcome_pp) ~ ",
      covar_string
    ))

    # fit Cox PH intention-to-treat model
    # itt_cox <- survival::coxph(itt_formula, data = matched_data, x = TRUE)
    itt_res <- coxphf::coxphf(
      itt_formula,
      data = matched_data,
      maxit = 500,
      maxstep = 0.25,
      pl = FALSE
    )
    # itt_res <- fit_cox_model(itt_formula, data = matched_data)
    # itt_coxph_test <- safe_cox_zph(itt_cox)

    # fit Cox PH per-protocol model
    # pp_cox <- survival::coxph(pp_formula, data = matched_data, x = TRUE)
    pp_res <- coxphf::coxphf(
      pp_formula,
      data = matched_data,
      maxit = 500,
      maxstep = 0.25,
      pl = FALSE
    )
    # pp_res <- fit_cox_model(pp_formula, data = matched_data)
    # pp_coxph_test <- safe_cox_zph(pp_cox)

    itt_out <- data.table(
      beta = itt_res$coefficients[["trt_drug"]],
      hr = exp(itt_res$coefficients[["trt_drug"]]),
      hr_lo = itt_res$ci.lower[["trt_drug"]],
      hr_hi = itt_res$ci.upper[["trt_drug"]],
      p = itt_res$prob[["trt_drug"]],
      n_obs = itt_res$n,
      n_cases = sum(matched_data[, trt_drug]),
      n_events = sum(matched_data[, outcome_itt]),
      n_case_events = sum(matched_data[trt_drug == 1, outcome_itt]),
      type = "ITT"
    )
    pp_out <- data.table(
      beta = pp_res$coefficients[["trt_drug"]],
      hr = exp(pp_res$coefficients[["trt_drug"]]),
      hr_lo = pp_res$ci.lower[["trt_drug"]],
      hr_hi = pp_res$ci.upper[["trt_drug"]],
      p = pp_res$prob[["trt_drug"]],
      n_obs = pp_res$n,
      n_cases = sum(matched_data[, trt_drug]),
      n_events = sum(matched_data[, outcome_pp]),
      n_case_events = sum(matched_data[trt_drug == 1, outcome_pp]),
      type = "PP"
    )

    return(rbindlist(list(itt_out, pp_out), use.names = TRUE, fill = TRUE))
  }
}

dt1 <- copy(dt)
dt1[, drug_age := round(as.numeric(drug_date - t2d_date) / 365.25, 1)]
phe_sub <- phe_post[person_id %in% dt1[, person_id], ]

psm_covariates2 <- c(
  "sex",
  "race_eth",
  "t2d_age",
  "drug_age",
  "bmi_last_cat",
  "drug_date",
  "metformin",
  "insulin",
  "cci_score"
)

res_male <- rbindlist(
  lapply(
    outcomes,
    \(phe_out) {
      cli_progress_step("{phe_out}")
      analysis_function(
        data = dt1[sex == "Male", ],
        phe_pre = phe_pre,
        phe_post = phe_post,
        psm_covs = psm_covariates2[!psm_covariates2 %in% "sex"],
        phewas_code = phe_out,
        trt_var = "trt_drug",
        event_var = "outcome",
        phecode_var = "phecode",
        .max_date_var = max_date_var,
        taus = 1:3,
        distance_approach = dist_approach,
        caliper_width = match_cal_width
      )[, phecode := phe_out]
    }
  ),
  use.names = TRUE,
  fill = TRUE
)[, stratum := "Male"][]

res_female <- rbindlist(
  lapply(
    outcomes,
    \(phe_out) {
      cli_progress_step("{phe_out}")
      analysis_function(
        data = dt1[sex == "Female", ],
        phe_pre = phe_pre,
        phe_post = phe_post,
        psm_covs = psm_covariates2[!psm_covariates2 %in% "sex"],
        phewas_code = phe_out,
        trt_var = "trt_drug",
        event_var = "outcome",
        phecode_var = "phecode",
        .max_date_var = max_date_var,
        taus = 1:3,
        distance_approach = dist_approach,
        caliper_width = match_cal_width
      )[, phecode := phe_out]
    }
  ),
  use.names = TRUE,
  fill = TRUE
)[, stratum := "Female"][]

res_nhw <- rbindlist(
  lapply(
    outcomes,
    \(phe_out) {
      cli_progress_step("{phe_out}")
      analysis_function(
        data = dt1[race_eth == "NHW", ],
        phe_pre = phe_pre,
        phe_post = phe_post,
        psm_covs = psm_covariates2[!psm_covariates2 %in% "race_eth"],
        phewas_code = phe_out,
        trt_var = "trt_drug",
        event_var = "outcome",
        phecode_var = "phecode",
        .max_date_var = max_date_var,
        taus = 1:3,
        distance_approach = dist_approach,
        caliper_width = match_cal_width
      )[, phecode := phe_out]
    }
  ),
  use.names = TRUE,
  fill = TRUE
)[, stratum := "NHW"][]

res_nhb <- rbindlist(
  lapply(
    outcomes,
    \(phe_out) {
      cli_progress_step("{phe_out}")
      tryCatch(
        {
          analysis_function(
            data = dt1[race_eth == "NHB", ],
            phe_pre = phe_pre,
            phe_post = phe_post,
            psm_covs = psm_covariates2[!psm_covariates2 %in% "race_eth"],
            phewas_code = phe_out,
            trt_var = "trt_drug",
            event_var = "outcome",
            phecode_var = "phecode",
            .max_date_var = max_date_var,
            taus = 1:3,
            distance_approach = dist_approach,
            caliper_width = 0.25
          )[, phecode := phe_out]
        },
        error = function(e) data.table()
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)[, stratum := "NHB"][]

res_hisp <- rbindlist(
  lapply(
    outcomes,
    \(phe_out) {
      cli_progress_step("{phe_out}")
      tryCatch(
        {
          analysis_function(
            data = dt1[race_eth == "HISP", ],
            phe_pre = phe_pre,
            phe_post = phe_post,
            psm_covs = psm_covariates2[!psm_covariates2 %in% "race_eth"],
            phewas_code = phe_out,
            trt_var = "trt_drug",
            event_var = "outcome",
            phecode_var = "phecode",
            .max_date_var = max_date_var,
            taus = 1:3,
            distance_approach = dist_approach,
            caliper_width = 0.25
          )[, phecode := phe_out]
        },
        error = function(e) data.table()
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)[, stratum := "HISP"][]

res_comb <- rbindlist(
  list(res_male, res_female, res_nhw, res_nhb, res_hisp),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  res_comb,
  file.path(
    results_directory,
    glue(
      "{make_clean_names(treatment)}_{make_clean_names(comparator)}_stratified_{version}.csv"
    )
  )
)

cli_alert_success("Analysis complete! 🥳")
