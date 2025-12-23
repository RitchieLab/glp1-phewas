# libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
    library(arrow)
    library(cli)
    library(glue)
    library(data.table)
    library(PheWAS)
    library(dplyr)
    library(comorbidity)
    library(MatchIt)
    library(survRM2)
    library(survival)
    library(survminer)
    library(furrr)
    library(janitor)
    library(optparse)
    library(coxphf)
})

# options ----------------------------------------------------------------------
option_list <- list(
    make_option(
        c("-p", "--project_directory"),
        type = "character",
        default = "/home/jupyter/workspaces/duplicateofglp1rapharmacogenomics",
        help = "Path to project directory [default = %default]"
    ),
    make_option(
        c("-v", "--version"),
        type = "character",
        default = "20250129",
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
        default = "logit",
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
options(datatable.print.class = TRUE)
options(datatable.print.trunc.cols = TRUE)

project_directory <- opt$project_directory
data_directory <- file.path(project_directory, "data")
raw_directory <- file.path(data_directory, "raw")
processed_directory <- file.path(data_directory, "processed")
results_directory <- file.path(project_directory, "results")
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
    "encounters",
    # lab counts
    "dbp",
    "glucose",
    "hba1c",
    "egfr",
    "crp",
    "ast"
)

# functions --------------------------------------------------------------------
## matching_patterns: identify variables that match a pattern
matching_patterns <- function(patterns, character_vector) {
    patterns[base::sapply(
        patterns,
        \(pattern) base::any(grepl(pattern, character_vector))
    )]
}

## rmst_diff: calculate restricted mean survival time differences
rmst_diff <- function(
    matched_data,
    time_var = "time",
    outcome_var = "outcome",
    .trt_var = "trt_drug",
    .taus = 1:3
) {
    data.table::rbindlist(
        lapply(
            .taus,
            function(i) {
                res <- survRM2::rmst2(
                    matched_data[[time_var]],
                    matched_data[[outcome_var]],
                    matched_data[[.trt_var]],
                    tau = i * 365
                )$unadjusted.result |>
                    data.table::as.data.table(keep.rownames = TRUE)
                res <- res[res[[1]] == "RMST (arm=1)-(arm=0)", ]
                data.table::setDT(res)
                res$time <- i
                res[]
            }
        )
    )
}

## coxph_extract: extract Cox PH model results
coxph_extract <- function(x) {
    x_sum <- summary(x)
    tmp_or <- x_sum$conf.int
    data.table::data.table(
        beta = x_sum$coefficients[1, 1],
        or = tmp_or[1, 1],
        or_lo = tmp_or[1, 3],
        or_hi = tmp_or[1, 4],
        p = x_sum$coefficients[1, 5]
    )
}

## coxphf_extract: extract Cox PH Firth model results
coxphf_extract <- function(x) {
    data.table::data.table(
        beta = x$coefficients[1],
        or = exp(x$coefficients[1]),
        or_lo = x$ci.lower[1],
        or_hi = x$ci.upper[1],
        p = x$prob[1]
    )
}

## transform_rmst_results: transform RMST difference results
transform_rmst_results <- function(x) {
    x |>
        select(
            time,
            rmst_est = Est.,
            rmst_lo = `lower .95`,
            rmst_hi = `upper .95`,
            rmst_p = p
        ) |>
        pivot_longer(
            cols = -time,
            names_to = "metric",
            values_to = "value"
        ) |>
        mutate(metric = paste0(metric, "_", time)) |>
        select(-time) |>
        pivot_wider(
            names_from = metric,
            values_from = value
        )
}

## fit_cox_model: fit Cox PH and, if necessary, Firth model
fit_cox_model <- function(formula, data) {
    # Initialize storage for models
    cox_model <- NULL
    cox_firth_model <- NULL
    warning_flag <- FALSE

    # Try fitting the Cox PH model
    cox_model <- tryCatch(
        {
            withCallingHandlers(
                coxph(formula, data = data, x = TRUE),
                warning = function(w) {
                    if (
                        grepl("loglik converg", w$message, ignore.case = TRUE)
                    ) {
                        message(
                            "Warning: Log-likelihood convergence issue detected in `coxph()`. Also fitting `coxphf()`."
                        )
                        warning_flag <<- TRUE # Set flag to indicate a warning occurred
                    }
                }
            )
        },
        error = function(e) {
            if (
                grepl(
                    "system is computationally singular",
                    e$message,
                    ignore.case = TRUE
                )
            ) {
                message(
                    "Error: Singular matrix detected in `coxph()`. Switching to `coxphf()` for penalized estimation."
                )
                return(NULL) # Return NULL to indicate failure
            }
            stop(e) # If another error, stop execution
        }
    )

    # If Cox PH model had a warning or failed, try fitting the Firth model
    if (is.null(cox_model) || warning_flag) {
        cox_firth_model <- tryCatch(
            {
                suppressWarnings(
                    coxphf(
                        formula,
                        data = data,
                        maxit = 500,
                        maxstep = 0.25,
                        pl = FALSE
                    )
                )
            },
            error = function(e) {
                if (
                    grepl(
                        "system is computationally singular",
                        e$message,
                        ignore.case = TRUE
                    )
                ) {
                    message(
                        "Error: Singular matrix detected in `coxphf()`. Returning NULL."
                    )
                    return(NULL) # Return NULL if Firth model also fails
                }
                stop(e) # If another error, stop execution
            }
        )
    }

    # Return both models and warning flag
    return(list(
        coxph = cox_model,
        coxphf = cox_firth_model,
        warning_in_coxph = warning_flag
    ))
}

## cox_extractor: extract Cox PH or Firth model results depending on which was fit
cox_extractor <- function(x) {
    if (!is.null(x[["coxphf"]])) {
        return(coxphf_extract(x[["coxphf"]])[, firth := TRUE])
    } else {
        return(coxph_extract(x[["coxph"]])[, firth := FALSE])
    }
}

## filter_predictors: filter predictors based on near-zero variance and high correlation
filter_predictors <- function(data, corr_cutoff = 0.9, ignore.vars = c()) {
    # Remove near-zero variance variables
    nzv <- caret::nearZeroVar(data, saveMetrics = TRUE)
    nzv_vars <- as.data.table(nzv, keep.rownames = "var")[nzv == TRUE, var]
    nzv_vars <- setdiff(nzv_vars, ignore.vars)

    data <- data |> dplyr::select(-tidyselect::any_of(nzv_vars))

    # Remove highly correlated variables
    num_vars <- names(data)[sapply(data, is.numeric)]
    cor_matrix <- cor(data[, ..num_vars], use = "pairwise.complete.obs")

    # Find columns to remove
    highly_correlated <- caret::findCorrelation(
        cor_matrix,
        cutoff = corr_cutoff,
        names = TRUE
    )
    hc_vars <- setdiff(highly_correlated, ignore.vars)

    data |> dplyr::select(-tidyselect::any_of(hc_vars))
}

## fn.cec: calculate case and event counts
fn.cec <- function(data, outcome_var = "outcome_itt") {
    data.table(
        n_cases = data[trt_drug == 1, .N],
        n_controls = data[trt_drug == 0, .N],
        n_case_events = data[trt_drug == 1 & get(outcome_var) == 1, .N],
        n_control_events = data[trt_drug == 0 & get(outcome_var) == 1, .N]
    )
}

## safe_cox_zph: run cox.zph with error handling
safe_cox_zph <- function(cox_model) {
    # Try running cox.zph, catch singular matrix errors
    cox_zph_result <- tryCatch(
        {
            cox.zph(cox_model) # Run proportional hazards test
        },
        error = function(e) {
            if (
                grepl(
                    "system is computationally singular",
                    e$message,
                    ignore.case = TRUE
                )
            ) {
                message(
                    "Error: Singular matrix detected in `cox.zph()`. Returning NA."
                )
                return(NULL) # Return NULL instead of crashing
            }
            message("Unexpected error in `cox.zph()`. Stopping execution.")
            stop(e) # Stop execution only for unknown errors
        }
    )

    return(cox_zph_result)
}

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
    yes <- dt[search_term == "semaglutide", unique(person_id)]
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

# load ccsr data
phe <- open_dataset(
    glue("{processed_directory}/CCSR_{version}.parquet")
) |>
    filter(person_id %in% dt[, person_id]) |>
    collect()
setDT(phe)
setnames(phe, "ccsr", "phecode")

setkey(dt, NULL)
setkey(phe, NULL)

phe2 <- data.table::merge.data.table(
    phe,
    dt[, .(person_id, drug_date)],
    by = "person_id"
)
phe2[, uniqueN(person_id)]

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
        return(data.table())
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
        itt_res <- fit_cox_model(itt_formula, data = matched_data)
        itt_coxph_test <- safe_cox_zph(itt_res[["coxph"]])

        # fit Cox PH per-protocol model
        pp_res <- fit_cox_model(pp_formula, data = matched_data)
        pp_coxph_test <- safe_cox_zph(pp_res[["coxph"]])

        # calculate RMST differences
        itt_rmst_diff <- rmst_diff(
            matched_data,
            .trt_var = trt_var,
            time_var = "time_itt",
            outcome_var = "outcome_itt",
            .taus = taus
        )
        pp_rmst_diff <- rmst_diff(
            matched_data,
            .trt_var = trt_var,
            time_var = "time_pp",
            outcome_var = "outcome_pp",
            .taus = taus
        )
        ### case and event counts
        cec_vars_itt <- c("trt_drug", "outcome_itt")
        cec_vars_pp <- c("trt_drug", "outcome_pp")
        if (length(adj_vars) > 0) {
            cec_vars_itt <- c(cec_vars_itt, adj_vars)
            cec_vars_pp <- c(cec_vars_pp, adj_vars)
        }
        cca_itt <- matched_data[
            complete.cases(matched_data[, ..cec_vars_itt]),
        ]
        cca_pp <- matched_data[complete.cases(matched_data[, ..cec_vars_pp]), ]

        # extract cox PH test p-values
        if (!is.null(itt_coxph_test)) {
            itt_coxph_test_p <- itt_coxph_test$table["GLOBAL", "p"]
        } else {
            itt_coxph_test_p <- NA
        }
        if (!is.null(pp_coxph_test)) {
            pp_coxph_test_p <- pp_coxph_test$table["GLOBAL", "p"]
        } else {
            pp_coxph_test_p <- NA
        }

        # create summary table
        quick_summary <- data.table::rbindlist(list(
            cbind(
                cox_extractor(itt_res)[, `:=`(
                    type = "ITT",
                    ph_test_p = itt_coxph_test_p
                )],
                fn.cec(cca_itt),
                transform_rmst_results(itt_rmst_diff)
            ),
            cbind(
                cox_extractor(pp_res)[, `:=`(
                    type = "PP",
                    ph_test_p = pp_coxph_test_p
                )],
                fn.cec(cca_pp, outcome_var = "outcome_pp"),
                transform_rmst_results(pp_rmst_diff)
            )
        ))[, adj_vars := paste0(adj_vars, collapse = ", ")]

        return(quick_summary[])
    }
}


dt1 <- copy(dt)
dt1[, drug_age := round(as.numeric(drug_date - t2d_date) / 365.25, 1)]
phe_sub <- phe_post[person_id %in% dt1[, person_id], ]
phecodes <- unique(phe_sub[, .(phecode, person_id)])[, .N, phecode][
    order(-N),
][N > 100, phecode]

# ## test
# trt_var <- "trt_drug"
# time_var <- "time"
# event_var <- "outcome"
# phewas_code <- "EM_236"
# phecode_var <- "phecode"
# taus <- 1:3
# analysis_function(
#     data = dt1,
#     phe_pre = phe_pre,
#     phe_post = phe_post,
#     psm_covs = psm_covariates,
#     phewas_code = phecodes[710],
#     trt_var = "trt_drug",
#     event_var = "outcome",
#     phecode_var = "phecode",
#     .max_date_var = max_date_var,
#     distance_approach = "glm",
#     cor_cut = 0.85,
#     taus = 1:2,
#     caliper_width = 0.25,
#     filter_pred = TRUE
# ) |>
#     glimpse()

results2 <- future_map(
    phecodes, # phecodes with more than 100 occurrences
    \(.phecode) {
        analysis_function(
            data = dt1,
            phe_pre = phe_pre,
            phe_post = phe_post,
            psm_covs = psm_covariates,
            phewas_code = .phecode,
            trt_var = "trt_drug",
            event_var = "outcome",
            phecode_var = "phecode",
            .max_date_var = max_date_var,
            taus = 1:3,
            distance_approach = dist_approach,
            caliper_width = match_cal_width
        )[, phecode := .phecode]
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
)

results3 <- results2 |>
    rbindlist(use.names = TRUE, fill = TRUE) |>
    select(phecode, everything())

fwrite(
    results3,
    file.path(
        results_directory,
        glue(
            "{make_clean_names(treatment)}_{make_clean_names(comparator)}_ccsrwas_{version}.csv"
        )
    )
)

cli_alert_success("Analysis complete! 🥳")
