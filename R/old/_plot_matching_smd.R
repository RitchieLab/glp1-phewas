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
    ),
    make_option(
        c("--phewas_code"),
        type = "character",
        default = "MB_293",
        help = "Phecode for plotting matched SMD [default = %default]"
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
.phewas_code <- opt$phewas_code

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
plot_match_msd <- function(
    match.obj,
    .treatment,
    .comparator,
    phewas_code = .phewas_code
) {
    test <- summary(match.obj)

    all <- as.data.table(test$sum.all, keep.rownames = "variable")
    matc <- as.data.table(test$sum.matched, keep.rownames = "variable")
    all[, sample := "Unmatched"]
    matc[, sample := "Matched"]

    stack <- rbindlist(list(all, matc), use.names = TRUE, fill = TRUE)
    stack <- clean_names(stack)
    stack[, sample := factor(sample, levels = c("Unmatched", "Matched"))]

    stack[,
        variable := factor(
            variable,
            levels = c("distance", setdiff(variable, "distance"))
        )
    ]
    stack[, variable := factor(variable, levels = rev(unique(variable)))]

    cols <- c(
        "Unmatched" = "#D55E00",
        "Matched" = "#009E73"
    )

    ss_fn <- function(y) {
        y[["nn"]] |>
            as.data.table(keep.rownames = "var") |>
            filter(var %in% c("All", "Matched")) |>
            mutate(total = Control + Treated) |>
            (\(x) {
                paste0(
                    "N: ",
                    format(x[var == "All", total], big.mark = ","),
                    " unmatched, ",
                    format(x[var == "Matched", total], big.mark = ","),
                    " matched"
                )
            })()
    }

    plot <- stack |>
        ggplot(aes(
            x = variable,
            y = abs(std_mean_diff),
            fill = sample,
            color = sample
        )) +
        geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey") +
        geom_vline(xintercept = stack[, .N / 2] - 0.5, color = "black") +
        geom_point(size = 3, alpha = 0.75) +
        coord_flip() +
        ylim(0, 1) +
        labs(
            subtitle = paste0(
                "Standardized mean differences before and after matching\n",
                ss_fn(test)
            ),
            title = paste0(
                .treatment,
                " vs. ",
                .comparator,
                " for ",
                phewas_code
            ),
            x = "",
            y = "Absolute standardized mean difference",
        ) +
        scale_fill_manual(values = cols) +
        scale_color_manual(values = cols) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0, face = "bold"),
            legend.position = "top",
            legend.title = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank()
        )
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
# data = dt1
# psm_covs = psm_covariates
# phewas_code = "MB_293"
# trt_var = "trt_drug"
# event_var = "outcome"
# phecode_var = "phecode"
# .max_date_var = "max_glp1_sglt2_date"
# cor_cut = 0.9
# taus = 1:3
# caliper_width = 0.1
# filter_pred = FALSE
# distance_approach = "logit"
# mah_vars = NULL

quick_match <- function(
    data,
    phe_pre,
    phe_post,
    psm_covs = psm_covariates,
    phewas_code = .phewas_code, # outome of interest
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
    return(match_obj)
}

dt1 <- copy(dt)
dt1[, drug_age := round(as.numeric(drug_date - t2d_date) / 365.25, 1)]

result <- quick_match(
    data = dt1,
    phe_pre = phe_pre,
    phe_post = phe_post,
    psm_covs = psm_covariates,
    phewas_code = .phewas_code,
    trt_var = "trt_drug",
    event_var = "outcome",
    phecode_var = "phecode",
    .max_date_var = max_date_var,
    taus = 1:3,
    distance_approach = dist_approach,
    caliper_width = match_cal_width
)

smd_plot <- plot_match_msd(
    result,
    .treatment = treatment,
    .comparator = comparator
)

ggsave(
    plot = smd_plot,
    filename = paste0(
        "plots/",
        make_clean_names(treatment),
        "_",
        make_clean_names(comparator),
        "_",
        make_clean_names(.phewas_code),
        "_",
        "match_plot.pdf"
    ),
    device = cairo_pdf,
    width = 8,
    height = 11,
    units = "in"
)

saveRDS(
    object = smd_plot,
    file = paste0(
        "plots/",
        make_clean_names(treatment),
        "_",
        make_clean_names(comparator),
        "_",
        make_clean_names(.phewas_code),
        "_",
        "match_plot.rds"
    )
)


cli_alert_success("Analysis complete! 🥳")
