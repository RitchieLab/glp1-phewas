# helper function for querying and saving data from BigQuery to Cloud Storage and then downloading it locally
query_helper <- function(query, data_name, project_path, version) {
  # generate the path to save the data
  query_path <- file.path(
    Sys.getenv("WORKSPACE_BUCKET"),
    "bq_exports",
    version,
    glue("{data_name}_{version}"),
    glue("{data_name}_{version}_*.parquet")
  )

  # save the data to Cloud Storage
  bigrquery::bq_table_save(
    bigrquery::bq_dataset_query(
      Sys.getenv("WORKSPACE_CDR"),
      query,
      billing = Sys.getenv("GOOGLE_PROJECT")
    ),
    query_path,
    destination_format = "parquet",
    compression = "zstd",
    compression_level = 12
  )

  # create the directory to save the data
  data_directory <- file.path(project_path, "data", "raw", data_name)
  dir.create(data_directory, recursive = TRUE, showWarnings = FALSE)

  # download the data from Cloud Storage to the project directory
  system2("gsutil", paste("-m cp ", query_path, data_directory))
}

### pull descendant ids
pull_descendant_ids <- function(concepts) {
  concepts_query <- glue::glue(
    "SELECT DISTINCT
        ca.descendant_id,
        b.name,
        b.code,
        b.type,
        b.full_text
    FROM `{Sys.getenv('WORKSPACE_CDR')}`.cb_criteria_ancestor ca
    JOIN (
        SELECT DISTINCT
          c.concept_id,
          c.type,
          c.code,
          c.name,
          c.full_text
        FROM `{Sys.getenv('WORKSPACE_CDR')}`.cb_criteria c
        JOIN (
            SELECT id
            FROM `{Sys.getenv('WORKSPACE_CDR')}`.cb_criteria
            WHERE concept_id IN ({paste(concepts, collapse = ', ')})
              AND full_text LIKE '%_rank1]%'
        ) anchor
        ON c.path LIKE CONCAT('%.', CAST(anchor.id AS STRING), '.%')
            OR c.path LIKE CONCAT('%.', CAST(anchor.id AS STRING))
            OR c.path LIKE CONCAT(CAST(anchor.id AS STRING), '.%')
            OR c.path = CAST(anchor.id AS STRING)
        WHERE c.is_standard = 1 AND c.is_selectable = 1
    ) b
    ON ca.ancestor_id = b.concept_id",
    .sep = ""
  )

  bigrquery::bq_table_download(
    bigrquery::bq_project_query(
      Sys.getenv("GOOGLE_PROJECT"),
      concepts_query
    )
  ) |>
    data.table::as.data.table()
}

# simplifies saving later
psave <- function(obj, file) {
  write_parquet(
    x = obj,
    sink = file,
    compression = "zstd",
    compression_level = 12
  )
}

# remove phecode sex-subject sex discordant observations
restrictPhecodesBySex <- function(phecodes, id_sex, phe_sex) {
  # merge in subject sex
  data <- merge.data.table(phecodes, id_sex, by = "person_id", all.x = TRUE)
  # merge in phecode sex
  phe_sex[,
    phe_sex := fcase(
      male_only,
      "M",
      female_only,
      "F",
      default = "B"
    )
  ]
  data <- merge.data.table(
    data,
    phe_sex[, .(phecode, phe_sex)],
    all.x = TRUE,
    by = "phecode"
  )
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
  id_var = "person_id",
  chunk_size = 10000,
  id_sex = NULL,
  phe_sex = NULL,
  progress = TRUE,
  phe_var = "phecode",
  n_var = "N",
  sex_var = "sex"
) {
  # generate ID chunks
  ids <- unique(phe_dat[[id_var]])
  id_chunks <- base::split(ids, ceiling(seq_along(ids) / chunk_size))

  # iterate over the chunks
  res_list <- lapply(
    cli_progress_along(seq_along(id_chunks)),
    \(i) {
      data_sub <- phe_dat[phe_dat[[id_var]] %in% id_chunks[[i]], ]
      dcast(
        data_sub,
        formula = get(id_var) ~ get(phe_var),
        value.var = n_var
      )
    }
  )
  cli_progress_done()

  res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  res[is.na(res)] <- 0

  if (!is.null(id_sex) && !is.null(phe_sex)) {
    male_codes <- phe_sex[male_only == TRUE, phecode]
    male_codes <- male_codes[male_codes %in% names(res)]
    female_codes <- phe_sex[female_only == TRUE, phecode]
    female_codes <- female_codes[female_codes %in% names(res)]

    male_ids <- id_sex[tolower(substr(id_sex[[sex_var]], 1, 1)) == "m", ][[
      id_var
    ]]
    female_ids <- id_sex[tolower(substr(id_sex[[sex_var]], 1, 1)) == "f", ][[
      id_var
    ]]

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
          cases = sum(pim[[x]] == 1, na.rm = TRUE),
          controls = sum(pim[[x]] == 0, na.rm = TRUE),
          missing = sum(is.na(pim[[x]]), na.rm = TRUE)
        )
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )
  result[, total := cases + controls + missing][]
}

# calculate cci
cci_score <- function(icd_data, id_var = "person_id") {
  # check
  if (!data.table::is.data.table(icd_data)) {
    message("coercing icd_data to data.table object")
    icd_data <- data.table::as.data.table(icd_data)
  }

  # subset
  icd10_data <- icd_data[
    tolower(icd_data[["vocabulary_id"]]) %chin% c("icd10", "icd10cm"),
  ]
  icd9_data <- icd_data[
    tolower(icd_data[["vocabulary_id"]]) %chin% c("icd9", "icd9cm"),
  ]

  # get comorbidities
  results_icd10 <- comorbidity::comorbidity(
    icd10_data,
    id = id_var,
    code = "code",
    map = "charlson_icd10_quan",
    assign0 = TRUE
  ) |>
    data.table::as.data.table()
  results_icd9 <- comorbidity::comorbidity(
    icd9_data,
    id = id_var,
    code = "code",
    map = "charlson_icd9_quan",
    assign0 = TRUE
  ) |>
    data.table::as.data.table()

  # deduplicate
  comorbid_names <- names(results_icd10)[names(results_icd10) != id_var]
  names(results_icd10) <- c(id_var, paste0(comorbid_names, "10"))
  names(results_icd9) <- c(id_var, paste0(comorbid_names, "9"))

  res_merge <- data.table::merge.data.table(
    results_icd10,
    results_icd9,
    by = id_var,
    all = TRUE
  )
  res_merge[is.na(res_merge)] <- 0
  for (cn in comorbid_names) {
    res_merge[[cn]] <- as.numeric(
      (res_merge[[paste0(cn, "10")]] + res_merge[[paste0(cn, "9")]]) > 0
    )
  }

  select_vars <- c(id_var, comorbid_names)
  res_merge <- res_merge[, select_vars, with = FALSE]
  class(res_merge) <- c("comorbidity", "data.frame")
  attr(res_merge, "map") <- "charlson_icd10_quan"

  cci_score <- comorbidity::score(res_merge, weights = "quan", assign0 = TRUE)

  res_merge$cci_score <- cci_score

  data.table::as.data.table(res_merge)
}

categorize_quartiles <- function(dt, var) {
  # Compute quartiles
  q1 <- quantile(dt[[var]], 0.25, na.rm = TRUE)
  q3 <- quantile(dt[[var]], 0.75, na.rm = TRUE)

  # Create category column
  dt[,
    paste0(var, "_cat") := fcase(
      get(var) <= q1,
      "low",
      get(var) > q3,
      "high",
      get(var) > q1 & get(var) <= q3,
      "middle",
      default = "Unknown"
    )
  ]
}


# create factor levels
low_mid_high <- function(x) {
  factor(
    fcase(
      x == "low",
      "Low 25%",
      x == "middle",
      "Middle 50%",
      x == "high",
      "High 25%",
      x == "Unknown",
      "Unknown"
    ),
    levels = c("Low 25%", "Middle 50%", "High 25%", "Unknown")
  )
}

# table formatting function ----------------------------------------------------
tab_sum <- function(x, by_var = NULL, mod_span_head = NULL) {
  gtsummary::tbl_summary(
    x,
    by = by_var,
    type = list(
      dplyr::where(is.numeric) ~ "continuous",
      gtsummary::all_dichotomous() ~ "dichotomous"
    ),
    statistic = list(
      gtsummary::all_continuous() ~ "{mean} ({sd})",
      gtsummary::all_categorical() ~ "{p}% ({n})"
    ),
    digits = gtsummary::all_continuous() ~ 2,
    missing_text = "Missing"
  ) |>
    (\(y) {
      if (!is.null(by_var)) {
        y |>
          gtsummary::add_overall() |>
          gtsummary::add_p(
            pvalue_fun = gtsummary::label_style_pvalue(digits = 2),
            test = list(
              gtsummary::all_continuous() ~ "kruskal.test",
              gtsummary::all_categorical() ~ "chisq.test"
            ),
            test.args = list(
              gtsummary::all_categorical() ~
                list(simulate.p.value = TRUE, B = 50000)
            )
          )
      } else {
        y
      }
    })() |>
    (\(z) {
      if (!is.null(mod_span_head)) {
        z |> gtsummary::modify_spanning_header(mod_span_head)
      } else {
        z
      }
    })() |>
    gtsummary::bold_labels()
}

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
    dplyr::select(
      time,
      rmst_est = Est.,
      rmst_lo = `lower .95`,
      rmst_hi = `upper .95`,
      rmst_p = p
    ) |>
    tidyr::pivot_longer(
      cols = -time,
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(metric = paste0(metric, "_", time)) |>
    dplyr::select(-time) |>
    tidyr::pivot_wider(
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
        survival::coxph(formula, data = data, x = TRUE),
        warning = function(w) {
          if (grepl("loglik converg", w$message, ignore.case = TRUE)) {
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
          coxphf::coxphf(
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
      survival::cox.zph(cox_model) # Run proportional hazards test
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
