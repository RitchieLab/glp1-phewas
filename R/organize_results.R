library(data.table)
library(glue)
library(cli)

options(datatable.print.class = TRUE)
options(datatable.print.trunc.cols = TRUE)

project_directory <- "/Users/maxsal/Library/CloudStorage/Box-Box/projects/GLP-1 RA/PheWAS/"
version <- "20251103_r2r"

phecodeX_labels <- fread(
  "https://github.com/PheWAS/PhecodeX/raw/main/phecodeX_R_labels.csv",
  showProgress = FALSE
)

sglt2 <- fread(file.path(
  project_directory,
  "results",
  paste0("glp_1_ra_sglt2i_phewas_", version, ".csv")
))[, `:=`(
  log10p = log10(p),
  n_total = n_cases + n_controls
)][]
dpp4 <- fread(file.path(
  project_directory,
  "results",
  paste0("glp_1_ra_dpp4i_phewas_", version, ".csv")
))[, `:=`(
  log10p = log10(p),
  n_total = n_cases + n_controls
)][]

sema_sglt2 <- fread(file.path(
  project_directory,
  "results",
  paste0("semaglutide_sglt2i_phewas_", version, ".csv")
))[, `:=`(
  log10p = log10(p),
  n_total = n_cases + n_controls
)][]
sema_dpp4 <- fread(file.path(
  project_directory,
  "results",
  paste0("semaglutide_dpp4i_phewas_", version, ".csv")
))[, `:=`(
  log10p = log10(p),
  n_total = n_cases + n_controls
)][]

organize_results <- function(
  x,
  phe_labs = phecodeX_labels,
  rm_these = c("EM_202", "EM_236", "NS_333", "EM_204.5"),
  report_summary = TRUE,
  detailed_summary = FALSE
) {
  x <- copy(x)

  # Progressive filtering with counts
  n_initial <- x[, uniqueN(phecode)]

  x <- x[
    phe_labs[, .(phecode = phenotype, description, group)],
    on = "phecode",
    nomatch = 0
  ]
  setnames(x, old = c("or", "or_lo", "or_hi"), new = c("hr", "hr_lo", "hr_hi"))

  # Tag insufficient cases BEFORE filtering
  x[is.na(beta), note := "Not tested: insufficient cases after matching"]
  n_insufficient <- sum(is.na(x$beta))

  # Identify exclusions
  excluded_a_priori <- grepl(
    paste0(rm_these, collapse = "|"),
    x[, phecode],
    ignore.case = TRUE
  )
  excluded_ph <- x$ph_test_p < 0.05 | is.na(x$ph_test_p)
  excluded_insufficient <- x[, is.na(beta)]

  # Count exclusions by type BEFORE filtering
  n_excluded_a_priori <- sum(excluded_a_priori & !excluded_insufficient) / 2
  n_excluded_ph_itt <- sum(
    excluded_ph & x$type == "ITT" & !excluded_insufficient,
    na.rm = TRUE
  )
  n_excluded_ph_pp <- sum(
    excluded_ph & x$type == "PP" & !excluded_insufficient,
    na.rm = TRUE
  )
  n_excluded_ph_na <- sum(is.na(x$ph_test_p) & !excluded_insufficient)
  n_excluded_ph_fail <- sum(
    x$ph_test_p < 0.05 & !excluded_insufficient,
    na.rm = TRUE
  )

  # Set exclusion notes (check for existing notes to avoid overwriting)
  x[
    is.na(ph_test_p) & is.na(note),
    note := "Excluded: PH test not performed (NA)"
  ]
  x[
    excluded_ph & excluded_a_priori & is.na(note),
    note := "Excluded: a priori AND PH violation"
  ]
  x[
    excluded_ph & !excluded_a_priori & is.na(note),
    note := "Excluded: PH assumption violation (p < 0.05)"
  ]
  x[
    excluded_a_priori & !excluded_ph & is.na(note),
    note := "Excluded: a priori (diabetes/obesity/sleep apnea)"
  ]

  # Split into exclusions and main analysis
  x_rm <- x[excluded_a_priori | excluded_ph | excluded_insufficient, ]
  x <- x[!excluded_a_priori & !excluded_ph & !excluded_insufficient, ]

  # Calculate total sample size per phenotype (for main analysis set only)
  x[, n_total := n_cases + n_controls]

  # Process ITT
  itt <- x[type == "ITT", ]
  n_itt <- nrow(itt)
  itt[, pwide_sig_threshold := 0.05 / .N]
  itt[,
    significant := fcase(
      p < pwide_sig_threshold,
      "Phenome-wide significant",
      p < 0.05,
      "Suggestive",
      default = "Not significant"
    )
  ]

  # Create formatted HR string for display
  itt[, hr_print := sprintf("%.2f (%.2f, %.2f)", hr, hr_lo, hr_hi)]

  # Calculate sample size stats for ITT
  sample_size_stats_itt <- list(
    min = min(itt$n_total, na.rm = TRUE),
    max = max(itt$n_total, na.rm = TRUE),
    median = median(itt$n_total, na.rm = TRUE),
    min_phecode = itt[
      n_total == min(n_total, na.rm = TRUE),
      .(phecode = phecode[1], description = description[1])
    ],
    max_phecode = itt[
      n_total == max(n_total, na.rm = TRUE),
      .(phecode = phecode[1], description = description[1])
    ]
  )

  # Calculate direction stats for ITT
  itt_sig_sug <- itt[
    significant %in% c("Phenome-wide significant", "Suggestive")
  ]
  n_itt_risk <- sum(itt_sig_sug$beta > 0, na.rm = TRUE)
  n_itt_prot <- sum(itt_sig_sug$beta < 0, na.rm = TRUE)

  # Group summaries for ITT
  itt_group_summary <- itt[
    significant %in% c("Phenome-wide significant", "Suggestive"),
    .(
      n_total = .N,
      n_risk = sum(beta > 0, na.rm = TRUE),
      n_protective = sum(beta < 0, na.rm = TRUE)
    ),
    by = group
  ][order(-n_total)]

  # Extract top results for ITT
  itt_phenome_wide <- itt[significant == "Phenome-wide significant"][order(p)]
  itt_top_suggestive <- itt[significant == "Suggestive"][order(p)][1:min(5, .N)]

  # Create display tables for ITT
  itt_phenome_wide_table <- itt_phenome_wide[, .(
    phecode,
    description,
    group,
    hr_print,
    p = sprintf("%.2e", p),
    n_total
  )]

  itt_top_suggestive_table <- itt_top_suggestive[, .(
    phecode,
    description,
    group,
    hr_print,
    p = sprintf("%.2e", p),
    n_total
  )]

  # Process PP
  pp <- x[type == "PP", ]
  n_pp <- nrow(pp)
  pp[, pwide_sig_threshold := 0.05 / .N]
  pp[,
    significant := fcase(
      p < pwide_sig_threshold,
      "Phenome-wide significant",
      p < 0.05,
      "Suggestive",
      default = "Not significant"
    )
  ]

  # Create formatted HR string for display
  pp[, hr_print := sprintf("%.2f (%.2f, %.2f)", hr, hr_lo, hr_hi)]

  # Calculate sample size stats for PP
  sample_size_stats_pp <- list(
    min = min(pp$n_total, na.rm = TRUE),
    max = max(pp$n_total, na.rm = TRUE),
    median = median(pp$n_total, na.rm = TRUE),
    min_phecode = pp[
      n_total == min(n_total, na.rm = TRUE),
      .(phecode = phecode[1], description = description[1])
    ],
    max_phecode = pp[
      n_total == max(n_total, na.rm = TRUE),
      .(phecode = phecode[1], description = description[1])
    ]
  )

  # Calculate direction stats for PP
  pp_sig_sug <- pp[significant %in% c("Phenome-wide significant", "Suggestive")]
  n_pp_risk <- sum(pp_sig_sug$beta > 0, na.rm = TRUE)
  n_pp_prot <- sum(pp_sig_sug$beta < 0, na.rm = TRUE)

  # Group summaries for PP
  pp_group_summary <- pp[
    significant %in% c("Phenome-wide significant", "Suggestive"),
    .(
      n_total = .N,
      n_risk = sum(beta > 0, na.rm = TRUE),
      n_protective = sum(beta < 0, na.rm = TRUE)
    ),
    by = group
  ][order(-n_total)]

  # Extract top results for PP
  pp_phenome_wide <- pp[significant == "Phenome-wide significant"][order(p)]
  pp_top_suggestive <- pp[significant == "Suggestive"][order(p)][1:min(5, .N)]

  # Create display tables for PP
  pp_phenome_wide_table <- pp_phenome_wide[, .(
    phecode,
    description,
    group,
    hr_print,
    p = sprintf("%.2e", p),
    n_total
  )]

  pp_top_suggestive_table <- pp_top_suggestive[, .(
    phecode,
    description,
    group,
    hr_print,
    p = sprintf("%.2e", p),
    n_total
  )]

  # Handle other types (shouldn't be any, but just in case)
  other <- x[!type %in% c("ITT", "PP"), ]

  # Combine all results
  out <- rbindlist(list(itt, pp, x_rm, other), use.names = TRUE, fill = TRUE)

  # Order columns
  order_cols <- c(
    "phecode",
    "description",
    "group",
    "type",
    "beta",
    "hr",
    "hr_lo",
    "hr_hi",
    "p",
    "pwide_sig_threshold",
    "significant",
    "note"
  )
  order_cols <- c(order_cols, setdiff(names(out), order_cols))

  # Report summary
  if (report_summary == TRUE) {
    cli_h1("organize_results() summary")

    cli_h2("Filtering pipeline:")
    cli_alert_info("Starting phenotypes: {.strong {n_initial}}")
    cli_alert_info("Excluded a priori: {.strong {n_excluded_a_priori}}")
    cli_alert_info(
      "Excluded (PH test failure p<0.05): {.strong {n_excluded_ph_fail}}"
    )
    cli_alert_info("Excluded (PH test NA): {.strong {n_excluded_ph_na}}")
    cli_alert_info("Excluded (insufficient cases): {.strong {n_insufficient}}")

    cli_h2("Final analysis set:")
    cli_alert_info(
      "ITT: {.strong {n_itt}} phenotypes, Bonferroni threshold = {.strong {formatC(0.05/n_itt, format='e', digits=2)}}"
    )
    cli_alert_info(
      "PP:  {.strong {n_pp}} phenotypes, Bonferroni threshold = {.strong {formatC(0.05/n_pp, format='e', digits=2)}}"
    )

    cli_h2("Sample size range:")
    cli_h3("ITT:")
    cli_alert_info(
      "Min: {.strong {sample_size_stats_itt$min}} ({.emph {sample_size_stats_itt$min_phecode$phecode}}: {sample_size_stats_itt$min_phecode$description})"
    )
    cli_alert_info(
      "Max: {.strong {sample_size_stats_itt$max}} ({.emph {sample_size_stats_itt$max_phecode$phecode}}: {sample_size_stats_itt$max_phecode$description})"
    )
    cli_alert_info("Median: {.strong {sample_size_stats_itt$median}}")

    cli_h3("PP:")
    cli_alert_info(
      "Min: {.strong {sample_size_stats_pp$min}} ({.emph {sample_size_stats_pp$min_phecode$phecode}}: {sample_size_stats_pp$min_phecode$description})"
    )
    cli_alert_info(
      "Max: {.strong {sample_size_stats_pp$max}} ({.emph {sample_size_stats_pp$max_phecode$phecode}}: {sample_size_stats_pp$max_phecode$description})"
    )
    cli_alert_info("Median: {.strong {sample_size_stats_pp$median}}")

    cli_h2("Phenome-wide significant results:")
    n_itt_sig <- nrow(itt_phenome_wide)
    n_pp_sig <- nrow(pp_phenome_wide)
    cli_alert_info(
      "ITT: {.strong {n_itt_sig}} ({.emph {round(n_itt_sig/n_itt*100, 1)}%})"
    )
    if (n_itt_sig > 0) {
      cli_h3("ITT phenome-wide significant:")
      print(itt_phenome_wide_table)
    }

    cli_alert_info(
      "PP:  {.strong {n_pp_sig}} ({.emph {round(n_pp_sig/n_pp*100, 1)}%})"
    )
    if (n_pp_sig > 0) {
      cli_h3("PP phenome-wide significant:")
      print(pp_phenome_wide_table)
    }

    cli_h2("Suggestive results (p < 0.05):")
    n_itt_sug <- sum(itt$significant == "Suggestive", na.rm = TRUE)
    n_pp_sug <- sum(pp$significant == "Suggestive", na.rm = TRUE)
    cli_alert_info(
      "ITT: {.strong {n_itt_sug}} ({.emph {round(n_itt_sug/n_itt*100, 1)}%})"
    )
    if (nrow(itt_top_suggestive_table) > 0) {
      cli_h3("ITT top 5 suggestive:")
      print(itt_top_suggestive_table)
    }

    cli_alert_info(
      "PP:  {.strong {n_pp_sug}} ({.emph {round(n_pp_sug/n_pp*100, 1)}%})"
    )
    if (nrow(pp_top_suggestive_table) > 0) {
      cli_h3("PP top 5 suggestive:")
      print(pp_top_suggestive_table)
    }

    cli_h2("Direction of effects (significant + suggestive):")
    total_itt <- n_itt_risk + n_itt_prot
    total_pp <- n_pp_risk + n_pp_prot
    if (total_itt > 0) {
      cli_alert_info(
        "ITT: {.strong {n_itt_risk}} risk-increasing ({.emph {round(n_itt_risk/total_itt*100, 1)}%}), {.strong {n_itt_prot}} protective ({.emph {round(n_itt_prot/total_itt*100, 1)}%})"
      )
    }
    if (total_pp > 0) {
      cli_alert_info(
        "PP:  {.strong {n_pp_risk}} risk-increasing ({.emph {round(n_pp_risk/total_pp*100, 1)}%}), {.strong {n_pp_prot}} protective ({.emph {round(n_pp_prot/total_pp*100, 1)}%})"
      )
    }

    # Detailed summary if requested
    if (detailed_summary == TRUE) {
      cli_h2("Group breakdown (ITT - significant + suggestive):")
      if (nrow(itt_group_summary) > 0) {
        print(itt_group_summary)
      }

      cli_h2("Group breakdown (PP - significant + suggestive):")
      if (nrow(pp_group_summary) > 0) {
        print(pp_group_summary)
      }
    }
  }

  # Store comprehensive statistics as attributes
  summary_stats <- list(
    filtering = list(
      n_initial = n_initial,
      n_excluded_a_priori = n_excluded_a_priori,
      n_excluded_ph_fail = n_excluded_ph_fail,
      n_excluded_ph_na = n_excluded_ph_na,
      n_excluded_ph = c(ITT = n_excluded_ph_itt, PP = n_excluded_ph_pp),
      n_insufficient_cases = n_insufficient
    ),
    n_tested = c(ITT = n_itt, PP = n_pp),
    bonferroni_thresholds = c(ITT = 0.05 / n_itt, PP = 0.05 / n_pp),
    sample_sizes = list(
      ITT = sample_size_stats_itt,
      PP = sample_size_stats_pp
    ),
    results = list(
      n_phenome_wide_sig = c(ITT = n_itt_sig, PP = n_pp_sig),
      n_suggestive = c(ITT = n_itt_sug, PP = n_pp_sug),
      n_risk_increasing = c(ITT = n_itt_risk, PP = n_pp_risk),
      n_protective = c(ITT = n_itt_prot, PP = n_pp_prot)
    ),
    group_summaries = list(
      ITT = itt_group_summary,
      PP = pp_group_summary
    ),
    # NEW: Top results tables
    top_results = list(
      ITT = list(
        phenome_wide_significant = itt_phenome_wide_table,
        top_5_suggestive = itt_top_suggestive_table
      ),
      PP = list(
        phenome_wide_significant = pp_phenome_wide_table,
        top_5_suggestive = pp_top_suggestive_table
      )
    )
  )

  setattr(out, "summary_stats", summary_stats)

  # Sort: non-excluded by p-value, then excluded at bottom
  out[order(is.na(significant), p), ..order_cols]
}


# Helper function to view top results (alternative to using attributes)
view_top_results <- function(organized_results, approach = "ITT") {
  stats <- attr(organized_results, "summary_stats")

  cli_h1("Top results for {approach}")

  # Phenome-wide significant
  pwide <- stats$top_results[[approach]]$phenome_wide_significant
  if (nrow(pwide) > 0) {
    cli_h2("Phenome-wide significant (n={nrow(pwide)}):")
    print(pwide)
  } else {
    cli_alert_warning("No phenome-wide significant results")
  }

  # Top 5 suggestive
  sug <- stats$top_results[[approach]]$top_5_suggestive
  if (nrow(sug) > 0) {
    cli_h2("Top 5 suggestive:")
    print(sug)
  } else {
    cli_alert_warning("No suggestive results")
  }

  invisible(list(phenome_wide = pwide, suggestive = sug))
}


# Helper function to extract top results for manuscript writing
get_top_results_text <- function(organized_results, approach = "ITT") {
  stats <- attr(organized_results, "summary_stats")

  pwide <- stats$top_results[[approach]]$phenome_wide_significant
  sug <- stats$top_results[[approach]]$top_5_suggestive

  list(
    n_phenome_wide = nrow(pwide),
    n_suggestive = stats$results$n_suggestive[approach],
    phenome_wide_table = pwide,
    top_suggestive_table = sug,
    # Formatted text snippets for manuscript
    phenome_wide_list = if (nrow(pwide) > 0) {
      paste(
        sprintf("%s (%s)", pwide$description, pwide$hr_print),
        collapse = ", "
      )
    } else {
      "none"
    }
  )
}

sglt2_res <- organize_results(sglt2)
dpp4_res <- organize_results(dpp4)

sema_sglt2_res <- organize_results(sema_sglt2)
sema_dpp4_res <- organize_results(sema_dpp4)

fwrite(
  x = sglt2_res,
  file = "results/glp1_sglt2_phewas_20251103_r2r_organized.csv"
)
fwrite(
  x = dpp4_res,
  file = "results/glp1_dpp4_phewas_20251103_r2r_organized.csv"
)

fwrite(
  x = sema_sglt2_res,
  file = "results/semaglutde_sglt2_phewas_20251103_r2r_organized.csv"
)
fwrite(
  x = sema_dpp4_res,
  file = "results/semaglutide_dpp4_phewas_20251103_r2r_organized.csv"
)
