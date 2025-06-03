source(
    "https://gitlab.com/maxsal/plot_phewasx/-/raw/main/plot_phewasx.R?ref_type=heads"
)

ss_print <- function(x) {
    paste0(
        glue::glue_collapse(x[["print_phecode"]], ", ", last = " and "),
        " (n=",
        formatC(x[, unique(n_total)], big.mark = ",", format = "d"),
        ")"
    )
}

summarize_phewas_results <- function(
    dt,
    .type = "ITT",
    remove = "EM_236|NS_333.1|EM_202|EM_204.5"
) {
    # HR + RMST
    dt1 <- merge(
        dt,
        .clean_phecode_labels()[, .(phecode, description, group)],
        by = "phecode",
        all.x = TRUE
    )[
        !is.na(beta),
    ][
        type == .type,
    ][
        ph_test_p > 0.05,
    ][
        !grepl(remove, phecode, ignore.case = TRUE),
    ][, `:=`(
        print_hr = paste0(
            trimws(format(round(or, 2), nsmall = 2)),
            " (",
            trimws(format(round(or_lo, 2), nsmall = 2)),
            ", ",
            trimws(format(round(or_hi, 2), nsmall = 2)),
            ")"
        ),
        print_p = formatC(p, format = "e", digits = 2),
        print_rmst_p = formatC(rmst_p_3, format = "e", digits = 2),
        print_rmst = paste0(
            trimws(format(round(rmst_est_3, 1), nsmall = 1)),
            " (",
            trimws(format(round(rmst_lo_3, 1), nsmall = 1)),
            ", ",
            trimws(format(round(rmst_hi_3, 1), nsmall = 1)),
            ")"
        ),
        print_phecode = paste0("[", phecode, "] ", description),
        risk_hr = fcase(
            or > 1 & p < 0.05,
            "GLP-1 risk",
            or < 1 & p < 0.05,
            "GLP-1 protective",
            default = "no association"
        ),
        risk_rmst = fcase(
            rmst_est_3 > 0 & rmst_p_3 < 0.05,
            "GLP-1 protective",
            rmst_est_3 < 0 & rmst_p_3 < 0.05,
            "GLP-1 risk",
            default = "no association"
        ),
        significance = fcase(
            p < 0.05 / .N,
            "Significant",
            p < 0.05,
            "Suggestive",
            default = "No association"
        )
    )]

    n_out <- dt1[, .N]
    tmp_min <- ss_print(dt1[
        n_total == min(n_total),
        .(print_phecode, description, n_total)
    ])
    tmp_max <- ss_print(dt1[
        n_total == max(n_total),
        .(print_phecode, description, n_total)
    ])

    # HR
    hr_sig <- dt1[p < 0.05 / .N, ]
    hr_sug <- dt1[p < 0.05, ]
    hr_sig_n <- hr_sig[, .N]
    hr_sug_n <- hr_sug[, .N]
    hr_sig[, .(phecode, description, print_hr)]
    hr_sug_risk <- hr_sug[, .N, risk_hr][, prop := round((N * 100 / sum(N)), 1)]
    hr_sug_risk_group <- hr_sug[, .N, .(risk_hr, group)][,
        .SD[N == max(N)],
        by = risk_hr
    ]

    # RMST
    rmst_sig <- dt1[rmst_p_3 < 0.05 / .N, ]
    rmst_sug <- dt1[rmst_p_3 < 0.05, ]
    rmst_sig_n <- rmst_sig[, .N]
    rmst_sug_n <- rmst_sug[, .N]
    rmst_sig[, .(phecode, description, print_rmst)]
    rmst_sug_risk <- rmst_sug[, .N, risk_rmst][,
        prop := round((N * 100 / sum(N)), 1)
    ]
    rmst_sug_risk_group <- rmst_sug[, .N, .(risk_rmst, group)][,
        .SD[N == max(N)],
        by = risk_rmst
    ]

    list(
        n_out = formatC(n_out, big.mark = ",", format = "d"),
        tmp_min = tmp_min,
        tmp_max = tmp_max,
        hr_sig = hr_sig[order(p), .(phecode, description, print_hr, print_p)],
        hr_sug = hr_sug[
            order(p),
            .(
                phecode,
                description,
                print_hr,
                print_p,
                significance
            )
        ],
        hr_sug_risk = hr_sug_risk,
        hr_sug_risk_group = hr_sug_risk_group,
        hr_sig_n = paste0(
            "There were ",
            hr_sig_n,
            " outcomes with significant HR."
        ),
        hr_sug_n = paste0(
            "There were ",
            hr_sug_n,
            " outcomes with suggestive (+sig) HR."
        ),
        rmst_sig = rmst_sig[
            order(rmst_p_3),
            .(
                phecode,
                description,
                print_rmst,
                print_rmst_p
            )
        ],
        rmst_sug = rmst_sug[
            order(rmst_p_3),
            .(
                phecode,
                description,
                print_rmst,
                print_rmst_p,
                significance
            )
        ],
        rmst_sig_n = paste0(
            "There were ",
            rmst_sig_n,
            " outcomes with significant RMST."
        ),
        rmst_sug_n = paste0(
            "There were ",
            rmst_sug_n,
            " outcomes with suggestive (+sig) RMST."
        ),
        rmst_sug_risk = rmst_sug_risk,
        rmst_sug_risk_group = rmst_sug_risk_group
    )
}
