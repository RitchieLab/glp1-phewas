# libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
    library(arrow)
    library(cli)
    library(glue)
    library(data.table)
    library(tidyverse)
    library(patchwork)
    library(ggtext)
})

options(datatable.print.class = TRUE)
options(datatable.print.trunc.cols = TRUE)

project_directory <- "/Users/maxsal/Library/CloudStorage/Box-Box/projects/GLP-1 RA/PheWAS/"
version <- "20251103_r2r"
treatment <- "glp1"
comparator <- "sglt2"

phecodeX_labels <- fread(
    "https://github.com/PheWAS/PhecodeX/raw/main/phecodeX_R_labels.csv",
    showProgress = FALSE
)

source(
    "https://gitlab.com/maxsal/plot_phewasx/-/raw/main/plot_phewasx.R?ref_type=heads"
)
source(
    "/Users/maxsal/Library/CloudStorage/Box-Box/projects/GLP-1\ RA/PheWAS/R/fn/plot_rmst_phewasx.R"
)

# functions --------------------------------------------------------------------
quick_plot <- function(
    dt,
    no_ob_osa = TRUE,
    treatment = "GLP-1 RA",
    comparator = "SGLT2i"
) {
    too_small <- dt[is.na(beta), .N]
    cap <- if (too_small > 0) {
        glue(
            "Note: {too_small} phenotypes excluded because of too few occurrences (<50)"
        )
    } else {
        NULL
    }

    dt <- dt[!grepl("EM_202", phecode, ignore.case = TRUE)]
    if (no_ob_osa) {
        dt <- dt[
            !grepl("EM_236|NS_333|EM_204.5", phecode, ignore.case = TRUE)
        ]
    }

    itt_plot <- dt[type == "ITT" & ph_test_p > 0.05, ] |>
        plot_phewasx(phe_var = "phecode", color_x_labels = FALSE) +
        labs(
            title = paste0(
                "A. ",
                treatment,
                " vs. ",
                comparator,
                ", intention-to-treat"
            ),
            subtitle = paste0(
                "n=",
                format(
                    dt[type == "ITT" & ph_test_p > 0.05, .N],
                    big.mark = ","
                ),
                " phenotypes with sample sizes ranging from ",
                format(
                    min(dt[type == "ITT" & ph_test_p > 0.05, n_total]),
                    big.mark = ","
                ),
                " to ",
                format(
                    max(dt[type == "ITT" & ph_test_p > 0.05, n_total]),
                    big.mark = ","
                )
            ),
            caption = cap
        )

    pp_plot <- dt[type == "PP" & ph_test_p > 0.05, ] |>
        plot_phewasx(phe_var = "phecode", color_x_labels = FALSE) +
        labs(
            title = paste0(
                "B. ",
                treatment,
                " vs. ",
                comparator,
                ", per-protocol"
            ),
            subtitle = paste0(
                "n=",
                format(dt[type == "PP" & ph_test_p > 0.05, .N], big.mark = ","),
                " phenotypes with sample sizes ranging from ",
                format(
                    min(dt[type == "PP" & ph_test_p > 0.05, n_total]),
                    big.mark = ","
                ),
                " to ",
                format(
                    max(dt[type == "PP" & ph_test_p > 0.05, n_total]),
                    big.mark = ","
                )
            ),
            caption = cap
        )

    patched <- wrap_plots(
        itt_plot,
        pp_plot +
            theme(
                axis.title.y = element_blank()
            ),
        nrow = 1,
        ncol = 2
    )

    return(
        list(
            itt = itt_plot,
            pp = pp_plot,
            patched = patched
        )
    )
}

quick_rmst_plot <- function(
    dt,
    no_ob_osa = TRUE,
    treatment = "GLP-1 RA",
    comparator = "SGLT2i"
) {
    too_small <- dt[is.na(beta), .N]
    cap <- if (too_small > 0) {
        glue(
            "Note: {too_small} phenotypes excluded because of too few occurrences (<50)"
        )
    } else {
        NULL
    }

    dt <- dt[!grepl("EM_202", phecode, ignore.case = TRUE)]
    if (no_ob_osa) {
        dt <- dt[
            !grepl("EM_236|NS_333|EM_204.5", phecode, ignore.case = TRUE)
        ]
    }

    itt_plot <- dt[type == "ITT" & ph_test_p > 0.05, ] |>
        plot_rmst_phewasx(
            phe_var = "phecode",
            rmst_diff_type = "lines",
            annotate = TRUE,
            annotate_label = treatment, color_x_labels = FALSE
        ) +
        labs(
            title = paste0(
                "A. ",
                treatment,
                " vs. ",
                comparator,
                ", intention-to-treat"
            ),
            subtitle = paste0(
                "n=",
                format(
                    dt[type == "ITT" & ph_test_p > 0.05, .N],
                    big.mark = ","
                ),
                " phenotypes with sample sizes ranging from ",
                format(
                    min(dt[type == "ITT" & ph_test_p > 0.05, n_total]),
                    big.mark = ","
                ),
                " to ",
                format(
                    max(dt[type == "ITT" & ph_test_p > 0.05, n_total]),
                    big.mark = ","
                )
            ),
            caption = cap
        )

    pp_plot <- dt[type == "PP" & ph_test_p > 0.05, ] |>
        plot_rmst_phewasx(
            phe_var = "phecode",
            rmst_diff_type = "lines",
            annotate = TRUE,
            annotate_label = treatment,
            color_dot_symbol = NULL, color_x_labels = FALSE
        ) +
        labs(
            title = paste0(
                "B. ",
                treatment,
                " vs. ",
                comparator,
                ", per-protocol"
            ),
            subtitle = paste0(
                "n=",
                format(dt[type == "PP" & ph_test_p > 0.05, .N], big.mark = ","),
                " phenotypes with sample sizes ranging from ",
                format(
                    min(dt[type == "PP" & ph_test_p > 0.05, n_total]),
                    big.mark = ","
                ),
                " to ",
                format(
                    max(dt[type == "PP" & ph_test_p > 0.05, n_total]),
                    big.mark = ","
                )
            ),
            caption = cap
        )

    patched <- wrap_plots(
        itt_plot,
        pp_plot +
            theme(
                axis.title.y = element_blank()
            ),
        nrow = 1,
        ncol = 2
    )

    return(
        list(
            itt = itt_plot,
            pp = pp_plot,
            patched = patched
        )
    )
}


hr_rmst_plot <- function(
    data,
    top_n = 20,
    title = NULL,
    order_var = "or",
    rev = FALSE,
    rmst_diff = 30,
    rmst_diff_type = "rect",
    rect_col = "gray90"
) {
    # initialize
    cols <- c(
        "Significant" = "#D55E00",
        "Suggestive" = "#0072B2",
        "Not significant" = "#CC79A7"
    )
    shps <- c(
        "Significant" = 15,
        "Suggestive" = 16,
        "Not significant" = 17
    )

    flip <- 1
    if (rev) {
        flip <- -1
    }

    order_text <- ""
    if (order_var %in% c("or", "hr")) {
        order_text <- "  - Estimates ordered by hazard ratio"
    } else if (order_var %in% c("p", "p_adj")) {
        order_text <- "  - Estimates ordered by p-value"
    }
    # hazard ratio plot
    p1 <- data[order(p), ][1:top_n, ] |>
        ggplot(aes(x = reorder(print, flip * abs(get(order_var))), y = or)) +
        geom_hline(yintercept = 1, linetype = 1, color = "black") +
        geom_pointrange(
            aes(ymin = or_lo, ymax = or_hi, shape = sig, color = sig),
            size = 0.5
        ) +
        scale_color_manual(values = cols) +
        scale_shape_manual(values = shps) +
        labs(
            title = title,
            subtitle = paste0(
                "Top ",
                top_n,
                " associations with smallest p-values"
            ),
            x = "",
            y = "Hazard ratio (95% CI)",
            caption = paste0(
                "Note(s): \n  - ",
                fcase(
                    rmst_diff_type == "lines",
                    paste0("Dashed lines at ±", rmst_diff, " days"),
                    grepl("rect", rmst_diff_type, ignore.case = TRUE),
                    paste0("Shaded area at ±", rmst_diff, " days"),
                    default = ""
                ),
                "\n",
                "  - Significance determined by hazard ratio p-value",
                "\n",
                order_text
            )
        ) +
        coord_flip() +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0, face = "bold"),
            plot.caption = element_text(hjust = 0),
            legend.position = "top",
            legend.title = element_blank()
        )

    # RMST difference plot
    p2 <- data[order(p), ][1:top_n, ] |>
        ggplot(aes(
            x = reorder(print, flip * abs(get(order_var))),
            y = rmst_est_3
        )) +
        (\(x)
            if (rmst_diff_type == "lines") {
                geom_hline(
                    yintercept = c(-rmst_diff, rmst_diff),
                    linetype = 2,
                    color = "black"
                )
            })() +
        (\(x)
            if (grepl("rect", rmst_diff_type, ignore.case = TRUE)) {
                geom_rect(
                    aes(
                        xmin = -Inf,
                        xmax = Inf,
                        ymin = -rmst_diff,
                        ymax = rmst_diff
                    ),
                    fill = rect_col,
                    alpha = 0.25
                )
            })() +
        geom_hline(yintercept = 0, linetype = 1, color = "black") +
        geom_pointrange(
            aes(ymin = rmst_lo_3, ymax = rmst_hi_3, shape = sig, color = sig),
            size = 0.5
        ) +
        scale_color_manual(values = cols) +
        scale_shape_manual(values = shps) +
        labs(
            x = "",
            y = "3-year RMST difference (95% CI), days"
        ) +
        coord_flip() +
        theme_minimal() +
        theme(
            axis.text.y = element_blank(),
            legend.position = "top",
            legend.title = element_blank()
        )

    # combine plots
    wrap_plots(
        p1,
        p2,
        design = "
        AAABBB
        ",
        widths = c(2, 1),
        guides = "collect"
    ) &
        plot_annotation(
            theme = theme(
                plot.title = element_text(hjust = 0, face = "bold"),
                legend.position = "top"
            )
        )
}


# load data --------------------------------------------------------------------
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

# plot -------------------------------------------------------------------------
## Manhattan plots
sglt2i_plots <- quick_plot(sglt2, treatment = "GLP-1 RA", comparator = "SGLT2i")
sglt2i_rmst_plots <- quick_rmst_plot(
    sglt2,
    treatment = "GLP-1 RA",
    comparator = "SGLT2i"
)

ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("glp1_sglt2_phewas_", version, "_plot.pdf")
    ),
    sglt2i_plots[["patched"]],
    width = 12,
    height = 6,
    device = cairo_pdf
)
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("glp1_sglt2_rmst_", version, "_plot.pdf")
    ),
    sglt2i_rmst_plots[["patched"]],
    width = 12,
    height = 6,
    device = cairo_pdf
)

dpp4i_plots <- quick_plot(dpp4, treatment = "GLP-1 RA", comparator = "DPP4i")
dpp4i_rmst_plots <- quick_rmst_plot(
    dpp4,
    treatment = "GLP-1 RA",
    comparator = "DPP4i"
)

ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("glp1_dpp4i_phewas_", version, "_plot.pdf")
    ),
    dpp4i_plots[["patched"]],
    width = 12,
    height = 6,
    device = cairo_pdf
)
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("glp1_dpp4i_rmst_", version, "_plot.pdf")
    ),
    dpp4i_rmst_plots[["patched"]],
    width = 12,
    height = 6,
    device = cairo_pdf
)

sema_sglt2_plots <- quick_plot(
    sema_sglt2,
    treatment = "Semaglutide",
    comparator = "SGLT2i"
)
sema_sglt2_rmst_plots <- quick_rmst_plot(
    sema_sglt2,
    treatment = "Semaglutide",
    comparator = "SGLT2i"
)
sema_dpp4_plots <- quick_plot(
    sema_dpp4,
    treatment = "Semaglutide",
    comparator = "DPP4i"
)
sema_dpp4_rmst_plots <- quick_rmst_plot(
    sema_dpp4,
    treatment = "Semaglutide",
    comparator = "DPP4i"
)

ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("sema_sglt2_phewas_", version, "_plot.pdf")
    ),
    sema_sglt2_plots[["patched"]],
    width = 12,
    height = 6,
    device = cairo_pdf
)
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("sema_sglt2_rmst_", version, "_plot.pdf")
    ),
    sema_sglt2_rmst_plots[["patched"]],
    width = 12,
    height = 6,
    device = cairo_pdf
)

ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("sema_dpp4i_phewas_", version, "_plot.pdf")
    ),
    sema_dpp4_plots[["patched"]],
    width = 12,
    height = 6,
    device = cairo_pdf
)
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("sema_dpp4i_rmst_", version, "_plot.pdf")
    ),
    sema_dpp4_rmst_plots[["patched"]],
    width = 12,
    height = 6,
    device = cairo_pdf
)

glp1_patched <- wrap_plots(
    sglt2i_plots[[1]] +
        labs(
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.text.x = element_blank()
        ),
    sglt2i_plots[[2]] +
        labs(
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.text.x = element_blank(),
            axis.title.y = element_blank()
        ),
    dpp4i_plots[[1]] +
        labs(
            title = "C. GLP-1 RA vs. DPP4i, intention-to-treat",
            subtitle = "",
            caption = ""
        ),
    dpp4i_plots[[2]] +
        labs(
            title = "D. GLP-1 RA vs. DPP4i, per-protocol",
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.title.y = element_blank()
        ),
    ncol = 2,
    nrow = 2
)
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("glp1_phewas_", version, "_plot.pdf")
    ),
    glp1_patched,
    width = 12,
    height = 8,
    device = cairo_pdf
)
glp1_rmst_patched <- wrap_plots(
    sglt2i_rmst_plots[[1]] +
        labs(
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.text.x = element_blank()
        ),
    sglt2i_rmst_plots[[2]] +
        labs(
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.text.x = element_blank(),
            axis.title.y = element_blank()
        ),
    dpp4i_rmst_plots[[1]] +
        labs(
            title = "C. GLP-1 RA vs. DPP4i, intention-to-treat",
            subtitle = "",
            caption = ""
        ),
    dpp4i_rmst_plots[[2]] +
        labs(
            title = "D. GLP-1 RA vs. DPP4i, per-protocol",
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.title.y = element_blank()
        ),
    ncol = 2,
    nrow = 2
)
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("glp1_rmst_", version, "_plot.pdf")
    ),
    glp1_rmst_patched,
    width = 12,
    height = 8,
    device = cairo_pdf
)

sema_patched <- wrap_plots(
    sema_sglt2_plots[[1]] +
        labs(
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.text.x = element_blank()
        ),
    sema_sglt2_plots[[2]] +
        labs(
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.text.x = element_blank(),
            axis.title.y = element_blank()
        ),
    sema_dpp4_plots[[1]] +
        labs(
            title = "C. Semaglutide vs. DPP4i, intention-to-treat",
            subtitle = "",
            caption = ""
        ),
    sema_dpp4_plots[[2]] +
        labs(
            title = "D. Semaglutide vs. DPP4i, per-protocol",
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.title.y = element_blank()
        ),
    ncol = 2,
    nrow = 2
)
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("sema_phewas_", version, "_plot.pdf")
    ),
    sema_patched,
    width = 12,
    height = 8,
    device = cairo_pdf
)

sema_rmst_patched <- wrap_plots(
    sema_sglt2_rmst_plots[[1]] +
        labs(
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.text.x = element_blank()
        ),
    sema_sglt2_rmst_plots[[2]] +
        labs(
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.text.x = element_blank(),
            axis.title.y = element_blank()
        ),
    sema_dpp4_rmst_plots[[1]] +
        labs(
            title = "C. Semaglutide vs. DPP4i, intention-to-treat",
            subtitle = "",
            caption = ""
        ),
    sema_dpp4_rmst_plots[[2]] +
        labs(
            title = "D. Semaglutide vs. DPP4i, per-protocol",
            subtitle = "",
            caption = ""
        ) +
        theme(
            axis.title.y = element_blank()
        ),
    ncol = 2,
    nrow = 2
)
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("sema_rmst_", version, "_plot.pdf")
    ),
    sema_rmst_patched,
    width = 12,
    height = 8,
    device = cairo_pdf
)

## HR and RMST side-by-side plots ----------------------------------------------
sglt2_hr_rmst_itt_plot <- merge(
    sglt2,
    phecodeX_labels[, .(phecode = phenotype, description)],
    by = "phecode",
    all.x = TRUE
)[
    !is.na(beta) & type == "ITT" & ph_test_p > 0.05,
][
    !grepl("EM_236|NS_333|EM_202|EM_204.5", phecode, ignore.case = TRUE),
] |>
    mutate(
        print = str_wrap(paste0("[", phecode, "] ", description), 30),
        sig = fcase(
            p < 0.05 / n(),
            "Significant",
            p < 0.05,
            "Suggestive",
            default = "Not significant"
        )
    ) |>
    hr_rmst_plot(
        title = "GLP-1 RA vs. SGLT2i, intention-to-treat",
        order_var = "or"
    )
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("glp1_sglt2_hr_rmst_", version, "_itt_plot.pdf")
    ),
    sglt2_hr_rmst_itt_plot,
    width = 10,
    height = 10,
    device = cairo_pdf
)

sglt2_hr_rmst_pp_plot <- merge(
    sglt2,
    phecodeX_labels[, .(phecode = phenotype, description)],
    by = "phecode",
    all.x = TRUE
)[
    !is.na(beta) & type == "PP" & ph_test_p > 0.05,
][
    !grepl("EM_236|NS_333|EM_202|EM_204.5", phecode, ignore.case = TRUE),
] |>
    mutate(
        print = str_wrap(paste0("[", phecode, "] ", description), 30),
        sig = fcase(
            p < 0.05 / n(),
            "Significant",
            p < 0.05,
            "Suggestive",
            default = "Not significant"
        )
    ) |>
    hr_rmst_plot(
        title = "GLP-1 RA vs. SGLT2i, per-protocol",
        order_var = "or"
    )
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("glp1_sglt2_hr_rmst_", version, "_pp_plot.pdf")
    ),
    sglt2_hr_rmst_pp_plot,
    width = 10,
    height = 10,
    device = cairo_pdf
)

dpp4_hr_rmst_itt_plot <- merge(
    dpp4,
    phecodeX_labels[, .(phecode = phenotype, description)],
    by = "phecode",
    all.x = TRUE
)[
    !is.na(beta) & type == "ITT" & ph_test_p > 0.05,
][
    !grepl("EM_236|NS_333|EM_202|EM_204.5", phecode, ignore.case = TRUE),
] |>
    mutate(
        print = str_wrap(paste0("[", phecode, "] ", description), 30),
        sig = fcase(
            p < 0.05 / n(),
            "Significant",
            p < 0.05,
            "Suggestive",
            default = "Not significant"
        )
    ) |>
    hr_rmst_plot(
        title = "GLP-1 RA vs. DPP4i, intention-to-treat",
        order_var = "or"
    )
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("glp1_dpp4_hr_rmst_", version, "_itt_plot.pdf")
    ),
    dpp4_hr_rmst_itt_plot,
    width = 10,
    height = 10,
    device = cairo_pdf
)

dpp4_hr_rmst_pp_plot <- merge(
    dpp4,
    phecodeX_labels[, .(phecode = phenotype, description)],
    by = "phecode",
    all.x = TRUE
)[
    !is.na(beta) & type == "PP" & ph_test_p > 0.05,
][
    !grepl("EM_236|NS_333|EM_202|EM_204.5", phecode, ignore.case = TRUE),
] |>
    mutate(
        print = str_wrap(paste0("[", phecode, "] ", description), 30),
        sig = fcase(
            p < 0.05 / n(),
            "Significant",
            p < 0.05,
            "Suggestive",
            default = "Not significant"
        )
    ) |>
    hr_rmst_plot(
        title = "GLP-1 RA vs. DPP4i, per-protocol",
        order_var = "or"
    )
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("glp1_dpp4_hr_rmst_", version, "_pp_plot.pdf")
    ),
    dpp4_hr_rmst_pp_plot,
    width = 10,
    height = 10,
    device = cairo_pdf
)

sema_sglt2_hr_rmst_itt_plot <- merge(
    sema_sglt2,
    phecodeX_labels[, .(phecode = phenotype, description)],
    by = "phecode",
    all.x = TRUE
)[
    !is.na(beta) & type == "ITT" & ph_test_p > 0.05,
][
    !grepl("EM_236|NS_333|EM_202|EM_204.5", phecode, ignore.case = TRUE),
] |>
    mutate(
        print = str_wrap(paste0("[", phecode, "] ", description), 30),
        sig = fcase(
            p < 0.05 / n(),
            "Significant",
            p < 0.05,
            "Suggestive",
            default = "Not significant"
        )
    ) |>
    hr_rmst_plot(
        title = "Semaglutide vs. SGLT2i, intention-to-treat",
        order_var = "or"
    )
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("sema_sglt2_hr_rmst_", version, "_itt_plot.pdf")
    ),
    sema_sglt2_hr_rmst_itt_plot,
    width = 10,
    height = 10,
    device = cairo_pdf
)

sema_sglt2_hr_rmst_pp_plot <- merge(
    sema_sglt2,
    phecodeX_labels[, .(phecode = phenotype, description)],
    by = "phecode",
    all.x = TRUE
)[
    !is.na(beta) & type == "PP" & ph_test_p > 0.05,
][
    !grepl("EM_236|NS_333|EM_202|EM_204.5", phecode, ignore.case = TRUE),
] |>
    mutate(
        print = str_wrap(paste0("[", phecode, "] ", description), 30),
        sig = fcase(
            p < 0.05 / n(),
            "Significant",
            p < 0.05,
            "Suggestive",
            default = "Not significant"
        )
    ) |>
    hr_rmst_plot(
        title = "Semaglutide vs. SGLT2i, per-protocol",
        order_var = "or"
    )
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("sema_sglt2_hr_rmst_", version, "_pp_plot.pdf")
    ),
    sema_sglt2_hr_rmst_pp_plot,
    width = 10,
    height = 10,
    device = cairo_pdf
)

sema_dpp4_hr_rmst_itt_plot <- merge(
    sema_dpp4,
    phecodeX_labels[, .(phecode = phenotype, description)],
    by = "phecode",
    all.x = TRUE
)[
    !is.na(beta) & type == "ITT" & ph_test_p > 0.05,
][
    !grepl("EM_236|NS_333|EM_202|EM_204.5", phecode, ignore.case = TRUE),
] |>
    mutate(
        print = str_wrap(paste0("[", phecode, "] ", description), 30),
        sig = fcase(
            p < 0.05 / n(),
            "Significant",
            p < 0.05,
            "Suggestive",
            default = "Not significant"
        )
    ) |>
    hr_rmst_plot(
        title = "Semaglutide vs. DPP4i, intention-to-treat",
        order_var = "or"
    )
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("sema_dpp4_hr_rmst_", version, "_itt_plot.pdf")
    ),
    sema_dpp4_hr_rmst_itt_plot,
    width = 10,
    height = 10,
    device = cairo_pdf
)

sema_dpp4_hr_rmst_pp_plot <- merge(
    sema_dpp4,
    phecodeX_labels[, .(phecode = phenotype, description)],
    by = "phecode",
    all.x = TRUE
)[
    !is.na(beta) & type == "PP" & ph_test_p > 0.05,
][
    !grepl("EM_236|NS_333|EM_202|EM_204.5", phecode, ignore.case = TRUE),
] |>
    mutate(
        print = str_wrap(paste0("[", phecode, "] ", description), 30),
        sig = fcase(
            p < 0.05 / n(),
            "Significant",
            p < 0.05,
            "Suggestive",
            default = "Not significant"
        )
    ) |>
    hr_rmst_plot(
        title = "Semaglutide vs. DPP4i, per-protocol",
        order_var = "or"
    )
ggsave(
    file.path(
        project_directory,
        "figures",
        paste0("sema_dpp4_hr_rmst_", version, "_pp_plot.pdf")
    ),
    sema_dpp4_hr_rmst_pp_plot,
    width = 10,
    height = 10,
    device = cairo_pdf
)
