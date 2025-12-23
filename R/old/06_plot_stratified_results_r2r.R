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

phecodeX_labels <- fread(glue("{project_directory}/tables/phecodeX_labels.csv"))

sglt2 <- fread(glue(
  "{project_directory}/results/glp_1_ra_sglt2i_stratified_{version}.csv"
))
sglt2 <- sglt2[
  phecodeX_labels[, .(phecode = phenotype, description)],
  on = "phecode",
  nomatch = 0
]
dpp4 <- fread(glue(
  "{project_directory}/results/glp_1_ra_dpp4i_stratified_{version}.csv"
))
dpp4 <- dpp4[
  phecodeX_labels[, .(phecode = phenotype, description)],
  on = "phecode",
  nomatch = 0
]

sema_sglt2 <- fread(glue(
  "{project_directory}/results/semaglutide_sglt2i_stratified_{version}.csv"
))
sema_sglt2 <- sema_sglt2[
  phecodeX_labels[, .(phecode = phenotype, description)],
  on = "phecode",
  nomatch = 0
]
sema_dpp4 <- fread(glue(
  "{project_directory}/results/semaglutide_dpp4i_stratified_{version}.csv"
))
sema_dpp4 <- sema_dpp4[
  phecodeX_labels[, .(phecode = phenotype, description)],
  on = "phecode",
  nomatch = 0
]

sex_cols <- c(
  "Female" = "#CC79A7",
  "Male" = "#0072B2"
)

race_eth_cols <- c(
  "NHW" = "#56B4E9",
  "NHB" = "#009E73",
  "HISP" = "#D55E00"
)

quick_plot <- function(x, strata = c("Male", "Female"), dodge_width = 0.75) {
  if (any(grepl("male", strata, ignore.case = TRUE))) {
    cols <- sex_cols
  } else {
    cols <- race_eth_cols
  }
  x[stratum %in% strata, ] |>
    mutate(
      Approach = fifelse(type == "", NA, type),
      Strata = stratum,
      hr = fifelse(hr == 1, NA, hr)
    ) |>
    ggplot(aes(
      x = stringr::str_wrap(description, width = 25),
      y = hr,
      color = Strata,
      shape = Approach,
      linetype = Approach
    )) +
    geom_hline(yintercept = 1) +
    geom_pointrange(
      aes(ymin = hr_lo, ymax = hr_hi),
      position = position_dodge(dodge_width)
    ) +
    scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits = c(0, 3)) +
    scale_color_manual(values = cols) +
    # coord_cartesian(ylim = c(0, 3)) +
    labs(
      x = "",
      y = "Hazard ratio (95% CI)"
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

# sex
sgs <- quick_plot(sglt2)
dps <- quick_plot(dpp4)
ssgs <- quick_plot(sema_sglt2)
sdps <- quick_plot(sema_dpp4)

sex_plots <- wrap_plots(
  sgs +
    labs(title = "A. GLP-1 RA vs. SGLT2i", y = "") +
    theme(axis.text.x = element_blank()),
  dps +
    labs(title = "B. GLP-1 RA vs. DPP4i", y = "") +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()),
  ssgs + labs(title = "C. Semaglutide vs. SGLT2i"),
  sdps +
    labs(title = "D. Semaglutide vs. DPP4i") +
    theme(axis.text.y = element_blank()),
  ncol = 2,
  nrow = 2,
  guides = "collect"
) &
  theme(legend.position = "bottom")
ggsave(
  filename = glue(
    "{project_directory}figures/sex_stratified_plots_{version}.pdf"
  ),
  plot = sex_plots,
  width = 8,
  height = 10,
  device = cairo_pdf
)

# race/ethnicity
sgr <- quick_plot(sglt2, strata = c("NHW", "NHB", "HISP"))
dpr <- quick_plot(dpp4, strata = c("NHW", "NHB", "HISP"))
ssgr <- quick_plot(sema_sglt2, strata = c("NHW", "NHB", "HISP"))
sdpr <- quick_plot(sema_dpp4, strata = c("NHW", "NHB", "HISP"))

race_eth_plots <- wrap_plots(
  sgr +
    labs(title = "A. GLP-1 RA vs. SGLT2i", y = "") +
    theme(axis.text.x = element_blank()),
  dpr +
    labs(title = "B. GLP-1 RA vs. DPP4i", y = "") +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank()),
  ssgr + labs(title = "C. Semaglutide vs. SGLT2i"),
  sdpr +
    labs(title = "D. Semaglutide vs. DPP4i") +
    theme(axis.text.y = element_blank()),
  ncol = 2,
  nrow = 2,
  guides = "collect"
) &
  theme(legend.position = "bottom")
ggsave(
  filename = glue(
    "{project_directory}figures/race_eth_stratified_plots_{version}.pdf"
  ),
  plot = race_eth_plots,
  width = 8,
  height = 10,
  device = cairo_pdf
)
