# libraries and such -----------------------------------------------------------
suppressMessages({
    library(data.table)
    library(tidyverse)
    library(arrow)
    library(cli)
    library(patchwork)
})

options(datatable.print.trunc.cols = TRUE)
options(datatable.print.class = TRUE)

p_dir <- "/home/jupyter/workspaces/duplicateofglp1rapharmacogenomics/"
# setwd(p_dir)

# load data --------------------------------------------------------------------
## analytic dataset
analytic_dataset <- read_parquet(
    paste0(p_dir, "data/processed/analytic_dataset_20250129.parquet")
)

## first occurrence of EM_236: "Overweight or obesity"
phe <- read_parquet(
    paste0(p_dir, "data/processed/phenome_first_20250129.parquet")
)[person_id %in% analytic_dataset[, unique(person_id)], ]
obese_phe <- phe[phecode == "EM_236.11", .(person_id, phe_date = date)]

## load bmi dataset
bmi <- open_dataset(
    paste0(p_dir, "data/raw/bmi/")
) |>
    mutate(person_id = as.character(person_id)) |>
    filter(
        standard_concept_name == "Body mass index (BMI) [Ratio]" &
            !is.na(value_as_number) &
            person_id %in% analytic_dataset[, unique(person_id)]
    ) |>
    select(
        person_id,
        bmi = value_as_number,
        date = measurement_datetime
    ) |>
    collect() |>
    as.data.table()
# convert datetime to date
bmi[, date := as.Date(date)]
# remove duplicates
bmi <- unique(bmi)

# part 1: BMI diagnosis relative to measurement --------------------------------
first_obese_bmi <- bmi[bmi >= 40, .SD[which.min(date)], by = person_id][, .(
    person_id,
    bmi,
    bmi_date = date
)]

meas_vs_dx <- merge(
    first_obese_bmi,
    obese_phe,
    by = "person_id",
    all.x = TRUE
)[, diff_days := as.numeric(phe_date - bmi_date)][]

meas_vs_dx[,
    diff_days_cat := fcase(
        diff_days < 0,
        "Diagnosis before measurement",
        diff_days == 0,
        "Diagnosis same day as measurement",
        diff_days > 0,
        "Diagnosis after measurement",
        default = NA
    )
]

meas_vs_dx[!is.na(diff_days), .N, diff_days_cat][, prop := N / sum(N)][order(
    -prop
)]

# plot of time from first obese measurement to diagnosis
obese_dx_plot <- meas_vs_dx |>
    ggplot(aes(x = diff_days)) +
    annotate(
        "text",
        label = "Phecode before severely obese\nBMI measurement",
        x = -5000,
        y = 250,
        hjust = 0.5,
        vjust = 0.5,
        size = 4,
        color = "gray60"
    ) +
    annotate(
        "text",
        label = "Phecode after severely obese\nBMI measurement",
        x = 5000,
        y = 250,
        hjust = 0.5,
        vjust = 0.5,
        size = 4,
        color = "gray60"
    ) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_x_continuous(
        breaks = seq(-7500, 7500, by = 2500),
        labels = scales::comma
    ) +
    scale_y_continuous(
        labels = scales::comma
    ) +
    coord_cartesian(xlim = c(-7500, 7500)) +
    labs(
        title = stringr::str_wrap(
            "Time from first severely obese BMI measurement to first diagnosis of morbid obesity",
            width = 50
        ),
        x = "Phecode date - BMI date (days)",
        y = "Count"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )
ggsave(
    filename = paste0(p_dir, "figures/obese_dx_plot.pdf"),
    plot = obese_dx_plot,
    width = 7,
    height = 5,
    device = cairo_pdf
)

min_phe <- phe[phe[, .I[which.min(date)], person_id]$V1, ]
bmi_idx <- bmi[, .I[which.min(date)], person_id]$V1
min_bmi <- bmi[bmi_idx, ]

min_dt <- merge(
    min_phe[, .(person_id, phe_date = date)],
    min_bmi[, .(person_id, bmi_date = date)],
    by = "person_id"
)
min_dt[, diff_years := round(as.numeric(phe_date - bmi_date) / 365.25, 1)]

# plot of time from first phecode to first measurement
bmi_dx_plot <- min_dt |>
    ggplot(aes(x = diff_years)) +
    annotate(
        "text",
        label = "Phecode before\nBMI measurement",
        x = -25,
        y = 1250,
        hjust = 0.5,
        vjust = 0.5,
        size = 4,
        color = "gray60"
    ) +
    annotate(
        "text",
        label = "Phecode after\nBMI measurement",
        x = 12.5,
        y = 1250,
        hjust = 0.5,
        vjust = 0.5,
        size = 4,
        color = "gray60"
    ) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    # scale_x_continuous(
    #     breaks = seq(-40, 10, by = 10)
    # ) +
    labs(
        title = stringr::str_wrap(
            "Time from first  BMI measurement to first phecode",
            width = 50
        ),
        x = "Phecode date - BMI date (years)",
        y = "Count"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )
ggsave(
    filename = paste0(p_dir, "figures/bmi_dx_plot.pdf"),
    plot = bmi_dx_plot,
    width = 7,
    height = 5,
    device = cairo_pdf
)

set.seed(2025)
# restrict to measurements after 2017-01-01
bmi <- bmi[date >= as.Date("2017-01-01"), ]
# sample three individuals with at least 5 BMI measurements
bmi_traj_dt <- merge(
    bmi[, .(person_id, bmi_meas = bmi, bmi_date = date)],
    analytic_dataset[
        drug_class == "GLP-1 RA",
        .(person_id, glp1_date = glp1_min_date)
    ],
    by = "person_id",
    all = FALSE
)
bmi_traj_dt[, diff_days := as.numeric(bmi_date - glp1_date)]

phe_full <- open_dataset(
    paste0(p_dir, "data/processed/phenome_first_20250129.parquet")
) |>
    filter(
        person_id %in%
            analytic_dataset[, unique(person_id)] &
            phecode == "EM_236.11"
    ) |>
    collect() |>
    merge(
        analytic_dataset[
            drug_class == "GLP-1 RA",
            .(person_id, glp1_date = glp1_min_date)
        ],
        by = "person_id",
        all.x = TRUE
    ) |>
    mutate(
        diff_days = as.numeric(date - glp1_date)
    )

post_ids <- phe_full[
    !is.na(glp1_date) &
        diff_days > 0 &
        person_id %in%
            bmi_traj_dt[between(diff_days, -365, 365), .N, person_id][
                N >= 3,
                person_id
            ],
    sample(person_id, 9)
]
phe_full[person_id %in% post_ids, ]

bmi_traj_plot <- bmi_traj_dt[person_id %in% post_ids, ] |>
    ggplot(aes(x = diff_days, y = bmi_meas, color = person_id)) +
    geom_rect(
        aes(
            xmin = -Inf,
            xmax = Inf,
            ymin = 40,
            ymax = Inf
        ),
        color = NA,
        fill = "gray95",
        alpha = 0.3
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1) +
    geom_vline(
        aes(xintercept = diff_days, color = person_id),
        data = phe_full[person_id %in% post_ids, ]
    ) +
    geom_line(alpha = 0.5, linewidth = 1) +
    geom_point(alpha = 0.5) +
    labs(
        title = stringr::str_wrap(
            "BMI trajectory before and after GLP-1 RA initiation",
            width = 50
        ),
        subtitle = "Vertical lines indicate morbid obesity phecode diagnosis",
        x = "Days from GLP-1 RA initiation",
        y = "BMI (kg/m^2)"
    ) +
    ggokabeito::scale_color_okabe_ito() +
    coord_cartesian(xlim = c(-180, 365), ylim = c(18.5, NA)) +
    theme_minimal() +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
    )
ggsave(
    filename = paste0(p_dir, "figures/bmi_traj_plot.pdf"),
    plot = bmi_traj_plot,
    width = 7,
    height = 5,
    device = cairo_pdf
)

samp <- bmi[, .N, person_id][N > 5, sample(person_id, 3)]

phe[person_id %in% samp, ]
