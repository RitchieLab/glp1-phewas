# libraries --------------------------------------------------------------------
library(data.table)
library(tidyverse)
library(ggtext)

project_directory <- "/Users/maxsalvatore/Library/CloudStorage/Box-Box/projects/GLP-1 RA/PheWAS/"

# data -------------------------------------------------------------------------
dt <- fread(file.path(project_directory, "data", "xie_et_al_2025_results.csv"))
ccsr_labs <- fread(file.path(
    project_directory,
    "tables",
    "ccsr_categories.csv"
))
setnames(dt, "organ_system", "group")
dt[, groupnum := as.integer(factor(group, levels = unique(group)))]
dt[order(groupnum, outcomes), order := 1:.N, by = comparator]

ccsr_labs <- unique(ccsr_labs[, .(abbrev, color)])
dt <- merge(dt, ccsr_labs, by = "abbrev", all.x = TRUE)

dt1 <- dt[comparator == "SGLT2i", ]

# plot function ----------------------------------------------------------------
plot_xie <- function(dt1) {
    dt1[,
        dir := fcase(
            est > 1,
            "Up",
            est < 1,
            "Down",
            default = "None"
        )
    ]
    dt1_cols <- unique(dt1[, .(group, color)])
    cols <- dt1_cols[, color]
    names(cols) <- dt1_cols[, group]

    dt1[, order := order + ((groupnum - 1) * 10)]

    dt1_breaks <- dt1[, .(x = mean(order)), group]
    dt1_breaks[, group := str_wrap(dt1_breaks[, group], 25)]

    dirs <- c(
        "Up" = 24,
        "Down" = 25,
        "None" = 21
    )

    dt1_top <- dt1[order(p), ][1:5, ]

    dt1 |>
        ggplot(aes(
            x = order,
            y = -log10(p),
            color = group,
            fill = group,
            shape = dir
        )) +
        geom_hline(
            yintercept = -log10(0.05),
            linetype = "dashed",
            color = "orange"
        ) +
        geom_hline(
            yintercept = -log10(0.05 / dt1[, .N]),
            linetype = "dashed",
            color = "red"
        ) +
        geom_point(size = 3, alpha = 0.5) +
        scale_color_manual(values = cols) +
        scale_fill_manual(values = cols) +
        scale_x_continuous(
            breaks = dt1_breaks[, x],
            labels = dt1_breaks[, group]
        ) +
        scale_shape_manual(values = dirs) +
        labs(
            x = ""
        ) +
        ggrepel::geom_label_repel(
            data = dt1_top,
            ggplot2::aes(
                x = order,
                y = -log10(p),
                label = outcomes
            ),
            color = "black",
            show.legend = FALSE,
            inherit.aes = FALSE,
            size = 3,
            max.time = 3,
            max.iter = 100000
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0, face = "bold"),
            axis.text.x = element_text(angle = -45, hjust = 0),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            legend.position = "none"
        )
}

# plot -------------------------------------------------------------------------
ggsave(
    file.path(
        project_directory,
        "figures",
        "glp1_sglt2_xie_plot.pdf"
    ),
    plot_xie(dt[comparator == "SGLT2i", ]),
    width = 7,
    height = 5,
    device = cairo_pdf
)
ggsave(
    file.path(
        project_directory,
        "figures",
        "glp1_dpp4_xie_plot.pdf"
    ),
    plot_xie(dt[comparator == "DPP4i", ]),
    width = 7,
    height = 5,
    device = cairo_pdf
)
