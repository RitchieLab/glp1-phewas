source(
    "https://gitlab.com/maxsal/plot_phewasx/-/raw/main/plot_phewasx.R?ref_type=heads"
)

plot_rmst_phewasx <- function(
    data,
    pheinfox = NULL,
    rmst_var = "rmst_est_3",
    p_var = "rmst_p_3",
    color_var = "color",
    phe_var = "phecode",
    order_reset = TRUE,
    remove_insig = FALSE,
    rmst_diff = 30,
    rmst_diff_type = "rect",
    rect_col = "gray90",
    y_breaks = NULL,
    annotate = FALSE,
    annotate_label = "treatment",
    group_space = 20,
    label_top = 5,
    title = NULL,
    genetic_offset = 15,
    color_dot_pt = 15,
    color_dot_symbol = "\u25CF",
    label_size = 3,
    base_theme = ggplot2::theme_minimal(),
    ...
) {
    # initialize
    if (is.null(pheinfox)) pheinfox <- .clean_phecode_labels()

    # prep
    if (!is.data.table(data)) data <- data.table::as.data.table(data)
    data2 <- data.table::copy(data)
    if ((rmst_var %in% names(data2)) & (p_var %in% names(data2))) {
        phewas <- data2[,
            direction := data.table::fcase(
                get(rmst_var) > 0 & get(p_var) < 0.05,
                "Positive",
                get(rmst_var) < 0 & get(p_var) < 0.05,
                "Negative",
                default = "None"
            )
        ]
        if (remove_insig) {
            phewas <- phewas[direction != "None", ]
        }
    } else {
        phewas <- data2[, direction := "None"]
    }

    plot_data <- merge.data.table(
        phewas,
        pheinfox[, .(
            phecode,
            group,
            groupnum,
            description,
            order
        )],
        by.x = phe_var,
        by.y = "phecode"
    )[order(order), ]

    # reset ordering, if requested
    if (order_reset) {
        plot_data[order(order), order := 1:.N]
    }

    # add spacing, if requested
    if (group_space > 0) {
        plot_data <- plot_data[,
            order := order + ((groupnum - 1) * group_space)
        ]
    }
    if (genetic_offset > 0) {
        plot_data <- plot_data[
            group == "Genetic",
            order := order + genetic_offset
        ]
    }

    vars <- c("group", color_var)
    plot_data_mean <- plot_data[,
        .(mean = mean(order, na.rm = TRUE)),
        by = group
    ]
    plot_data_mean <- merge.data.table(
        plot_data_mean,
        unique(pheinfox[, ..vars]),
        by = "group"
    )
    phe_colors <- plot_data_mean[[color_var]]
    names(phe_colors) <- plot_data_mean[, group]

    # plot
    plot <- plot_data |>
        ggplot2::ggplot(aes(
            x = order,
            y = get(rmst_var),
            fill = group,
            color = group
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
                    color = NA,
                    alpha = 0.25,
                    show.legend = FALSE
                )
            })() +
        (\(x)
            if (annotate) {
                annotate(
                    "text",
                    x = max(plot_data[, order]),
                    y = max(abs(plot_data[[rmst_var]])) * 0.2,
                    label = paste0("Decreased risk with ", annotate_label),
                    hjust = 1,
                    size = 5,
                    color = "gray50"
                )
            })() +
        (\(x)
            if (annotate) {
                annotate(
                    "text",
                    x = max(plot_data[, order]),
                    y = -max(abs(plot_data[[rmst_var]])) * 0.2,
                    label = paste0("Increased risk with ", annotate_label),
                    hjust = 1,
                    size = 5,
                    color = "gray50"
                )
            })() +
        ggplot2::geom_point(
            ggplot2::aes(shape = direction),
            size = 2,
            alpha = 0.5,
            show.legend = FALSE
        ) +
        ggplot2::geom_hline(
            yintercept = 0,
            linewidth = 1
        ) +
        (\()
            if (!is.null(color_dot_symbol)) {
                ggplot2::scale_x_continuous(
                    breaks = plot_data_mean[, mean],
                    labels = paste0(
                        "<span style=\"color: ",
                        plot_data_mean[, color],
                        "\"><span style=\"font-size: ",
                        color_dot_pt,
                        "pt\">",
                        color_dot_symbol,
                        "</span></span> ",
                        plot_data_mean[, group]
                    )
                )
            } else {
                ggplot2::scale_x_continuous(
                    breaks = plot_data_mean[, mean],
                    labels = paste0(
                        "<span style=\"color: ",
                        plot_data_mean[, color],
                        "\">",
                        plot_data_mean[, group],
                        "</span>"
                    )
                )
            })() +
        (\(x)
            if (!is.null(y_breaks)) {
                ggplot2::scale_y_continuous(
                    breaks = y_breaks
                )
            })() +
        ggplot2::labs(
            x = "",
            y = "3-year RMST difference, days"
        ) +
        ggplot2::scale_fill_manual(values = phe_colors) +
        ggplot2::scale_color_manual(values = phe_colors) +
        ggplot2::scale_shape_manual(
            values = c("Positive" = 24, "Negative" = 25, "None" = 19)
        ) +
        theme_phewas(base_theme = base_theme, ...) +
        ggplot2::theme(
            axis.text.x = ggtext::element_markdown(angle = -45, hjust = 0)
        )

    if (!is.null(label_top)) {
        plot <- plot +
            ggrepel::geom_label_repel(
                data = plot_data[order(-abs(get(rmst_var))), ][1:label_top, ],
                ggplot2::aes(
                    x = order,
                    y = get(rmst_var),
                    label = description
                ),
                color = "black",
                show.legend = FALSE,
                inherit.aes = FALSE,
                size = label_size,
                max.time = 3,
                max.iter = 100000
            )
    }

    if (!is.null(title)) {
        plot <- plot + ggplot2::labs(title = title)
    }

    return(plot)
}
