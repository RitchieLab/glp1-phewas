source(
  "https://gitlab.com/maxsal/plot_phewasx/-/raw/main/plot_phewasx.R?ref_type=heads"
)

#' Prepare data for PhecodeX PheWAS RMST plot
#'
#' This function handles all data transformations (merging with phecode metadata,
#' reordering, spacing, etc.) independently of plotting. This enables:
#' - Reproducible extraction of minimal plot data
#' - Separation of data prep from visualization
#' - Alternative visualizations using the same prepared data
#'
#' @param data data of RMST PheWAS results
#' @param phe_var name of column containing phecode values
#' @param rmst_var name of column containing RMST estimates
#' @param p_var name of column containing RMST p-values
#' @param color_var name of column to store color values (default: "color")
#' @param pheinfox data of PheWAS codes, descriptions, groups, etc.
#' @param order_reset reset ordering of phecodes (TRUE/FALSE, default: TRUE)
#' @param remove_insig remove insignificant results (TRUE/FALSE, default: FALSE)
#' @param group_space space between phenotype groups (default: 20)
#' @param genetic_offset offset for genetic phecodes (default: 15)
#'
#' @return A data.table with prepared columns: phecode, order, group, groupnum,
#'         description, RMST values, p-values, color, and direction
#'
#' @importFrom data.table copy
#' @importFrom data.table as.data.table
#' @importFrom data.table is.data.table
#' @importFrom data.table fcase
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare data for plotting
#' phewasx_rmst_data <- phewasx_rmst_plot_data(
#'     data = my_rmst_results,
#'     phe_var = "phecode",
#'     rmst_var = "rmst_est_3",
#'     p_var = "rmst_p_3"
#' )
#'
#' # Save for reproducibility
#' fwrite(phewasx_rmst_data, "plot_data_rmst_minimal.csv")
#'
#' # Use in plot
#' plot_rmst_phewasx(phewasx_rmst_data, order_reset = FALSE, group_space = 0, genetic_offset = 0)
#' }

phewasx_rmst_plot_data <- function(
    data,
    phe_var        = "phecode",
    rmst_var       = "rmst_est_3",
    p_var          = "rmst_p_3",
    color_var      = "color",
    pheinfox       = NULL,
    order_reset    = TRUE,
    remove_insig   = FALSE,
    group_space    = 20,
    genetic_offset = 15
) {
    # Initialize pheinfox if not provided
    if (is.null(pheinfox)) pheinfox <- .clean_phecode_labels()

    # prep
    if (!data.table::is.data.table(data)) data <- data.table::as.data.table(data)
    data2 <- data.table::copy(data)

    # Create direction variable (with significance threshold)
    if ((rmst_var %in% names(data2)) & (p_var %in% names(data2))) {
        phewas <- data2[,
            direction := data.table::fcase(
                get(rmst_var) > 0 & get(p_var) < 0.05, "Positive",
                get(rmst_var) < 0 & get(p_var) < 0.05, "Negative",
                default = "None"
            )
        ]
        if (remove_insig) {
            phewas <- phewas[direction != "None", ]
        }
    } else {
        phewas <- data2[, direction := "None"]
    }

    # Merge with phecode metadata
    plot_data <- merge.data.table(
        phewas,
        pheinfox[, .SD, .SDcols = c(
            "phecode",
            "group",
            "groupnum",
            "description",
            "order",
            "color"
        )],
        by.x = phe_var,
        by.y = "phecode",
        all.x = TRUE
    )[order(order), ]

    # Reset ordering if requested
    if (order_reset) {
        plot_data[order(order), order := 1:.N]
    }

    # Add group spacing if requested
    if (group_space > 0) {
        plot_data[, order := order + ((groupnum - 1) * group_space)]
    }

    # Add genetic offset if requested
    if (genetic_offset > 0) {
        plot_data[group == "Genetic", order := order + genetic_offset]
    }

    # Return selected columns
    plot_data[, .SD, .SDcols = c(
        phe_var,
        "order",
        "group",
        "groupnum",
        "description",
        rmst_var,
        p_var,
        "color",
        "direction"
    )]
}

#' Plot PhecodeX PheWAS RMST results
#'
#' @param data data of RMST PheWAS results (raw or pre-prepared with phewasx_rmst_plot_data())
#' @param pheinfox data of PheWAS codes, descriptions, groups, etc.
#' @param rmst_var name of column containing RMST estimates
#' @param p_var name of column containing RMST p-values
#' @param color_var name of column containing color values
#' @param phe_var name of column containing phecode values
#' @param order_reset reset ordering of phecodes
#' @param remove_insig remove insignificant results
#' @param rmst_diff RMST difference threshold for visualization
#' @param rmst_diff_type type of visualization: "rect" (rectangle) or "lines" (horizontal lines)
#' @param rect_col color of rectangle background
#' @param y_breaks custom y-axis breaks
#' @param annotate add text annotations for treatment effect direction
#' @param annotate_label label for treatment in annotations
#' @param group_space space between groups
#' @param label_top number of top phenotypes to label
#' @param title title of plot
#' @param genetic_offset offset for genetic phecodes
#' @param color_x_labels Colorize the x-axis text at all (TRUE [default]) or not (FALSE)
#' @param color_dot_pt size of color dots
#' @param color_dot_symbol symbol for color dots
#' @param label_size size of labels
#' @param base_theme base ggplot2 theme
#' @param return_plot_data logical; if TRUE, return prepared data instead of plot
#'
#' @return a ggplot2 object (or data.table if return_plot_data = TRUE)
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom data.table is.data.table
#' @importFrom data.table as.data.table
#' @importFrom data.table copy
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggtext element_markdown
#' @export
#' @examples
#' \dontrun{
#' # Plot raw RMST data
#' plot_rmst_phewasx(my_rmst_results)
#'
#' # Or prepare data first for reproducibility
#' prepared <- phewasx_rmst_plot_data(my_rmst_results)
#' plot_rmst_phewasx(prepared, order_reset = FALSE, group_space = 0, genetic_offset = 0)
#' }

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
    color_x_labels = TRUE,
    color_dot_pt = 15,
    color_dot_symbol = "\u25CF",
    label_size = 3,
    base_theme = ggplot2::theme_minimal(),
    return_plot_data = FALSE,
    ...
) {
    # Initialize pheinfox
    if (is.null(pheinfox)) pheinfox <- .clean_phecode_labels()

    # Prepare plot data
    plot_data <- phewasx_rmst_plot_data(
        data = data,
        phe_var = phe_var,
        rmst_var = rmst_var,
        p_var = p_var,
        color_var = color_var,
        pheinfox = pheinfox,
        order_reset = order_reset,
        remove_insig = remove_insig,
        group_space = group_space,
        genetic_offset = genetic_offset
    )

    # Return data if requested
    if (return_plot_data) return(plot_data)

    # Compute group means for x-axis labels
    plot_data_mean <- plot_data[,
        .(mean = mean(order, na.rm = TRUE)),
        by = group
    ]
    plot_data_mean <- merge.data.table(
        plot_data_mean,
        unique(pheinfox[, .SD, .SDcols = c("group", "color")]),
        by = "group"
    )
    phe_colors <- plot_data_mean[["color"]]
    names(phe_colors) <- plot_data_mean[["group"]]

    # Create plot
    plot <- plot_data |>
        ggplot2::ggplot(ggplot2::aes(
            x = order,
            y = get(rmst_var),
            fill = group,
            color = group
        )) +
        (\(x)
            if (rmst_diff_type == "lines") {
                ggplot2::geom_hline(
                    yintercept = c(-rmst_diff, rmst_diff),
                    linetype = 2,
                    color = "black"
                )
            })() +
        (\(x)
            if (grepl("rect", rmst_diff_type, ignore.case = TRUE)) {
                ggplot2::geom_rect(
                    ggplot2::aes(
                        xmin = -Inf,
                        xmax = Inf,
                        ymin = -rmst_diff,
                        ymax = rmst_diff
                    ),
                    fill = rect_col,
                    color = NA,
                    alpha = 0.25,
                    show.legend = FALSE,
                    inherit.aes = FALSE
                )
            })() +
        (\(x)
            if (annotate) {
                ggplot2::annotate(
                    "text",
                    x = max(plot_data[["order"]]),
                    y = max(abs(plot_data[[rmst_var]])) * 0.2,
                    label = paste0("Decreased risk with ", annotate_label),
                    hjust = 1,
                    size = 5,
                    color = "gray50"
                )
            })() +
        (\(x)
            if (annotate) {
                ggplot2::annotate(
                    "text",
                    x = max(plot_data[["order"]]),
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
            if (color_x_labels == TRUE) {
                if (!is.null(color_dot_symbol)) {
                    ggplot2::scale_x_continuous(
                        breaks = plot_data_mean[["mean"]],
                        labels = paste0(
                            "<span style=\"color: ",
                            plot_data_mean[["color"]],
                            "\"><span style=\"font-size: ",
                            color_dot_pt,
                            "pt\">",
                            color_dot_symbol,
                            "</span></span> ",
                            plot_data_mean[["group"]]
                        )
                    )
                } else {
                    ggplot2::scale_x_continuous(
                        breaks = plot_data_mean[["mean"]],
                        labels = paste0(
                            "<span style=\"color: ",
                            plot_data_mean[["color"]],
                            "\">",
                            plot_data_mean[["group"]],
                            "</span>"
                        )
                    )
                }
            } else {
                ggplot2::scale_x_continuous(
                    breaks = plot_data_mean[["mean"]],
                    labels = plot_data_mean[["group"]]
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

    # Add top labels if requested
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

    # Add title if provided
    if (!is.null(title)) {
        plot <- plot + ggplot2::labs(title = title)
    }

    return(plot)
}
