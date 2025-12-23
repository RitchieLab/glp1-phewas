# load packages
library(tidyverse)
library(quarto)
library(data.table)

template_dir <- "/Users/maxsalvatore/Library/CloudStorage/Box-Box/projects/GLP-1 RA/PheWAS/reports/"
output_dir <- "/Users/maxsalvatore/Library/CloudStorage/Box-Box/projects/GLP-1 RA/PheWAS/reports/"

# create combinations
grid <- data.table(expand.grid(
    treatment_name = c("GLP-1 RA", "Semaglutide"),
    comparator_name = c("SGLT2i", "DPP4i"),
    type_name = c("intention-to-treat", "per-protocol")
)) |>
    mutate_if(is.factor, as.character)

to_file_name <- function(x) {
    paste0(
        fcase(
            x$treatment_name == "GLP-1 RA",
            "glp1",
            x$treatment_name == "Semaglutide",
            "sema"
        ),
        "_",
        fcase(
            x$comparator_name == "SGLT2i",
            "sglt2",
            x$comparator_name == "DPP4i",
            "dpp4"
        ),
        "_",
        fcase(
            x$type_name == "intention-to-treat",
            "itt",
            x$type_name == "per-protocol",
            "pp"
        ),
        "_phewas_report.pdf"
    )
}

grid_list <- lapply(split(grid, seq(nrow(grid))), \(x) as.list(x))

reports <- tibble(
    input = paste0(template_dir, "phewas_report.qmd"),
    output_file = unlist(lapply(grid_list, \(y) to_file_name(y))),
    execute_params = unname(grid_list)
)

reports |>
    pwalk(quarto_render, .progress = TRUE)

# # for testing
# quarto_render(
#     input = paste0(template_dir, "phewas_report.qmd"),
#     output_file = "phewas_report.pdf",
#     execute_params = list(
#         treatment_name = "GLP-1 RA",
#         comparator_name = "SGLT2i",
#         type_name = "intention-to-treat"
#     )
# )
