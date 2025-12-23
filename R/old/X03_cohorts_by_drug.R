# libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
    library(arrow)
    library(cli)
    library(glue)
    library(data.table)
    library(PheWAS)
    library(dplyr)
    library(purrr)
    library(ggVennDiagram)
})

options(datatable.print.class = TRUE)
options(datatable.print.trunc.cols = TRUE)

project_directory   <- "/home/jupyter/workspaces/glp1rapharmacogenomics"
data_directory      <- file.path(project_directory, "data")
processed_directory <- file.path(data_directory, "processed")
version             <- "20250129"

dir.create(processed_directory, showWarnings = FALSE)

start_date <- as.Date("2018-01-01")

# data -------------------------------------------------------------------------
t2d <- read_parquet(
    file.path(processed_directory, paste0("t2d_demo_covariates_", version, ".parquet"))
)

drug_codes  <- fread("tables/all_drug_codes.csv")[class %in% c("sglt2", "dpp4", "glp1"), ]  
glp1_codes  <- drug_codes[class == "glp1", concept_id]
semaglutide_codes <- drug_codes[search_term == "semaglutide", concept_id]
sglt2_codes <- drug_codes[class == "sglt2", concept_id]
dpp4_codes  <- drug_codes[class == "dpp4", concept_id]

drug <- open_dataset(file.path(data_directory, "drug")) |>
    filter(id %in% t2d[, id] & drug_concept_id %in% drug_codes[, concept_id]) |>
    rename(date = drug_exposure_start_datetime) |>
    collect()
setDT(drug)
drug[, date := as.Date(date)]
setkeyv(drug, c("id", "date"))

drug <- merge(
    drug, t2d[, .(id, t2d_date, contra_date)], by = "id", all.x = TRUE
)

drug[, min_date := min(date), by = id]
drug[min_date  < contra_date]

first_glp1 <- drug[drug_concept_id %in% glp1_codes, .(
        min_date = min(date),
        max_date = max(date),
        count    = .N
    ),by = id]

drug

drug[(is.na(contra_date) | date < contra_date) & date >= min_date & date > t2d_date, ]

t2d[id %in% drug[date >= min_date, id], ]

# individuals with t2d with demographics (t2d)
t2d[, .N] # 52641

# individuals with glp1ra/sglt2i/dpp4i exposure, anytime (drug)
drug[, uniqueN(id)] # 14329

# individuals in both datasets (drug2)
drug2 <- merge(drug, t2d[, .(id, t2d_date, contra_date)], by = "id", all = FALSE)
drug2[, uniqueN(id)] # 14329

# individuals with drug exposure after t2d diagnosis (post_t2d)
post_t2d <- drug2[date >= t2d_date, ]
post_t2d[, uniqueN(id)] # 13891

# individuals with drug exposure after t2d diagnosis and before contraindication (post_t2d_pre_contra)
post_t2d_pre_contra <- post_t2d[date < contra_date | is.na(contra_date), ]
post_t2d_pre_contra[, uniqueN(id)] # 12681

# individuals with drug exposure after t2d diagnosis and before contraindication prescription since 2018-01-01 (post_t2d_pre_contra_since_2018)
post_t2d_pre_contra_since_2018 <- post_t2d_pre_contra[date >= start_date, ]
post_t2d_pre_contra_since_2018[, uniqueN(id)] # 11045

glp1  <- drug2[id %in% post_t2d_pre_contra_since_2018[, id] & drug_concept_id %in% setdiff(glp1_codes, semaglutide_codes), ]
sema  <- drug2[id %in% post_t2d_pre_contra_since_2018[, id] & drug_concept_id %in% semaglutide_codes, ]
sglt2 <- drug2[id %in% post_t2d_pre_contra_since_2018[, id] & drug_concept_id %in% sglt2_codes, ]
dpp4  <- drug2[id %in% post_t2d_pre_contra_since_2018[, id] & drug_concept_id %in% dpp4_codes, ]

glp1_ids <- glp1[, unique(id)]
sema_ids <- sema[, unique(id)]
sglt2_ids <- sglt2[, unique(id)]
dpp4_ids <- dpp4[, unique(id)]

venn <- list(
    glp1 = glp1_ids,
    sema = sema_ids,
    sglt2 = sglt2_ids,
    dpp4 = dpp4_ids
) |>
    ggVennDiagram(force_upset = TRUE, order.set.by = "name", order.intersect.by = "none")

ggsave(
    plot = venn,
    filename = file.path(project_directory, "plots", paste0("glp1_sema_sglt2_dpp4_venn_", version, ".pdf")),
    width = 10,
    height = 10,
    device = cairo_pdf
)

# filter individuals who have been exposed to an alternate drug class within 2 years prior to the first exposure to the drug of interest
