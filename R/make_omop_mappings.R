# libraaries -------------------------------------------------------------------
library(data.table)
library(dplyr)
library(arrow)
library(stringr)

box_directory <- file.path(
    "/Users",
    "maxsalvatore",
    "Library",
    "CloudStorage",
    "Box-Box",
    "projects",
    "GLP-1 RA",
    "PheWAS"
)

athena <- open_dataset(file.path(box_directory, "data", "athena")) |>
    collect()
setDT(athena)

# comorbidities
comorb <- fread(file.path(
    box_directory,
    "tables",
    "charlson_comorbidity_codes.csv"
))
setnames(comorb, "code", "concept_code")

icd <- athena[vocabulary_id %in% c("ICD9CM", "ICD10CM"), ]

comorb_omop <- merge(
    comorb,
    icd[, .(concept_id, concept_code, vocabulary_id, concept_name)],
    by = c("vocabulary_id", "concept_code")
)

comorb_omop <- comorb_omop[, .(
    study_name = comorbidity,
    concept_name,
    concept_id,
    vocabulary_id,
    concept_code
)][, group := "comorbidity"][]

# type 2 diabetes
t2d <- fread(
    "https://raw.githubusercontent.com/PheWAS/PhecodeX/refs/heads/main/phecodeX_R_map.csv",
    showProgress = FALSE
)[phecode == "EM_202.2", ]

t2d_omop <- merge(
    t2d[, .(concept_code = code, vocabulary_id)],
    icd[, .(concept_id, concept_code, vocabulary_id, concept_name)],
    by = c("vocabulary_id", "concept_code")
)[, study_name := "Type 2 diabetes"][]

t2d_omop <- t2d_omop[, .(
    study_name,
    concept_name,
    concept_id,
    vocabulary_id,
    concept_code
)][, group := "t2d"][]

# drugs
drug <- athena[vocabulary_id %in% c("ATC", "RxNorm"), ]

drug_strings <- c(
    "glp1|glp-1",
    "semaglutide",
    "dulaglutide",
    "exenatide",
    "liraglutide",
    "lixisenatide",
    "albiglutide",
    "beinaglutide",
    "tirzepatide",
    "dpp4|dpp-4",
    "sitagliptin",
    "saxagliptin",
    "linagliptin",
    "alogliptin",
    "vildagliptin",
    "gemigliptin",
    "teneligliptin",
    "sglt2|sglt-2",
    "dapagliflozin",
    "canagliflozin",
    "empagliflozin",
    "ertugliflozin",
    "bexagliflozin",
    "luseogliflozin",
    "ipragliflozin",
    "sotagliflozin"
)

drugs_omop <- rbindlist(lapply(
    drug_strings,
    \(x) {
        drug[grepl(x, concept_name, ignore.case = TRUE), ][, `:=`(
            study_name = str_to_sentence(str_extract(concept_name, "[^;]*")),
            search_term = x
        )]
    }
))[,
    drug_class := fcase(
        search_term %in%
            c(
                "semaglutide",
                "dulaglutide",
                "exenatide",
                "liraglutide",
                "lixisenatide",
                "albiglutide",
                "beinaglutide",
                "tirzepatide"
            ),
        "GLP-1 RA",
        search_term %in%
            c(
                "sitagliptin",
                "saxagliptin",
                "linagliptin",
                "alogliptin",
                "vildagliptin",
                "gemigliptin",
                "teneligliptin"
            ),
        "DPP4i",
        search_term %in%
            c(
                "dapagliflozin",
                "canagliflozin",
                "empagliflozin",
                "ertugliflozin",
                "bexagliflozin",
                "luseogliflozin",
                "ipragliflozin",
                "sotagliflozin"
            ),
        "SGLT2i"
    )
][]

drugs_omop <- drugs_omop[, .(
    study_name,
    concept_name,
    concept_id,
    vocabulary_id,
    concept_code,
    drug_class,
    search_term
)][, group := "drug"][]

# contraindications
map <- fread(
    "https://raw.githubusercontent.com/PheWAS/PhecodeX/refs/heads/main/phecodeX_R_map.csv",
    showProgress = FALSE
)
labs <- fread(
    "https://raw.githubusercontent.com/PheWAS/PhecodeX/refs/heads/main/phecodeX_R_labels.csv",
    showProgress = FALSE
)
icd_loinc <- athena[vocabulary_id %in% c("ICD10CM", "ICD9CM", "LOINC"), ]

contra_codes <- c(
    "MEN type II" = list(c("GE_961.22", "GE_961.23")),
    "Gastroparesis" = list("GI_516.4"),
    "Dialysis or kidney transplant" = list(c("GU_582.3", "SS_847.2"))
)

contra_omop <- lapply(
    seq_along(contra_codes),
    \(x) {
        data.table(
            study_name = names(contra_codes)[x],
            phecode = unlist(contra_codes[[x]])
        )
    }
) |>
    rbindlist() |>
    merge(
        map,
        by = "phecode",
        all.x = TRUE
    ) |>
    merge(
        labs[, .(phecode = phenotype, phecode_description = description)],
        by = "phecode",
        all.x = TRUE
    ) |>
    list(
        data.table(
            study_name = c(
                "Medullary thyroid cancer",
                "Medullary thyroid cancer",
                "Hypoglycemia with coma",
                "Hypoglycemia with coma",
                "Hypoglycemia with coma",
                "eGFR"
            ),
            code = c("C73", "Z85.850", "E15", "E13.641", "E11.641", "77147-7"),
            vocabulary_id = c(rep("ICD10CM", 5), "LOINC")
        )
    ) |>
    rbindlist(use.names = TRUE, fill = TRUE) |>
    merge(
        icd_loinc[, .(
            concept_id,
            concept_name,
            code = concept_code,
            vocabulary_id
        )],
        by = c("code", "vocabulary_id"),
        all.x = TRUE
    )

contra_omop <- contra_omop[, .(
    study_name,
    concept_name,
    concept_id,
    vocabulary_id,
    concept_code = code,
    group = "contraindication"
)]

# labs -------------------------------------------------------------------------
lab_codes <- fread(file.path(box_directory, "tables", "labs_codes.csv"))

lab_omop <- merge(
    lab_codes,
    icd_loinc[, .(concept_id, concept_code, vocabulary_id, concept_name)],
    by = c("vocabulary_id", "concept_code")
)[, group := "lab"][]

# stack and save
omop <- rbindlist(
    list(comorb_omop, t2d_omop, drugs_omop, contra_omop, lab_omop),
    use.names = TRUE,
    fill = TRUE
)

fwrite(
    omop,
    file.path(box_directory, "tables", "omop_mappings.csv")
)
