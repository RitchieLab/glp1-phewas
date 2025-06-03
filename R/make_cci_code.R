library(data.table)

source(
    "https://raw.githubusercontent.com/ellessenne/comorbidity/refs/heads/master/data-raw/make-mapping.R"
)

cci_icd9_codes <- .maps[["charlson_icd9_quan"]]
cci_icd10_codes <- .maps[["charlson_icd10_quan"]]

nl_to_dt <- function(nl, vocab_id) {
    rbindlist(lapply(names(nl), function(name) {
        data.table(abbrev = name, code = nl[[name]])
    }))[, vocabulary_id := vocab_id][]
}

cci_icd9_tab <- nl_to_dt(cci_icd9_codes, "ICD9CM")
cci_icd10_tab <- nl_to_dt(cci_icd10_codes, "ICD10CM")

format_icd9 <- function(codes) {
    # Input validation
    if (!is.vector(codes) || !is.character(codes)) {
        stop("Input must be a character vector")
    }

    # Function to format a single ICD9 code
    format_single_code <- function(code) {
        # Remove any whitespace and convert to uppercase
        code <- toupper(trimws(code))

        # Skip NA or empty values
        if (is.na(code) || code == "") {
            return(NA_character_)
        }

        # Remove any dots or other separators
        clean_code <- gsub("[^A-Z0-9V]", "", code)

        # Check if the code is empty after cleaning
        if (clean_code == "") {
            return(NA_character_)
        }

        # Check if we have a V code (V01-V91)
        if (substr(clean_code, 1, 1) == "V") {
            if (nchar(clean_code) <= 3) {
                # Short V codes don't need dots
                return(clean_code)
            } else {
                # Add dot after the third character for V codes
                return(paste0(
                    substr(clean_code, 1, 3),
                    ".",
                    substr(clean_code, 4, nchar(clean_code))
                ))
            }
        }

        # Check if we have an E code (E800-E999)
        if (substr(clean_code, 1, 1) == "E") {
            if (nchar(clean_code) <= 4) {
                # Short E codes don't need dots
                return(clean_code)
            } else {
                # Add dot after the fourth character for E codes
                return(paste0(
                    substr(clean_code, 1, 4),
                    ".",
                    substr(clean_code, 5, nchar(clean_code))
                ))
            }
        }

        # For numeric codes (000-999)
        if (nchar(clean_code) <= 3) {
            # Codes with 3 or fewer digits don't need dots
            return(clean_code)
        } else {
            # Add dot after the third digit
            return(paste0(
                substr(clean_code, 1, 3),
                ".",
                substr(clean_code, 4, nchar(clean_code))
            ))
        }
    }

    # Apply formatting to each code in the vector
    sapply(codes, format_single_code)
}

cci_icd9_tab[, code := format_icd9(code)]

format_icd10 <- function(codes) {
    # Input validation
    if (!is.vector(codes) || !is.character(codes)) {
        stop("Input must be a character vector")
    }

    # Function to format a single ICD10 code
    format_single_code <- function(code) {
        # Remove any whitespace and convert to uppercase
        code <- toupper(trimws(code))

        # Skip NA or empty values
        if (is.na(code) || code == "") {
            return(NA_character_)
        }

        # Remove any dots or other separators
        clean_code <- gsub("[^A-Z0-9]", "", code)

        # Check if the code matches the ICD10 pattern (letter followed by digits)
        if (!grepl("^[A-Z][0-9]", clean_code)) {
            warning(paste(
                "Code",
                code,
                "does not appear to be a valid ICD10 code format"
            ))
            return(clean_code)
        }

        # Format: Letter followed by 2 digits, then dot, then remaining digits
        # Example: I22 -> I22, I252 -> I25.2
        letter <- substr(clean_code, 1, 1)
        nums <- substr(clean_code, 2, nchar(clean_code))

        if (nchar(nums) <= 2) {
            # If we have 2 or fewer digits, no dot needed
            formatted <- clean_code
        } else {
            # For 3+ digits, add dot after first 2 digits
            formatted <- paste0(
                letter,
                substr(nums, 1, 2),
                ".",
                substr(nums, 3, nchar(nums))
            )
        }

        return(formatted)
    }

    # Apply formatting to each code in the vector
    sapply(codes, format_single_code)
}

cci_icd10_tab[, code := format_icd10(code)]

cci_tab <- rbindlist(
    list(cci_icd9_tab, cci_icd10_tab)
)[, `:=`(
    comorbidity = fcase(
        abbrev == "mi",
        "Myocardial infarction",
        abbrev == "chf",
        "Congestive heart failure",
        abbrev == "pvd",
        "Peripheral vascular disease",
        abbrev == "cevd",
        "Cerebrovascular disease",
        abbrev == "dementia",
        "Dementia",
        abbrev == "cpd",
        "Chronic pulmonary disease",
        abbrev == "rheumd",
        "Rheumatic disease",
        abbrev == "pud",
        "Peptic ulcer disease",
        abbrev == "mld",
        "Mild liver disease",
        abbrev == "diab",
        "Diabetes without chronic complications",
        abbrev == "diabwc",
        "Diabetes with chronic complications",
        abbrev == "hp",
        "Hemiplegia or paraplegia",
        abbrev == "rend",
        "Renal disease",
        abbrev == "canc",
        "Any malignancy, including lymphoma and leukemia, except malignant neoplasm of skin",
        abbrev == "msld",
        "Moderate or severe liver disease",
        abbrev == "metacanc",
        "Metastatic solid tumor",
        abbrev == "aids",
        "AIDS/HIV"
    )
)][]

fwrite(
    cci_tab[, .(comorbidity, abbrev, code, vocabulary_id)],
    "/Users/maxsalvatore/Library/CloudStorage/Box-Box/projects/GLP-1 RA/PheWAS/tables/charlson_comorbidity_codes.csv"
)
