# libraries --------------------------------------------------------------------
source("R/config.R")
source("R/setup.R")
source("R/fn/helper-functions.R")
load_libraries(c(pkg_groups$core, pkg_groups$bigquery))

project_directory <- config$project_directory

dir.create(
  file.path(config$project_directory, config$paths$raw),
  showWarnings = FALSE,
  recursive = TRUE
)
dir.create(
  file.path(config$project_directory, config$paths$results),
  showWarnings = FALSE,
  recursive = TRUE
)

# query data -------------------------------------------------------------------
## demographics data -----------------------------------------------------------
cli_progress_step("querying demographic data")
demo_sql_query <- glue::glue(
  "
    SELECT
        person.person_id,
        person.gender_concept_id,
        p_gender_concept.concept_name as gender,
        person.birth_datetime as date_of_birth,
        person.race_concept_id,
        p_race_concept.concept_name as race,
        person.ethnicity_concept_id,
        p_ethnicity_concept.concept_name as ethnicity,
        person.sex_at_birth_concept_id,
        p_sex_at_birth_concept.concept_name as sex_at_birth 
    FROM
        `person` person 
    LEFT JOIN
        `concept` p_gender_concept 
            ON person.gender_concept_id = p_gender_concept.concept_id 
    LEFT JOIN
        `concept` p_race_concept 
            ON person.race_concept_id = p_race_concept.concept_id 
    LEFT JOIN
        `concept` p_ethnicity_concept 
            ON person.ethnicity_concept_id = p_ethnicity_concept.concept_id 
    LEFT JOIN
        `concept` p_sex_at_birth_concept 
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id",
  .sep = ""
)

query_helper(
  demo_sql_query,
  glue("demo_{config$version}"),
  project_path = project_directory,
  version = config$version
)

## drug data -------------------------------------------------------------------
cli_progress_step("querying drug data")

glp1_concepts <- unique(rbindlist(lapply(
  list.files("tables/codes/glp1ra", full.names = TRUE),
  fread
))[use == 1, concept_id])

sglt2_concepts <- unique(rbindlist(lapply(
  list.files("tables/codes/sglt2i", full.names = TRUE),
  fread
))[use == 1, concept_id])

dpp4_concepts <- unique(rbindlist(lapply(
  list.files("tables/codes/dpp4i", full.names = TRUE),
  fread
))[use == 1, concept_id])

all_concepts <- unique(c(glp1_concepts, sglt2_concepts, dpp4_concepts))

drug_sql_query <- glue::glue(
  "
    SELECT
        d_exposure.person_id,
        d_exposure.drug_concept_id,
        d_standard_concept.concept_name as standard_concept_name,
        d_standard_concept.concept_code as standard_concept_code,
        d_standard_concept.vocabulary_id as standard_vocabulary,
        d_exposure.drug_exposure_start_datetime,
        d_exposure.drug_exposure_end_datetime,
        d_exposure.verbatim_end_date,
        d_exposure.drug_type_concept_id,
        d_type.concept_name as drug_type_concept_name,
        d_exposure.stop_reason,
        d_exposure.refills,
        d_exposure.quantity,
        d_exposure.days_supply,
        d_exposure.sig,
        d_exposure.route_concept_id,
        d_route.concept_name as route_concept_name,
        d_exposure.lot_number,
        d_exposure.visit_occurrence_id,
        d_visit.concept_name as visit_occurrence_concept_name,
        d_exposure.drug_source_value,
        d_exposure.drug_source_concept_id,
        d_source_concept.concept_name as source_concept_name,
        d_source_concept.concept_code as source_concept_code,
        d_source_concept.vocabulary_id as source_vocabulary,
        d_exposure.route_source_value,
        d_exposure.dose_unit_source_value 
    FROM
        ( SELECT
            * 
        FROM
            `drug_exposure` d_exposure 
        WHERE
            (
                drug_concept_id IN (SELECT
                    DISTINCT ca.descendant_id 
                FROM
                    `cb_criteria_ancestor` ca 
                JOIN
                    (SELECT
                        DISTINCT c.concept_id       
                    FROM
                        `cb_criteria` c       
                    JOIN
                        (SELECT
                            CAST(cr.id as string) AS id             
                        FROM
                            `cb_criteria` cr             
                        WHERE
                            concept_id IN ({paste0(all_concepts, collapse = ', ')})             
                            AND full_text LIKE '%_rank1]%'       ) a 
                            ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                            OR c.path LIKE CONCAT('%.', a.id) 
                            OR c.path LIKE CONCAT(a.id, '.%') 
                            OR c.path = a.id) 
                    WHERE
                        is_standard = 1 
                        AND is_selectable = 1) b 
                        ON (ca.ancestor_id = b.concept_id)))) d_exposure 
        LEFT JOIN
            `concept` d_standard_concept 
                ON d_exposure.drug_concept_id = d_standard_concept.concept_id 
        LEFT JOIN
            `concept` d_type 
                ON d_exposure.drug_type_concept_id = d_type.concept_id 
        LEFT JOIN
            `concept` d_route 
                ON d_exposure.route_concept_id = d_route.concept_id 
        LEFT JOIN
            `visit_occurrence` v 
                ON d_exposure.visit_occurrence_id = v.visit_occurrence_id 
        LEFT JOIN
            `concept` d_visit 
                ON v.visit_concept_id = d_visit.concept_id 
        LEFT JOIN
            `concept` d_source_concept 
                ON d_exposure.drug_source_concept_id = d_source_concept.concept_id",
  .sep = ""
)

query_helper(
  drug_sql_query,
  glue("drug_{config$version}"),
  project_path = project_directory,
  version = config$version
)

glp1_descendants <- pull_descendant_ids(glp1_concepts)[,
  drug_class := "GLP-1 RA"
]
sglt2_descendants <- pull_descendant_ids(sglt2_concepts)[,
  drug_class := "SGLT2i"
]
dpp4_descendants <- pull_descendant_ids(dpp4_concepts)[, drug_class := "DPP4i"]

all_descendants <- rbindlist(list(
  glp1_descendants,
  sglt2_descendants,
  dpp4_descendants
))

fwrite(
  x = all_descendants,
  file = file.path(
    project_directory,
    "data",
    "raw",
    glue("drug_concept_ids_{config$version}.tsv")
  ),
  sep = "\t"
)

## metformin and insulins data -------------------------------------------------
cli_progress_step("querying metformin and insulins data")
metformin_concepts <- unique(fread("tables/codes/metformin_codes.txt")[
  use == 1,
  concept_id
])
insulins_concepts <- unique(fread("tables/codes/insulin_codes.txt")[
  use == 1,
  concept_id
])

mni_concepts <- unique(c(metformin_concepts, insulins_concepts))

metformin_and_insulins_sql_query <- glue::glue(
  "
    SELECT
        d_exposure.person_id,
        d_exposure.drug_concept_id,
        d_standard_concept.concept_name as standard_concept_name,
        d_standard_concept.concept_code as standard_concept_code,
        d_standard_concept.vocabulary_id as standard_vocabulary,
        d_exposure.drug_exposure_start_datetime,
        d_exposure.drug_exposure_end_datetime,
        d_exposure.verbatim_end_date,
        d_exposure.drug_type_concept_id,
        d_type.concept_name as drug_type_concept_name,
        d_exposure.stop_reason,
        d_exposure.refills,
        d_exposure.quantity,
        d_exposure.days_supply,
        d_exposure.sig,
        d_exposure.route_concept_id,
        d_route.concept_name as route_concept_name,
        d_exposure.lot_number,
        d_exposure.visit_occurrence_id,
        d_visit.concept_name as visit_occurrence_concept_name,
        d_exposure.drug_source_value,
        d_exposure.drug_source_concept_id,
        d_source_concept.concept_name as source_concept_name,
        d_source_concept.concept_code as source_concept_code,
        d_source_concept.vocabulary_id as source_vocabulary,
        d_exposure.route_source_value,
        d_exposure.dose_unit_source_value 
    FROM
        ( SELECT
            * 
        FROM
            `drug_exposure` d_exposure 
        WHERE
            (
                drug_concept_id IN (SELECT
                    DISTINCT ca.descendant_id 
                FROM
                    `cb_criteria_ancestor` ca 
                JOIN
                    (SELECT
                        DISTINCT c.concept_id       
                    FROM
                        `cb_criteria` c       
                    JOIN
                        (SELECT
                            CAST(cr.id as string) AS id             
                        FROM
                            `cb_criteria` cr             
                        WHERE
                            concept_id IN ({paste0(mni_concepts, collapse = ', ')})             
                            AND full_text LIKE '%_rank1]%'       ) a 
                            ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                            OR c.path LIKE CONCAT('%.', a.id) 
                            OR c.path LIKE CONCAT(a.id, '.%') 
                            OR c.path = a.id) 
                    WHERE
                        is_standard = 1 
                        AND is_selectable = 1) b 
                        ON (ca.ancestor_id = b.concept_id)))) d_exposure 
        LEFT JOIN
            `concept` d_standard_concept 
                ON d_exposure.drug_concept_id = d_standard_concept.concept_id 
        LEFT JOIN
            `concept` d_type 
                ON d_exposure.drug_type_concept_id = d_type.concept_id 
        LEFT JOIN
            `concept` d_route 
                ON d_exposure.route_concept_id = d_route.concept_id 
        LEFT JOIN
            `visit_occurrence` v 
                ON d_exposure.visit_occurrence_id = v.visit_occurrence_id 
        LEFT JOIN
            `concept` d_visit 
                ON v.visit_concept_id = d_visit.concept_id 
        LEFT JOIN
            `concept` d_source_concept 
                ON d_exposure.drug_source_concept_id = d_source_concept.concept_id",
  .sep = ""
)

query_helper(
  metformin_and_insulins_sql_query,
  glue("metformin_and_insulins_{config$version}"),
  project_path = project_directory,
  version = config$version
)

metformin_descendants <- pull_descendant_ids(metformin_concepts)[,
  drug_class := "Metformin"
]
insulins_descendants <- pull_descendant_ids(insulins_concepts)[,
  drug_class := "Insulins"
]
mni_descendants <- rbindlist(
  list(metformin_descendants, insulins_descendants),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  x = mni_descendants,
  file = file.path(
    project_directory,
    "data",
    "raw",
    glue("metformin_and_insulins_ids_{config$version}.tsv")
  ),
  sep = "\t"
)

## EHR (ICD code) data ---------------------------------------------------------
cli_progress_step("querying ICD data")
icd_sql_query <- glue::glue(
  "
    SELECT
        c_occurrence.person_id,
        c_occurrence.condition_concept_id,
        c_standard_concept.concept_name as standard_concept_name,
        c_standard_concept.concept_code as standard_concept_code,
        c_standard_concept.vocabulary_id as standard_vocabulary,
        c_occurrence.condition_start_datetime,
        c_occurrence.condition_end_datetime,
        c_occurrence.condition_type_concept_id,
        c_type.concept_name as condition_type_concept_name,
        c_occurrence.stop_reason,
        c_occurrence.visit_occurrence_id,
        visit.concept_name as visit_occurrence_concept_name,
        c_occurrence.condition_source_value,
        c_occurrence.condition_source_concept_id,
        c_source_concept.concept_name as source_concept_name,
        c_source_concept.concept_code as source_concept_code,
        c_source_concept.vocabulary_id as source_vocabulary,
        c_occurrence.condition_status_source_value,
        c_occurrence.condition_status_concept_id,
        c_status.concept_name as condition_status_concept_name
    FROM
        `condition_occurrence` c_occurrence
    LEFT JOIN
        `concept` c_standard_concept
            ON c_occurrence.condition_concept_id = c_standard_concept.concept_id
    LEFT JOIN
        `concept` c_type
            ON c_occurrence.condition_type_concept_id = c_type.concept_id
    LEFT JOIN
        `visit_occurrence` v
            ON c_occurrence.visit_occurrence_id = v.visit_occurrence_id
    LEFT JOIN
        `concept` visit
            ON v.visit_concept_id = visit.concept_id
    LEFT JOIN
        `concept` c_source_concept
            ON c_occurrence.condition_source_concept_id = c_source_concept.concept_id
    LEFT JOIN
        `concept` c_status
            ON c_occurrence.condition_status_concept_id = c_status.concept_id
    WHERE c_source_concept.vocabulary_id IN ('ICD9CM', 'ICD10CM')",
  .sep = ""
)

query_helper(
  icd_sql_query,
  glue("icd_{config$version}"),
  project_path = project_directory,
  version = config$version
)

## bmi data --------------------------------------------------------------------
cli_progress_step("querying BMI data")
bmi_concept <- fread("tables/omop_mappings.csv")[
  vocabulary_id == "LOINC" & group == "measurement",
  unique(concept_id)
]

bmi_sql_query <- glue::glue(
  "
    SELECT
        measurement.person_id,
        measurement.measurement_concept_id,
        m_standard_concept.concept_name as standard_concept_name,
        m_standard_concept.concept_code as standard_concept_code,
        m_standard_concept.vocabulary_id as standard_vocabulary,
        measurement.measurement_datetime,
        measurement.measurement_type_concept_id,
        m_type.concept_name as measurement_type_concept_name,
        measurement.operator_concept_id,
        m_operator.concept_name as operator_concept_name,
        measurement.value_as_number,
        measurement.value_as_concept_id,
        m_value.concept_name as value_as_concept_name,
        measurement.unit_concept_id,
        m_unit.concept_name as unit_concept_name,
        measurement.range_low,
        measurement.range_high,
        measurement.visit_occurrence_id,
        m_visit.concept_name as visit_occurrence_concept_name,
        measurement.measurement_source_value,
        measurement.measurement_source_concept_id,
        m_source_concept.concept_name as source_concept_name,
        m_source_concept.concept_code as source_concept_code,
        m_source_concept.vocabulary_id as source_vocabulary,
        measurement.unit_source_value,
        measurement.value_source_value 
    FROM
        ( SELECT
            * 
        FROM
            `measurement` measurement 
        WHERE
            (
                measurement_concept_id IN (SELECT
                    DISTINCT c.concept_id 
                FROM
                    `cb_criteria` c 
                JOIN
                    (SELECT
                        CAST(cr.id as string) AS id       
                    FROM
                        `cb_criteria` cr       
                    WHERE
                        concept_id IN ({paste0(bmi_concept, collapse = ', ')})       
                        AND full_text LIKE '%_rank1]%'      ) a 
                        ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                        OR c.path LIKE CONCAT('%.', a.id) 
                        OR c.path LIKE CONCAT(a.id, '.%') 
                        OR c.path = a.id) 
                WHERE
                    is_standard = 1 
                    AND is_selectable = 1)
            )) measurement 
    LEFT JOIN
        `concept` m_standard_concept 
            ON measurement.measurement_concept_id = m_standard_concept.concept_id 
    LEFT JOIN
        `concept` m_type 
            ON measurement.measurement_type_concept_id = m_type.concept_id 
    LEFT JOIN
        `concept` m_operator 
            ON measurement.operator_concept_id = m_operator.concept_id 
    LEFT JOIN
        `concept` m_value 
            ON measurement.value_as_concept_id = m_value.concept_id 
    LEFT JOIN
        `concept` m_unit 
            ON measurement.unit_concept_id = m_unit.concept_id 
    LEFT JOIN
        `visit_occurrence` v 
            ON measurement.visit_occurrence_id = v.visit_occurrence_id 
    LEFT JOIN
        `concept` m_visit 
            ON v.visit_concept_id = m_visit.concept_id 
    LEFT JOIN
        `concept` m_source_concept 
            ON measurement.measurement_source_concept_id = m_source_concept.concept_id",
  .sep = ""
)

query_helper(
  bmi_sql_query,
  glue("bmi_{config$version}"),
  project_path = project_directory,
  version = config$version
)

## ses data --------------------------------------------------------------------
cli_progress_step("querying SES data")
ses_sql_query <- glue::glue(
  "
    SELECT
        observation.person_id,
        observation.observation_datetime,
        zip_code.zip3_as_string as zip_code,
        zip_code.fraction_assisted_income as assisted_income,
        zip_code.fraction_high_school_edu as high_school_education,
        zip_code.median_income,
        zip_code.fraction_no_health_ins as no_health_insurance,
        zip_code.fraction_poverty as poverty,
        zip_code.fraction_vacant_housing as vacant_housing,
        zip_code.deprivation_index,
        zip_code.acs as american_community_survey_year
    FROM
        `zip3_ses_map` zip_code
    JOIN
        `observation` observation
            ON CAST(SUBSTR(observation.value_as_string,
        0,
        STRPOS(observation.value_as_string,
        '*') - 1) AS INT64) = zip_code.zip3
        AND observation_source_concept_id = 1585250
        AND observation.value_as_string NOT LIKE 'Res%'",
  .sep = ""
)

query_helper(
  ses_sql_query,
  glue("ses_{config$version}"),
  project_path = project_directory,
  version = config$version
)

## death data ------------------------------------------------------------------
cli_progress_step("querying death data")
death_sql_query <- glue::glue(
  "
    SELECT DISTINCT
        person_id,
        death_date,
        death_datetime,
        death_type_concept_id,
        cause_concept_id,
        cause_source_value,
        cause_source_concept_id
    FROM
        `death`",
  .sep = ""
)

query_helper(
  death_sql_query,
  glue("death_{config$version}"),
  project_path = project_directory,
  version = config$version
)

## lab data --------------------------------------------------------------------
cli_progress_step("querying lab data")
lab_concepts <- fread("tables/omop_mappings.csv")[
  vocabulary_id == "LOINC" & group == "lab",
  unique(concept_id)
]

lab_sql_query <- glue::glue(
  "
    SELECT
        measurement.person_id,
        measurement.measurement_concept_id,
        m_standard_concept.concept_name as standard_concept_name,
        m_standard_concept.concept_code as standard_concept_code,
        m_standard_concept.vocabulary_id as standard_vocabulary,
        measurement.measurement_datetime,
        measurement.measurement_type_concept_id,
        m_type.concept_name as measurement_type_concept_name,
        measurement.operator_concept_id,
        m_operator.concept_name as operator_concept_name,
        measurement.value_as_number,
        measurement.value_as_concept_id,
        m_value.concept_name as value_as_concept_name,
        measurement.unit_concept_id,
        m_unit.concept_name as unit_concept_name,
        measurement.range_low,
        measurement.range_high,
        measurement.visit_occurrence_id,
        m_visit.concept_name as visit_occurrence_concept_name,
        measurement.measurement_source_value,
        measurement.measurement_source_concept_id,
        m_source_concept.concept_name as source_concept_name,
        m_source_concept.concept_code as source_concept_code,
        m_source_concept.vocabulary_id as source_vocabulary,
        measurement.unit_source_value,
        measurement.value_source_value 
    FROM
        ( SELECT
            * 
        FROM
            `measurement` measurement 
        WHERE
            (
                measurement_concept_id IN (SELECT
                    DISTINCT c.concept_id 
                FROM
                    `cb_criteria` c 
                JOIN
                    (SELECT
                        CAST(cr.id as string) AS id       
                    FROM
                        `cb_criteria` cr       
                    WHERE
                        concept_id IN ({paste0(lab_concepts, collapse = ', ')})       
                        AND full_text LIKE '%_rank1]%'      ) a 
                        ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                        OR c.path LIKE CONCAT('%.', a.id) 
                        OR c.path LIKE CONCAT(a.id, '.%') 
                        OR c.path = a.id) 
                WHERE
                    is_standard = 1 
                    AND is_selectable = 1)
            )) measurement 
    LEFT JOIN
        `concept` m_standard_concept 
            ON measurement.measurement_concept_id = m_standard_concept.concept_id 
    LEFT JOIN
        `concept` m_type 
            ON measurement.measurement_type_concept_id = m_type.concept_id 
    LEFT JOIN
        `concept` m_operator 
            ON measurement.operator_concept_id = m_operator.concept_id 
    LEFT JOIN
        `concept` m_value 
            ON measurement.value_as_concept_id = m_value.concept_id 
    LEFT JOIN
        `concept` m_unit 
            ON measurement.unit_concept_id = m_unit.concept_id 
    LEFT JOIN
        `visit_occurrence` v 
            ON measurement.visit_occurrence_id = v.visit_occurrence_id 
    LEFT JOIN
        `concept` m_visit 
            ON v.visit_concept_id = m_visit.concept_id 
    LEFT JOIN
        `concept` m_source_concept 
            ON measurement.measurement_source_concept_id = m_source_concept.concept_id",
  .sep = ""
)

query_helper(
  lab_sql_query,
  glue("lab_{config$version}"),
  project_path = project_directory,
  version = config$version
)

## visit data ------------------------------------------------------------------
cli_progress_step("querying visit data")
visit_sql_query <- glue::glue(
  "
    SELECT
        visit.PERSON_ID,
        visit.visit_concept_id,
        v_standard_concept.concept_name as standard_concept_name,
        v_standard_concept.concept_code as standard_concept_code,
        v_standard_concept.vocabulary_id as standard_vocabulary,
        visit.visit_start_datetime,
        visit.visit_end_datetime,
        visit.visit_type_concept_id,
        v_type.concept_name as visit_type_concept_name,
        visit.visit_source_value,
        visit.visit_source_concept_id,
        v_source_concept.concept_name as source_concept_name,
        v_source_concept.concept_code as source_concept_code,
        v_source_concept.vocabulary_id as source_vocabulary
    FROM
        `visit_occurrence` visit 
    LEFT JOIN
        `concept` v_standard_concept 
            ON visit.visit_concept_id = v_standard_concept.concept_id 
    LEFT JOIN
        `concept` v_type 
            ON visit.visit_type_concept_id = v_type.concept_id 
    LEFT JOIN
        `concept` v_source_concept 
            ON visit.visit_source_concept_id = v_source_concept.concept_id 
    LEFT JOIN
        `concept` v_admitting_source_concept 
            ON visit.admitting_source_concept_id = v_admitting_source_concept.concept_id 
    LEFT JOIN
        `concept` v_discharge 
            ON visit.discharge_to_concept_id = v_discharge.concept_id 
    ",
  .sep = ""
)

query_helper(
  visit_sql_query,
  glue("visit_{config$version}"),
  project_path = project_directory,
  version = config$version
)

cli_progress_done()

cli_alert_success("Data successfully queried and saved 🥳")
