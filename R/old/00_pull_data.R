# libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
    library(arrow)
    library(cli)
    library(glue)
    library(data.table)
    library(bigrquery)
    library(dplyr)
})

project_directory <- "/home/jupyter/workspaces/duplicateofglp1rapharmacogenomics/"
version           <- "20250129"

dir.create(file.path(project_directory, "data", "raw"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(project_directory, "results"), showWarnings = FALSE, recursive = TRUE)

# functions --------------------------------------------------------------------
query_helper <- function(query, data_name) {
    # generate the path to save the data
    query_path <- file.path(
        Sys.getenv("WORKSPACE_BUCKET"),
        "bq_exports",
        # Sys.getenv("OWNER_EMAIL"),
        version,
        glue("{data_name}_{version}"),
        glue("{data_name}_{version}_*.parquet")
    )
    
    # save the data to Cloud Storage
    bq_table_save(
        bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), query, billing = Sys.getenv("GOOGLE_PROJECT")),
        query_path,
        destination_format = "parquet",
        compression        = "zstd",
        compression_lelvel = 12
    )

    # create the directory to save the data
    data_directory <- file.path(project_directory, "data", "raw", data_name)
    dir.create(data_directory, recursive = TRUE, showWarnings = FALSE)

    # download the data from Cloud Storage to the project directory
    system2("gsutil", paste("-m cp ", query_path, data_directory))
}

# query data -------------------------------------------------------------------
# demographics data
cli_progress_step("querying demographic data")
demo_sql_query <- paste("
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
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id", sep="")

query_helper(demo_sql_query, "demo")

# drug data
cli_progress_step("querying drug data")
drug_sql_query <- paste("
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
                            concept_id IN (1123618, 1123627, 1580747, 1583722, 21600783, 40166035, 40170911, 40239216, 43013884, 43526465, 44506754, 44785829, 44816332, 45774435, 45774751, 779705, 793143, 793293)             
                            AND full_text LIKE '%_rank1]%'       ) a 
                            ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                            OR c.path LIKE CONCAT('%.', a.id) 
                            OR c.path LIKE CONCAT(a.id, '.%') 
                            OR c.path = a.id) 
                    WHERE
                        is_standard = 1 
                        AND is_selectable = 1) b 
                        ON (ca.ancestor_id = b.concept_id)))  
                    AND (d_exposure.PERSON_ID IN (SELECT
                        distinct person_id  
                FROM
                    `cb_search_person` cb_search_person  
                WHERE
                    cb_search_person.person_id IN (SELECT
                        criteria.person_id 
                    FROM
                        (SELECT
                            DISTINCT person_id, entry_date, concept_id 
                        FROM
                            `cb_search_all_events` 
                        WHERE
                            (concept_id IN(SELECT
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
                                        concept_id IN (21600783, 1123618, 1123627, 779705)             
                                        AND full_text LIKE '%_rank1]%'       ) a 
                                        ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                                        OR c.path LIKE CONCAT('%.', a.id) 
                                        OR c.path LIKE CONCAT(a.id, '.%') 
                                        OR c.path = a.id) 
                                WHERE
                                    is_standard = 1 
                                    AND is_selectable = 1) b 
                                    ON (ca.ancestor_id = b.concept_id)) 
                                AND is_standard = 1)) criteria ) ))) d_exposure 
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
                ON d_exposure.drug_source_concept_id = d_source_concept.concept_id", sep="")

query_helper(drug_sql_query, "drug")

# short-read whole genome sequencing data (srWGS)
cli_progress_step("querying short-read whole genome sequencing data")
srWGS_sql_query <- paste("
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
            ON person.sex_at_birth_concept_id = p_sex_at_birth_concept.concept_id  
    WHERE
        person.PERSON_ID IN (SELECT
            distinct person_id  
        FROM
            `cb_search_person` cb_search_person  
        WHERE
            cb_search_person.person_id IN (SELECT
                person_id 
            FROM
                `cb_search_person` p 
            WHERE
                has_whole_genome_variant = 1 ) 
            AND cb_search_person.person_id IN (SELECT
                criteria.person_id 
            FROM
                (SELECT
                    DISTINCT person_id, entry_date, concept_id 
                FROM
                    `cb_search_all_events` 
                WHERE
                    (concept_id IN(SELECT
                        DISTINCT c.concept_id 
                    FROM
                        `cb_criteria` c 
                    JOIN
                        (SELECT
                            CAST(cr.id as string) AS id       
                        FROM
                            `cb_criteria` cr       
                        WHERE
                            concept_id IN (4274025)       
                            AND full_text LIKE '%_rank1]%'      ) a 
                            ON (c.path LIKE CONCAT('%.', a.id, '.%') 
                            OR c.path LIKE CONCAT('%.', a.id) 
                            OR c.path LIKE CONCAT(a.id, '.%') 
                            OR c.path = a.id) 
                    WHERE
                        is_standard = 1 
                        AND is_selectable = 1) 
                    AND is_standard = 1 )) criteria ) )", sep="")

query_helper(srWGS_sql_query, "srWGS")

# EHR (ICD code) data
cli_progress_step("querying ICD data")
icd_sql_query <- paste("
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
    WHERE c_source_concept.vocabulary_id IN ('ICD9CM', 'ICD10CM')", sep = "")

query_helper(icd_sql_query, "icd")

# bmi data
cli_progress_step("querying BMI data")
bmi_sql_query <- paste("
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
                        concept_id IN (3038553)       
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
    LEFT JOIn
        `visit_occurrence` v 
            ON measurement.visit_occurrence_id = v.visit_occurrence_id 
    LEFT JOIN
        `concept` m_visit 
            ON v.visit_concept_id = m_visit.concept_id 
    LEFT JOIN
        `concept` m_source_concept 
            ON measurement.measurement_source_concept_id = m_source_concept.concept_id", sep="")

query_helper(bmi_sql_query, "bmi")

# ses data
cli_progress_step("querying SES data")
ses_sql_query <- paste("
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
        AND observation.value_as_string NOT LIKE 'Res%'", sep = "")

query_helper(ses_sql_query, "ses")

# death data 
cli_progress_step("querying death data")
death_sql_query <- paste("
    SELECT DISTINCT
        person_id,
        death_date,
        death_datetime,
        death_type_concept_id,
        cause_concept_id,
        cause_source_value,
        cause_source_concept_id
    FROM
        `death`", sep = "")

query_helper(death_sql_query, "death")

# lab data
cli_progress_step("querying lab data")
lab_sql_query <- paste("
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
                        concept_id IN (3028288, 3001802, 3006923, 3016771, 3013721, 3020460, 3007070, 3027114, 3016723, 3004501, 3004410, 46236952, 3012888, 3004249)       
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
    LEFT JOIn
        `visit_occurrence` v 
            ON measurement.visit_occurrence_id = v.visit_occurrence_id 
    LEFT JOIN
        `concept` m_visit 
            ON v.visit_concept_id = m_visit.concept_id 
    LEFT JOIN
        `concept` m_source_concept 
            ON measurement.measurement_source_concept_id = m_source_concept.concept_id", sep="")

query_helper(lab_sql_query, "lab")

cli_progress_done()

cli_alert_success("Data successfully queried and saved 🥳")
