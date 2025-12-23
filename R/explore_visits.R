# libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
  library(arrow)
  library(cli)
  library(glue)
  library(data.table)
  library(PheWAS)
  library(dplyr)
  library(purrr)
  library(janitor)
})

options(datatable.print.class = TRUE)
options(datatable.print.trunc.cols = TRUE)

project_directory <- "/home/jupyter/workspaces/duplicateofglp1rapharmacogenomics"
data_directory <- file.path(project_directory, "data")
raw_directory <- file.path(data_directory, "raw")
processed_directory <- file.path(data_directory, "processed")
version <- "20251103_r2r"

dir.create(processed_directory, showWarnings = FALSE)

# DEMOGRAPHICS
demo <- open_dataset(file.path(raw_directory, "demo")) |>
  collect()
setDT(demo)
demo[, `:=`(
  person_id = as.character(person_id),
  dob = as.Date(date_of_birth),
  sex = fcase(
    sex_at_birth == "Male",
    "Male",
    sex_at_birth == "Female",
    "Female"
  )
)]
setkey(demo, person_id)
demo <- demo[!is.na(sex), ]

# ICD CODE
icd <- open_dataset(file.path(raw_directory, "icd")) |>
  mutate(person_id = as.character(person_id)) |>
  filter(person_id %in% demo[, unique(person_id)]) |>
  select(person_id, date = condition_start_datetime) |>
  collect()
setDT(icd)
icd[, `:=`(
  person_id = as.character(person_id),
  date = as.Date(date)
)]
setkey(icd, person_id)

icd_n <- unique(icd[, .(person_id, date)])[,
  .(n_unique_dx_days = .N),
  person_id
]

rm("icd")

# VISITS
demo_and_icd_ids <- intersect(demo[, unique(person_id)], icd_n[, unique(person_id)])

visits <- open_dataset(file.path(raw_directory, "visit")) |>
  select(person_id = PERSON_ID, date = visit_start_datetime, visit_type = standard_concept_name) |>
  mutate(person_id = as.character(person_id)) |>
  filter(person_id %in% demo_and_icd_ids) |>
  # filter(visit_concept_id %in% c(9203, 581477, 38004515)) |>
  unique() |>
  collect()
setDT(visits)
visits[, `:=`(
  person_id = as.character(person_id),
  date = as.Date(date)
)]
setkey(visits, person_id)

# VISIT TYPES
n_visit_type <- visits[, .N, visit_type]

visit_type_plot <- n_visit_type[N > 100] |>
  ggplot(aes(x = reorder(visit_type, N), y = N)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Number of visits by type",
    subtitle = paste0("N = ", format(visits[, uniqueN(person_id)], big.mark = ","), " participants; 100+ occurrences"),
    x = "",
    y = "Number of visits"
  ) +
  scale_x_discrete(labels = scales::label_wrap(width = 50)) +
  scale_y_log10(labels = scales::comma) +
  coord_flip() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0, color = "gray40", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )
ggsave(
  plot = visit_type_plot,
  filename = "plots/visit_type_plot.pdf",
  width = 8, height = 8 * 1.618, device = cairo_pdf
)

# identify visit types with 1000 occurrences
n_visit_type_1000 <- n_visit_type[N >= 1000, visit_type]

# count number of visit types per person
n_visit_type_person <- visits[, .N, .(person_id, visit_type)]

# limit to visits with at least 1000 occurrences
n_visit_type_person <- n_visit_type_person[visit_type %in% n_visit_type_1000, ]

# fix visit_type variable for transformation
vec.vt <- make_clean_names(n_visit_type_1000)
names(vec.vt) <- n_visit_type_1000
vec.vt <- factor(vec.vt)
n_visit_type_person[, visit_type := as.factor(visit_type)]
n_visit_type_person[, visit_type2 := vec.vt[match(visit_type, names(vec.vt))]]

n_visit_type_person_wide <- dcast(
  n_visit_type_person,
  person_id ~ visit_type2,
  value.var = "N",
  fill = 0
)

n_total_visits <- unique(visits[, .(person_id, date)])[, .(total = .N), person_id]

n_visit_type_person_wide <- merge(
  n_visit_type_person_wide,
  n_total_visits,
  by = "person_id",
  all = TRUE
)

these_visit_types <- names(n_visit_type_person_wide)
these_visit_types <- these_visit_types[!these_visit_types %in% c("person_id", "NA")]

n_visit_type_person_wide <- merge(
  n_visit_type_person_wide,
  icd_n,
  by = "person_id",
  all.x = TRUE
)


visit_cor <- data.table(
  visit_type = these_visit_types
)

visit_cor$cor <- sapply(visit_cor$visit_type, \(x) cor(n_visit_type_person_wide[[x]], n_visit_type_person_wide[["n_unique_dx_days"]]))


visit_cor_plot <- visit_cor[visit_type != "total", ] |>
  ggplot(aes(x = reorder(visit_type, cor), y = cor)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Correlation between visit type and unique Dx days",
    subtitle = paste0("Correlation for total visits: ", format(round(visit_cor[visit_type == "total", cor], 3), nsmall = 3)),
    x = "",
    y = "Correlation"
  ) +
  scale_x_discrete(labels = scales::label_wrap(width = 50)) +
  scale_y_continuous(labels = scales::comma) +
  coord_flip() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0, color = "gray40", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )
ggsave(
  plot = visit_cor_plot,
  filename = "plots/visit_cor_plot.pdf",
  width = 12, height = 8, device = cairo_pdf
)

# total visits scatterplot
total_visit_scatter <- n_visit_type_person_wide |>
  ggplot(aes(x = n_unique_dx_days, y = total)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "gam", se = TRUE, color = "blue") +
  labs(
    title = "Diagnosis days by number of visits",
    subtitle = paste0("Correlation: ", format(round(visit_cor[visit_type == "total", cor], 3), nsmall = 3)),
    x = "Unique Dx days",
    y = "Total number of visits"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0, color = "gray40", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave(
  plot = total_visit_scatter,
  filename = "plots/total_visit_scatter.png",
  width = 8, height = 6, units = "in", dpi = 320
)
