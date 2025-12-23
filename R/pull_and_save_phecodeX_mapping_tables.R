#downloaded and saved on December 9, 2025
library(data.table)

phecodeX_labels <- fread(
  "https://github.com/PheWAS/PhecodeX/raw/main/phecodeX_R_labels.csv",
  showProgress = FALSE
)
phecodeX_rollup_map <- fread(
  "https://github.com/PheWAS/PhecodeX/raw/main/phecodeX_R_rollup_map.csv",
  showProgress = FALSE
)
phecodeX_map <- fread(
  "https://github.com/PheWAS/PhecodeX/raw/main/phecodeX_R_map.csv",
  showProgress = FALSE
) ## if you are using ICD-10 (not CM), load phecodeX_R_map_ICD_10_WHO.csv instead
phecodeX_sex <- fread(
  "https://github.com/PheWAS/PhecodeX/raw/main/phecodeX_R_sex.csv",
  showProgress = FALSE
)

fwrite(
  x = phecodeX_labels,
  file = "tables/phecodeX_labels.csv"
)
fwrite(
  x = phecodeX_rollup_map,
  file = "tables/phecodeX_rollup_map.csv"
)
fwrite(
  x = phecodeX_map,
  file = "tables/phecodeX_map.csv"
)
fwrite(
  x = phecodeX_sex,
  file = "tables/phecodeX_sex.csv"
)
