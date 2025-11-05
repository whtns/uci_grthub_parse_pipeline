library(readxl)
library(janitor)
library(tidyverse)

meta <- read_delim("data/sample_list.txt", delim = " ", col_names = c("sample", "wells")) |> 
  dplyr::mutate(condition = case_when(
    str_detect(sample, "differentiated") ~ "differentiated",
    str_detect(sample, "ctrl") ~ "ctrl",
    TRUE ~ "Unknown"
  ))  |> 
  identity()

write_csv(meta, "data/metadata.csv")
