library(readxl)
library(janitor)
library(tidyverse)

meta = sample_sheet = readxl::read_excel("data/SY 0825I-03_SparN_Recollect_Sample_Loading_Table_09042025 (CMY_FINAL)_withExperimental Explanations_NI.xlsm", sheet = 5)  |> 
    clean_names()

write_csv(meta, "data/metadata.csv")
