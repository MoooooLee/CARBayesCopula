#### A 1 Model Test ####

# 1 Set Seed and WD -------------------------------------------------------

set.seed(1234)
WD <- getwd()

# 2 Install and Lib Packages ----------------------------------------------

library(tidyverse)
library(spdep)


# 3 Data and Functions Input ----------------------------------------------

data("aedc_sa3")
data("aedc_W_sa3")


# 4 Dataset Modification --------------------------------------------------

# Remove Areas without auxilllary covariates data
# Substitute trial numbers from 0 to 1
SA3_Dataset %<>%
  mutate(across(contains("_N_"),
                ~replace_na(.,1))) %>%
  mutate(across(contains("_N_"),
                ~replace(., . ==0, 1)))

# End of this section -----------------------------------------------------

