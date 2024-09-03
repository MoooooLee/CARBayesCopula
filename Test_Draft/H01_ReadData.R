#### A 1 Model Test ####

# 1 Set Seed and WD -------------------------------------------------------

set.seed(1234)
WD <- getwd()

# 2 Install and Lib Packages ----------------------------------------------

library(tidyverse)
library(spdep)


# 3 Data and Functions Input ----------------------------------------------

data("SA3_Dataset")
data("SA3_W")


# 4 Dataset Modification --------------------------------------------------

# Remove Areas without auxilllary covariates data
# Substitute trial numbers from 0 to 1
SA3_Dataset %<>%
  mutate(across(contains("_N_"),
                ~replace_na(.,1))) %>%
  mutate(across(contains("_N_"),
                ~replace(., . ==0, 1)))

SA3_Trials <- cbind(SA3_Dataset$Health_N_21,
                    SA3_Dataset$Social_N_21)

# End of this section -----------------------------------------------------

