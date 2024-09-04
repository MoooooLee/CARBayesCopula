#### H02 EDA ####

# 1 Set Seed and WD -------------------------------------------------------

set.seed(1234)
WD <- getwd()

# 2 Install and Lib Packages ----------------------------------------------

library(tidyverse)
library(spdep)
library(copula)
source("Test_Draft/H01_ReadData.R")


# 3 Residual Analysis -----------------------------------------------------

SA3_Dataset_GLM <- SA3_Dataset %>%
  drop_na(Health_Vulnerable_21,
          Social_Vulnerable_21,
          Emotional_Vulnerable_21,
          Language_Vulnerable_21,
          Communication_Vulnerable_21) %>%
  pivot_longer(cols = contains(c("Vulnerable","_N_","AtRisk","OnTrack")),
               names_to = "var",
               values_to = "Count") %>%
  separate(var, into = c("Domain","var","Year"), sep = "_") %>%
  filter(var %in% c("Vulnerable","N"),
         Year == "21") %>%
  pivot_wider(names_from = "var",
              values_from = "Count")

# 3.1 Fit GLM -------------------------------------------------------------

Residual <- matrix(NA, nrow = 329, ncol = 5) %>%
  as.data.frame()
U_Residual <- matrix(NA, nrow = 329, ncol = 5) %>%
  as.data.frame()
i <- 2
for (i in 1:5) {
  Domain_Name <- c("Health", "Social","Emotional","Language","Communication")[i]
  Residual[,i] <- SA3_Dataset_GLM %>%
    filter(Year == "21",
           Domain == Domain_Name) %>%
    glm(cbind(Vulnerable, N-Vulnerable) ~ IRSD,
        family = binomial,
        data = .) %>%
    residuals()
  U_Residual[,i] <- pobs(Residual[,i])
  colnames(Residual)[i] <- Domain_Name
  colnames(U_Residual)[i] <- Domain_Name
}

# 3.2 Plot Residuals ------------------------------------------------------

U_Residual %>%
  ggplot() +
  geom_density_2d_filled(aes(x = Health,
                             y = Social))



