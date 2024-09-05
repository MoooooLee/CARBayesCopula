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
  as_tibble() %>%
  pivot_longer(cols = c(Social, Emotional, Language, Communication),
               names_to = "Domain",
               values_to = "Residual") %>%
  ggplot() +
  geom_density_2d_filled(aes(x = Health,
                             y = Residual)) +
  coord_fixed(ratio = 1) +
  facet_wrap(~Domain)  +
  labs(title = "Empirical Copulas of Residuals",
       subtitle = "between Health and each of other variables") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis labels
    axis.text.y = element_blank(),  # Hide y-axis labels
    axis.ticks.x = element_blank(), # Hide x-axis ticks
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

# 4 Copula Graph ----------------------------------------------------------
n <- 329
normalCopula(0.7, dim = 2) %>%
  rCopula(n,.) %>%
  as_tibble() %>%
  mutate(Copula = "Normal") %>%
  bind_rows(
    frankCopula(5, dim = 2) %>%
      rCopula(n,.) %>%
      as_tibble() %>%
      mutate(Copula = "Frank"),
    claytonCopula(5, dim = 2) %>%
      rCopula(n,.) %>%
      as_tibble() %>%
      mutate(Copula = "Clayton"),
    gumbelCopula(5, dim = 2) %>%
      rCopula(n,.) %>%
      as_tibble() %>%
      mutate(Copula = "Gumbel"),
    joeCopula(5, dim = 2) %>%
      rCopula(n,.) %>%
      as_tibble() %>%
      mutate(Copula = "Joe")
  ) %>%
  ggplot() +
  geom_density_2d_filled(aes(x = V1,
                             y = V2)) +
  coord_fixed(ratio = 1) +
  facet_wrap(~Copula)  +
  labs(title = "Copulas") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis labels
    axis.text.y = element_blank(),  # Hide y-axis labels
    axis.ticks.x = element_blank(), # Hide x-axis ticks
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )




