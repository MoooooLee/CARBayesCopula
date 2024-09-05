#### H03 Model Fit ####
#### Test Function Code ####

# 1 Data and Input --------------------------------------------------------

source("Test_Draft/H01_ReadData.R")

formula <-  cbind(Health_Vulnerable_21,
                  Social_Vulnerable_21) ~ IRSD

SA3_Trials <- cbind(SA3_Dataset$Health_N_21,
                    SA3_Dataset$Social_N_21)

data <- SA3_Dataset
trials <- SA3_Trials
W <- SA3_W

# 3 Create list of inputs ------------------------------------------------

input_copula <- c("gaussian",
                  "clayton", "inverseclayton",
                  "joe", "inversejoe",
                  "gumbel", "inversegumbel",
                  "frank")
input_Y <- c("Health", "Social","Emotional", "Language", "Communication") %>%
  crossing(Var1 = ., Var2 = .) %>%
  filter(Var1 <= Var2) %>%
  mutate(formula = str_c(Var1, "_N_21 + ", Var2, "_N_21 ~ IRSD"))

# 4 Run the function -----------------------------------------------------

for (i in 1:length(Input_copula)) {
  fit_CARleroux_copula(formula = formula,
                       trials = trials,
                       W = W,
                       data = SA3_Dataset,
                       copula = "gaussian",
                       burnin = 100,
                       n_sample = 150,
                       thin = 5) %>%
    saveRDS(file = paste0("Test_Draft/H03_model_fit/model_fit_",
                          Input_copula[i], ".rds"))
}


