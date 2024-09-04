#### H03 Model Fit ####
#### Test Function Code ####

# 1 Data and Input --------------------------------------------------------

source("Test_Draft/H01_ReadData.R")

formula <-  cbind(Health_Vulnerable_21,
                  Social_Vulnerable_21) ~ IRSD
data <- SA3_Dataset
trials <- SA3_Trials
W <- SA3_W

# 3 Create list of inputs ------------------------------------------------

Input_copula <- c("gaussian",
                  "clayton", "inverseclayton",
                  "joe", "inversejoe",
                  "gumbel", "inversegumbel",
                  "frank")

# 4 Run the function -----------------------------------------------------
output <- list()
for (i in 1:length(Input_copula)) {
  output[[i]] <-  fit_CARleroux_copula(formula = formula,
                                       trials = trials,
                                       W = W,
                                       data = SA3_Dataset,
                                       copula = "gaussian",
                                       burnin = 10000,
                                       n_sample = 15000,
                                       thin = 5)
}


