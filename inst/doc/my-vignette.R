## -----------------------------------------------------------------------------
#| label: A00-install
#| include: FALSE
#| eval: FALSE
## devtools::install_github("MoooooLee/CARBayesCopula")
## 


## -----------------------------------------------------------------------------
#| label: A01-install
#| include: TRUE
#| eval: FALSE
## devtools::install_github("MoooooLee/CARBayesCopula")
## 


## -----------------------------------------------------------------------------
#| label: A02-library
#| include: TRUE
#| output: FALSE
library(CARBayesCopula)
library(tidyverse)
library(spdep)
library(copula)
library(flextable)
library(ggpubr)



## -----------------------------------------------------------------------------
#| label: B01-read_data
#| include: TRUE

# Load the dataset
data("aedc_sa3")
SA3_dataset_GLM <- aedc_sa3 %>%
  st_drop_geometry() %>%
  select(SA3_NAME21, IRSD,
         Social_Vulnerable_21, Language_Vulnerable_21,
         Social_N_21, Language_N_21) %>%
  rename(social_vulnerable = Social_Vulnerable_21,
         language_vulnerable = Language_Vulnerable_21,
         social_trial = Social_N_21,
         language_trial = Language_N_21) %>%
  drop_na(social_vulnerable, language_vulnerable) 

head(SA3_dataset_GLM) %>%
  flextable()


## -----------------------------------------------------------------------------
#| label: B02-res_pobsres
#| include: TRUE
res <- tibble(
  res_social = (cbind(social_vulnerable, social_trial - social_vulnerable) ~ IRSD) %>%
    glm(family = binomial, data = SA3_dataset_GLM) %>%
    residuals(type = "pearson"),
  res_language = (cbind(language_vulnerable, language_trial - language_vulnerable) ~ IRSD) %>%
    glm(family = binomial, data = SA3_dataset_GLM) %>%
    residuals(type = "pearson")
) 

pobsres <- tibble(
  pobsres_social = pobs(res$res_social),
  pobsres_language = pobs(res$res_language)
) 

res_pobsres <- bind_cols(res, pobsres)
res_pobsres_samples <- res_pobsres %>% slice_head(n = 5)

plot1 <- res_pobsres %>%
  ggplot(aes(x = res_social, y = res_language)) +
  geom_density_2d_filled(alpha = 1, bins = 9) +
  geom_point(size = 0.6) +
  geom_segment(data = res_pobsres_samples,
               aes(x = res_social, xend = 10,
                   y = res_language, yend = res_language),
               linewidth = 0.5, color = "red") +
  geom_segment(data = res_pobsres_samples,
               aes(x = res_social, xend = res_social,
                   y = -10, yend = res_language),
               linewidth = 0.5, color = "red") +
  geom_point(data = res_pobsres_samples,
             size = 1.5, shape = 15,
             color = "red") + 
  geom_rug(sides = "br") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = " ", y = "Residuals in Language") +
  theme_minimal() +
  theme(legend.position = "none")

plot2 <- res_pobsres %>%
  ggplot(aes(x = pobsres_language, y = res_language)) +
  geom_segment(data = res_pobsres_samples,
               aes(x = -1, xend = pobsres_language,
                   y = res_language, yend = res_language),
               linewidth = 0.5, color = "red") +
  geom_segment(data = res_pobsres_samples,
               aes(x = pobsres_language, xend = pobsres_language,
                   y = -10, yend = res_language),
               linewidth = 0.5, color = "red") +
  geom_point(size = 0.6) +
  geom_point(data = res_pobsres_samples,
             size = 1.5, shape = 15,
             color = "red") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-5, 5)) +
  labs(x = " ", y = " ") +
  theme_minimal() +
  theme(legend.position = "none")

plot3 <- res_pobsres %>%
  ggplot(aes(x = res_social, y = pobsres_social)) +
  geom_segment(data = res_pobsres_samples,
               aes(x = res_social, xend = 10,
                   y = pobsres_social, yend = pobsres_social),
               linewidth = 0.5, color = "red") +
  geom_segment(data = res_pobsres_samples,
               aes(x = res_social, xend = res_social,
                   y = pobsres_social, yend = 10),
               linewidth = 0.5, color = "red") +
  geom_point(size = 0.6) +
  geom_point(data = res_pobsres_samples,
             size = 1.5, shape = 15,
             color = "red") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 1)) +
  labs(x = "Residuals in Social", y = "Probability integral transformed\nresiduals in Social ") +
  theme_minimal() +
  theme(legend.position = "none")

plot4 <- res_pobsres %>%
  ggplot(aes(x = pobsres_language, y = pobsres_social)) +
  geom_density_2d_filled(alpha = 1, bins = 9) +
  geom_point(size = 0.6) + 
  geom_segment(data = res_pobsres_samples,
               aes(x = -10, xend = pobsres_language,
                   y = pobsres_social, yend = pobsres_social),
               linewidth = 0.5, color = "red") +
  geom_segment(data = res_pobsres_samples,
               aes(x = pobsres_language, xend = pobsres_language,
                   y = pobsres_social, yend = 10),
               linewidth = 0.5, color = "red") +
  geom_point(data = res_pobsres_samples,
             size = 1.5, shape = 15,
             color = "red") +
  geom_rug(sides = "tl") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_fill_brewer(palette = "Blues") +
  labs(x = "Probability integral transformd\nresiduals in Language", y = " ") +
  theme_minimal() +
  theme(legend.position = "none")

ggarrange(plot1, plot2, plot3, plot4,
          ncol = 2, nrow = 2) %>%
  annotate_figure(
    top = text_grob("Residuals and Probability Integral Transforms",
                    face = "bold", size = 14),
    bottom = text_grob("between Social and Langauage Vulnerable",
                       face = "bold", size = 12)
  )


## -----------------------------------------------------------------------------
#| warning: false
#| label: B03-copulas
#| include: TRUE

n <- 2000
normalCopula(0.7, dim = 2) %>%
  rCopula(n,.) %>%
  as_tibble() %>%
  mutate(Copula = "Gaussian (alpha=0.7)") %>%
  bind_rows(
    frankCopula(3, dim = 2) %>%
      rCopula(n,.) %>%
      as_tibble() %>%
      mutate(Copula = "Frank (alpha=3)"),
    claytonCopula(1.5, dim = 2) %>%
      rCopula(n,.) %>%
      as_tibble() %>%
      mutate(Copula = "Clayton (alpha=1.5)"),
    gumbelCopula(2, dim = 2) %>%
      rCopula(n,.) %>%
      as_tibble() %>%
      mutate(Copula = "Gumbel (alpha=2)"),
    joeCopula(2, dim = 2) %>%
      rCopula(n,.) %>%
      as_tibble() %>%
      mutate(Copula = "Joe (alpha=2)")
  ) %>%
  ggplot() +
  geom_density_2d_filled(aes(x = V1,
                             y = V2),
                         bins = 9) +
  coord_fixed(ratio = 1) +
  scale_fill_brewer(palette = "Blues") +
  facet_wrap(~Copula)  +
  labs(title = "Density plot of copulas") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis labels
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(), # Hide x-axis ticks
    axis.text.y = element_blank(),  # Hide y-axis labels
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )


## -----------------------------------------------------------------------------
#| label: C01-BHM
#| include: TRUE

options(width = 120) # Set the output width to a larger value
data("aedc_W_sa3") # Load the neighborhood matrix

input_formula <- cbind(Social_Vulnerable_21, Language_Vulnerable_21) ~ IRSD
input_trials <- cbind(aedc_sa3$Social_N_21, aedc_sa3$Language_N_21)

set.seed(1234) # Set the seed for reproducibility
model_fit_clayton <- fit_CARleroux_copula(formula = input_formula,
                                         trials = input_trials,
                                         W = aedc_W_sa3,
                                         data = aedc_sa3,
                                         copula = "clayton",
                                         burnin = 200,
                                         n_sample = 300,
                                         thin = 5,
                                         verbose = FALSE)
model_fit_invgumbel <- fit_CARleroux_copula(formula = input_formula,
                                          trials = input_trials,
                                          W = aedc_W_sa3,
                                          data = aedc_sa3,
                                          copula = "inversegumbel",
                                         burnin = 200,
                                         n_sample = 300,
                                         thin = 5,
                                         verbose = FALSE)


## -----------------------------------------------------------------------------
#| label: C02-summary
#| include: TRUE
summary(model_fit_clayton)
summary(model_fit_invgumbel)


## -----------------------------------------------------------------------------
#| label: C02-load_rds
#| include: FALSE
#| eval: TRUE

main_table <- readRDS("RDS/A05_main_table.RDS")
perf_table_SA3 <- readRDS("RDS/A08_perf_table_SA3.RDS")
perf_table_SA3byIRSD <- readRDS("RDS/A08_perf_table_SA3byIRSD.RDS")

table_ic <- main_table %>%
  select(var1, var2, copula, ic) %>%
  unnest(ic) %>% 
  pivot_wider(names_from = ic, values_from = value) %>%
  nest(data = -(var1:var2)) 


## -----------------------------------------------------------------------------
#| label: C02-comparison
#| echo: false
table_ic$data[[10]] %>%
    as_tibble() %>%
    flextable(.)  %>%
    colformat_double(digits = 2) %>%
    color(i = ~LOOIC == min(LOOIC),
          j = ~LOOIC,
          color = "red") %>%
    color(i = ~WAIC == min(WAIC),
          j = ~WAIC,
          color = "orange") %>%
    color(i = ~loglikelihood == min(loglikelihood),
          j = ~loglikelihood,
          color = "blue") %>%
    bold(i = ~ (LOOIC == min(LOOIC)|
                  WAIC == min(WAIC))) %>%
    add_footer_lines(str_c(table_ic$var2[10],
                           " - ",
                           table_ic$var1[10])) %>%
    autofit()


## -----------------------------------------------------------------------------
#| label: C03-visualization
#| include: TRUE

est_model <- model_fit_invgumbel$fitted_values
colnames(est_model) <- c("social_est_model", "language_est_model")

aedc_sa3_results <- aedc_sa3 %>%
  select(SA3_NAME21, SA4_NAME21, GCC_NAME21, STE_NAME21,
         IRSD, 
         Social_Vulnerable_21, Language_Vulnerable_21,
         Social_N_21, Language_N_21) %>%
  bind_cols(est_model) %>%
  mutate(social_est_dir = Social_Vulnerable_21 / Social_N_21,
         language_est_dir = Language_Vulnerable_21 / Language_N_21,
         social_est_model = social_est_model / Social_N_21,
         language_est_model = language_est_model / Language_N_21) %>%
  pivot_longer(cols = contains(c("est_dir","est_model")),
               names_to = "domain_method",
               values_to = "est") %>%
  group_by(domain_method) %>%
  mutate(est_cat = factor(ntile(est, 10)))



## -----------------------------------------------------------------------------
#| label: C03-maps_dir
#| include: TRUE

aedc_sa3_results %>%
  filter(domain_method == "social_est_dir") %>%
  ggplot() +
  geom_sf(aes(fill = est_cat),
          color = "black") +
  scale_fill_brewer(palette = "RdYlGn",
                    direction = -1,
                    na.value = "grey") + 
  labs(title = "Vulnerability prevalence in Social",
       subtitle = "Estimated by direct proportion",
       fill = "Vulnerability Levels (Deciles)") +
  theme_void() +
  theme(legend.position = "bottom")

aedc_sa3_results %>%
  filter(domain_method == "language_est_dir") %>%
  ggplot() +
  geom_sf(aes(fill = est_cat),
          color = "black") +
  scale_fill_brewer(palette = "RdYlGn",
                    direction = -1,
                    na.value = "grey") + 
  labs(title = "Vulnerability prevalence in Language",
       subtitle = "Estimated by direct proportion",
       fill = "Vulnerability Levels (Deciles)") +
  theme_void() +
  theme(legend.position = "bottom")



## -----------------------------------------------------------------------------
#| label: C03-maps_inversegumbel
#| include: TRUE
aedc_sa3_results %>%
  filter(domain_method == "social_est_model") %>%
  ggplot() +
  geom_sf(aes(fill = est_cat),
          color = "black") +
  scale_fill_brewer(palette = "RdYlGn",
                    direction = -1,
                    na.value = "grey") + 
  labs(title = "Vulnerability prevalence in Social",
       subtitle = "Estimated by model with inverse Gumbel copula",
       fill = "Vulnerability Levels (Deciles)") +
  theme_void() +
  theme(legend.position = "bottom")

aedc_sa3_results %>%
  filter(domain_method == "language_est_model") %>%
  ggplot() +
  geom_sf(aes(fill = est_cat),
          color = "black") +
  scale_fill_brewer(palette = "RdYlGn",
                    direction = -1,
                    na.value = "grey") + 
  labs(title = "Vulnerability prevalence in Language",
       subtitle = "Estimated by model with inverse Gumbel copula",
       fill = "Vulnerability Levels (Deciles)") +
  theme_void() +
  theme(legend.position = "bottom")


