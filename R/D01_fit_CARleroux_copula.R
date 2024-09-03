#' Fit a Bivariate Copula CAR Bayesian Model
#'
#' This function fits a bivariate copula CAR Bayesian model with a Leroux structure
#' for a Binomial distribution.
#'
#' @param formula The model formula.
#' @param copula The type of copula to use (e.g., "gaussian", "clayton").
#' @param data The dataset to be used.
#' @param trials The number of trials for the Binomial model.
#' @param W Spatial weight matrix.
#' @param burnin Number of burn-in iterations for the MCMC.
#' @param n_sample Total number of MCMC samples.
#' @param thin Thinning interval for the MCMC.
#' @param prior_beta_mean Prior mean for the regression coefficients.
#' @param prior_beta_var Prior variance for the regression coefficients.
#' @param prior_tau2 Prior for the variance parameter.
#' @param alpha Copula parameter alpha.
#' @param rho Spatial autocorrelation parameter.
#' @param verbose Logical; if TRUE, prints progress during the MCMC sampling.
#' @return A model object containing the results of the fitted model.
#' @import coda loo
#' @export

fit_CARleroux_copula <- function(formula,
                                 copula = "gaussian",
                                 data = NULL,
                                 trials,
                                 W,
                                 burnin = 5000,
                                 n_sample = 10000,
                                 thin = 5,
                                 prior_beta_mean = NULL,
                                 prior_beta_var = NULL,
                                 prior_tau2 = NULL,
                                 alpha = NULL,
                                 rho = NULL,
                                 verbose = TRUE){
  if(is.null(copula)){
    stop("the copula argument is missing",
         call.=FALSE)
  }

  if(is.null(trials)){
    stop("a binomial model was specified but the trials arugment was not specified",
         call.=FALSE)
  }

  if(copula == "gaussian") {
    model <- fit_CARleroux_copula_gaussian(formula = formula,
                                        data = data,
                                        trials = trials,
                                        W = W,
                                        burnin = burnin,
                                        n_sample = n_sample,
                                        thin = thin,
                                        prior_beta_mean = prior_beta_mean,
                                        prior_beta_var = prior_beta_var,
                                        prior_tau2 = prior_tau2,
                                        alpha = alpha,
                                        rho = rho,
                                        verbose = verbose)
  }else if(copula == "frank") {
    model <- fit_CARleroux_copula_frank(formula = formula,
                                        data = data,
                                        trials = trials,
                                        W = W,
                                        burnin = burnin,
                                        n_sample = n_sample,
                                        thin = thin,
                                        prior_beta_mean = prior_beta_mean,
                                        prior_beta_var = prior_beta_var,
                                        prior_tau2 = prior_tau2,
                                        alpha = alpha,
                                        rho = rho,
                                        verbose = verbose)
  }else if(copula == "clayton") {
    model <- fit_CARleroux_copula_clayton(formula = formula,
                                       data = data,
                                       trials = trials,
                                       W = W,
                                       burnin = burnin,
                                       n_sample = n_sample,
                                       thin = thin,
                                       prior_beta_mean = prior_beta_mean,
                                       prior_beta_var = prior_beta_var,
                                       prior_tau2 = prior_tau2,
                                       alpha = alpha,
                                       rho = rho,
                                       verbose = verbose)
  }else if(copula == "inverseclayton") {
    model <- fit_CARleroux_copula_invclayton(formula = formula,
                                              data = data,
                                              trials = trials,
                                              W = W,
                                              burnin = burnin,
                                              n_sample = n_sample,
                                              thin = thin,
                                              prior_beta_mean = prior_beta_mean,
                                              prior_beta_var = prior_beta_var,
                                              prior_tau2 = prior_tau2,
                                              alpha = alpha,
                                              rho = rho,
                                              verbose = verbose)
  }else if(copula == "joe") {
    model <- fit_CARleroux_copula_joe(formula = formula,
                                   data = data,
                                   trials = trials,
                                   W = W,
                                   burnin = burnin,
                                   n_sample = n_sample,
                                   thin = thin,
                                   prior_beta_mean = prior_beta_mean,
                                   prior_beta_var = prior_beta_var,
                                   prior_tau2 = prior_tau2,
                                   alpha = alpha,
                                   rho = rho,
                                   verbose = verbose)
  }else if(copula == "inversejoe") {
    model <- fit_CARleroux_copula_invjoe(formula = formula,
                                          data = data,
                                          trials = trials,
                                          W = W,
                                          burnin = burnin,
                                          n_sample = n_sample,
                                          thin = thin,
                                          prior_beta_mean = prior_beta_mean,
                                          prior_beta_var = prior_beta_var,
                                          prior_tau2 = prior_tau2,
                                          alpha = alpha,
                                          rho = rho,
                                          verbose = verbose)
  }else if(copula == "gumbel") {
    model <- fit_CARleroux_copula_gumbel(formula = formula,
                                      data = data,
                                      trials = trials,
                                      W = W,
                                      burnin = burnin,
                                      n_sample = n_sample,
                                      thin = thin,
                                      prior_beta_mean = prior_beta_mean,
                                      prior_beta_var = prior_beta_var,
                                      prior_tau2 = prior_tau2,
                                      alpha = alpha,
                                      rho = rho,
                                      verbose = verbose)
  }else if(copula == "inversegumbel"){
    model <- fit_CARleroux_copula_invgumbel(formula = formula,
                                             data = data,
                                             trials = trials,
                                             W = W,
                                             burnin = burnin,
                                             n_sample = n_sample,
                                             thin = thin,
                                             prior_beta_mean = prior_beta_mean,
                                             prior_beta_var = prior_beta_var,
                                             prior_tau2 = prior_tau2,
                                             alpha = alpha,
                                             rho = rho,
                                             verbose = verbose)
  }else {
    stop("the copula argument is not valid",
         call.=FALSE)
  }

  return(model)
}
