#' Fit a Bivariate Copula CAR Bayesian Model
#'
#' This function fits a bivariate copula CAR Bayesian model with a BYM structure
#' for a Binomial distribution.
#' Currently, it supports only the Gaussian copula.
#'
#' @param formula A symbolic description of the model to be fitted.
#' @param copula The type of copula to use (default is "gaussian"). Currently, only "gaussian" is supported.
#' @param data An optional data frame containing the variables in the model.
#' @param trials The number of trials for the Binomial model.
#' @param W A spatial weight matrix.
#' @param burnin The number of burn-in iterations for the MCMC sampling (default is 5000).
#' @param n_sample The total number of MCMC samples (default is 10000).
#' @param thin The thinning interval for the MCMC samples (default is 5).
#' @param prior_beta_mean Prior mean for the regression coefficients (optional).
#' @param prior_beta_var Prior variance for the regression coefficients (optional).
#' @param prior_tau2 Prior for the variance parameter tau^2 (optional).
#' @param prior_sigma2 Prior for the variance parameter sigma^2 (optional).
#' @param tau2_AB Hyperparameters for the prior distribution of tau^2 (optional).
#' @param sigma2_AB Hyperparameters for the prior distribution of sigma^2 (optional).
#' @param verbose Logical; if TRUE, prints progress during the MCMC sampling (default is TRUE).
#' @return A model object containing the results of the fitted CAR BYM model with the specified copula.
#' @import coda loo
#' @export

fit_CARbym_copula  <- function(formula,
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
                               prior_sigma2 = NULL,
                               tau2_AB = NULL,
                               sigma2_AB = NULL,
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
    model <- fit_CARbym_copula_gaussian(formula = formula,
                                        data = data,
                                        trials = trials,
                                        W = W,
                                        burnin = burnin,
                                        n_sample = n_sample,
                                        thin = thin,
                                        prior_beta_mean = prior_beta_mean,
                                        prior_beta_var = prior_beta_var,
                                        prior_tau2 = prior_tau2,
                                        prior_sigma2 = prior_sigma2,
                                        tau2_AB = tau2_AB,
                                        sigma2_AB = sigma2_AB,
                                        verbose = verbose)
  }else {
    stop("the copula argument is not valid",
         call.=FALSE)
  }

  return(model)
}
