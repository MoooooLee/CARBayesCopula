#' Fit a Bivariate spatial generalized linear mixed model with Leroux structured random effects with copula
#'
#' Fit a bivariate spatial generalized linear mixed model to areal unit data, where the response variable follows binomial distribution and the random effect modeled by copula. The linear predictor is modeled by known covariates and a vector of random effects. The latter account for both spatial and between variable association, via a marginal distribution combined by copula. Spatial correlation is captured by the conditional autoregressive (CAR) prior proposed by Leroux et al. (2000), and between variable correlation is captured by a predefined copula function with no fixed structure. This is a special type of multivariate conditional autoregressive (MCAR) model. Further details are given in the vignette accompanying this package. Independent (over space) random effects can be obtained by setting \code{rho = 0}, while the intrinsic MCAR model can be obtained by setting \code{rho = 1}. Inference is conducted in a Bayesian setting using Markov chain Monte Carlo (MCMC) simulation. Missing (\code{NA}) values are allowed in the response, and posterior predictive distributions are created for the missing values using data augmentation. These are saved in the "samples" argument in the output of the function and are denoted by "\code{Y}". For a full model specification see the vignette accompanying this package.
#'
#' @param formula A formula for the covariate part of the model using the syntax of the \code{lm()} function. Offsets can be included here using the \code{offset()} function.  The covariates should each be a \eqn{K \times 1} vector, where \eqn{K} is the number of spatial units. The response variable can contain missing values (\code{NA}).
#' @param copula The type of copula to use. Available options include:
#' \code{"gaussian", "clayton", "inverseclayton", "frank", "joe", "inversejoe", "gumbel", "inversegumbel"}.
#' Each copula type may have a different valid range for the copula parameter alpha, as follows:
#'
#' \describe{
#'   \item{\code{"gaussian"}}{alpha can be taken from the interval \eqn{(-1,1)}, with alpha = 0 indicating independence.}
#'   \item{\code{"clayton"}}{alpha can be taken from \eqn{(0,+\infty)}, with values close to 0 indicating independence.}
#'   \item{\code{"inverseclayton"}}{alpha can be taken from \eqn{(0,+\infty)}, representing the inverse of the Clayton copula.}
#'   \item{\code{"frank"}}{alpha can take any real value other than 0, with closing to zero indicating independence.}
#'   \item{\code{"joe"}}{alpha can be taken from \eqn{[1,+\infty)}, with values close to 1 indicating weak dependence.}
#'   \item{\code{"inversejoe"}}{alpha can be taken from \eqn{[1,+\infty)}, representing the inverse of the Joe copula.}
#'   \item{\code{"gumbel"}}{alpha can be taken from \eqn{[1,+\infty)}, with alpha = 1 indicating independence.}
#'   \item{\code{"inversegumbel"}}{alpha can be taken from \eqn{[1,+\infty)}, representing the inverse of the Gumbel copula.}
#' }
#' @param data An optional data frame containing the variables in the model. If not found in data, the variables are taken from the environment from which `fit_CARleroux_copula` is called.
#' @param trials A matrix representing the number of trials for the binomial model, with dimensions matching those of the response variable.
#' @param W A spatial weight matrix representing neighborhood structures. which should be a non-negative \eqn{K \times K} neighborhood matrix (where \eqn{K} is the number of spatial units). Typically a binary specification is used, where the \eqn{(j,k)}th element equals one if areas \eqn{(j, k)} are spatially close (e.g. share a common border) and is zero otherwise. The matrix can be non-binary, but each row must contain at least one non-zero entry.
#' @param burnin The number of burn-in iterations for the MCMC sampling (default is \code{5000}).
#' @param n_sample The total number of MCMC samples (default is \code{10000}).
#' @param thin The thinning interval for the MCMC samples (default is \code{5}).
#' @param prior_beta_mean A vector of prior means for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector of zeros.
#' @param prior_beta_var A vector of prior variances for the regression parameters beta (Gaussian priors are assumed). Defaults to a vector with values \code{100000}.
#' @param prior_tau2 The prior shape and scale in the form of c(shape, scale) for an Inverse-Gamma(shape, scale) prior for tau2. Defaults to \code{c(1, 0.01)}.
#' @param alpha The copula parameter that controls the dependence structure between two marginal distributions of the random effects. The range of \code{alpha} varies depending on the copula type:
#' \itemize{
#'   \item Gaussian: \eqn{\alpha \in (-1, 1)}
#'   \item Clayton: \eqn{\alpha \in (0, \infty)}
#'   \item Inverse Clayton: \eqn{\alpha \in (-\infty, 0)}
#'   \item Frank: \eqn{\alpha \in (-\infty, \infty)}
#'   \item Joe: \eqn{\alpha \in (1, \infty)}
#'   \item Inverse Joe: \eqn{\alpha \in (-\infty, -1)}
#'   \item Gumbel: \eqn{\alpha \in [1, \infty)}
#'   \item Inverse Gumbel: \eqn{\alpha \in (-\infty, 1]}
#' }
#' If \code{alpha} is \code{NULL}, it will be estimated during model fitting. If provided, the specified value will be used as a fixed parameter.
#' @param rho The value in the interval \eqn{[0, 1]} that the spatial dependence parameter rho is fixed at if it should not be estimated. If this argument is \code{NULL} then rho is estimated in the model.
#' @param verbose Logical; if \code{TRUE}, prints progress during the MCMC sampling (default is \code{TRUE}).

#' @return A list containing the model fit, including posterior samples, fitted values, residuals, and model summary.
#' \describe{
#'  \item{\code{summary.results} }{A summary table of the parameters.}
#'  \item{\code{samples} }{A list containing the MCMC samples, as well as the loglik of each samples from the model.}
#'  \item{\code{fitted.values} }{A matrix of fitted values for each area and responsevariable.}
#'  \item{\code{residuals} }{A list with 2 elements, where each element is a vector of a type of residuals. The types of residual are "response" (raw), and "pearson".}
#'  \item{\code{modelfit} }{Model fit criteria including the Deviance Information Criterion(DIC) and its corresponding estimated effective number of parameters (p.d), the Log Marginal Predictive Likelihood (LMPL), the Watanabe-Akaike Information Criterion (WAIC) and its corresponding estimated number of effective parameters (p.w), the Leave-One-Out Information Criterion (LOOIC), the Expected Log Predictive Density (ELPD), and the loglikelihood.}
#'  \item{\code{accept} }{The acceptance probabilities for the parameters.}
#'  \item{\code{localised.structure} }{NULL, for compatability with other models.}
#'  \item{\code{formula} }{The formula (as a text string) for the response, covariate and offset parts of the model}
#'  \item{\code{model} }{A text string describing the model fit.}
#'  \item{\code{X} }{The design matrix of covariates.}
#'}
#'@references
#' \describe{
#' Gelfand, A and Vounatsou, P (2003). Proper multivariate conditional autoregressive models for spatial data analysis, Biostatistics, 4, 11-25.
#'
#' Kavanagh, L., D. Lee, and G. Pryce (2016). Is Poverty Decentralising? Quantifying Uncertainty in the Decentralisation of Urban Poverty, Annals of the American Association of Geographers, 106, 1286-1298.
#'
#' Leroux B, Lei X, Breslow N (2000). "Estimation of Disease Rates in SmallAreas: A New Mixed Model for Spatial Dependence." In M Halloran, D Berry (eds.), \emph{Statistical Models in Epidemiology, the Environment and Clinical Trials},pp. 179-191. Springer-Verlag, New York.
#' }
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
