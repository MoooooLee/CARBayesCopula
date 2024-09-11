#' Fitted Values for CARBayes Models
#'
#' This function returns the fitted values for objects of class \code{CARBayes}.
#'
#' @param object An object of class \code{CARBayes}.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return A matrix or list of fitted values, depending on the model type (univariate or multivariate).
#'
#' @export
fitted.CARBayes <- function(object,...)
{
  #### Return the fitted values
  return(object$fitted_values)
}
