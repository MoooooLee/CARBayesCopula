#' Residuals for CARBayes Models
#'
#' This function computes and returns residuals for objects of class \code{CARBayes}.
#'
#' @param object An object of class \code{CARBayes}.
#' @param type The type of residuals to return. Options are \code{"pearson"} (default) and \code{"response"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The residuals of the specified type. For multivariate models, the function returns a list of residuals; for univariate models, it returns a matrix of residuals.
#'
#' @export

residuals.CARBayes <- function(object, type="pearson", ...)
{
  residuals <- object$residuals


  #### The multivariate models provides lists the univariate models provide matrices

  if(inherits(residuals, "list"))
  {
    #### Return one of two types of residuals
    if(type=="response")
    {
      return(residuals$response)
    }else if(type=="pearson")
    {
      return(residuals$pearson)
    }else
    {
      return("Error. That is not one of the allowable residual types.")
    }

  }else
  {
    #### Return one of two types of residuals
    if(type=="response")
    {
      return(residuals$response)
    }else if(type=="pearson")
    {
      return(residuals$pearson)
    }else
    {
      return("Error. That is not one of the allowable residual types.")
    }
  }
}
