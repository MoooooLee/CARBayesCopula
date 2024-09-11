#' Summarize a CARBayes Model
#'
#' This function provides a summary of a fitted CARBayes model, including model details,
#' posterior quantities, Deviance Information Criterion (DIC), and Log Marginal Predictive Likelihood (LMPL).
#'
#' @param x A CARBayes model object, typically the result of fitting a model with CARBayes.
#' @param ... Additional arguments (currently not used).
#'
#' @return Invisibly returns the input `x`.
#'
#' @details The function prints the model formula, number of missing observations, and a summary of the posterior
#' results, including DIC, p.d, and LMPL. If the model includes a localized structure, it also prints the number of stepchanges
#' or clusters identified in the random effect surface.
#'
#' @export

summary.CARBayes <- function(x,...)
{
  #### Check for missingness
  if(length(x$samples$Y)==1)
  {
    n_miss <- 0
  }else
  {
    n_miss <- ncol(x$samples$Y)
  }

  if(class(x$localised_structure)=="list")
  {
    #### Print out the model fitted
    cat("\n#################\n")
    cat("#### Model fitted\n")
    cat("#################\n")
    cat(x$model)
    cat("Regression equation - ")
    print(x$formula)
    cat("Number of missing observations - ")
    cat(n_miss)
    cat("\n")

    #### Print out the results
    cat("\n############\n")
    cat("#### Results\n")
    cat("############\n")
    cat("Posterior quantities and DIC\n\n")
    print(x$summary_results[ ,-c(4,5)])
    cat("\nDIC = ", x$modelfit[1], "   ",
        "p.d = ", x$modelfit[2], "   ",
        "LMPL = ", round(x$modelfit[5],2), "   ",
        "LOOIC = ", round(x$modelfit[7],2),
        "\n")

    if(length(x$localised_structure[[2]])>1)
    {
      cat("\nThe number of stepchanges identified in the random effect surface\n")
      temp <- x$localised_structure[[1]][!is.na(x$localised_structure[[1]])]
      tab <- array(NA, c(1,2))
      tab[1, ] <- c(sum(temp)/2, length(temp)/2- sum(temp)/2)
      colnames(tab) <- c("no stepchange", "stepchange")
      print(tab)
    }else
    {}
  }else if(class(x$localised_structure)=="numeric")
  {
    #### Print out the model fitted
    cat("\n#################\n")
    cat("#### Model fitted\n")
    cat("#################\n")
    cat(x$model)
    cat("Regression equation - ")
    print(x$formula)
    cat("Number of missing observations - ")
    cat(n_miss)
    cat("\n")

    #### Print out the results
    cat("\n############\n")
    cat("#### Results\n")
    cat("############\n")
    cat("Posterior quantities and DIC\n\n")
    print(x$summary_results[ ,-c(4,5)])
    cat("\nDIC = ", x$modelfit[1], "   ",
        "p.d = ", x$modelfit[2], "   ",
        "LMPL = ", round(x$modelfit[5],2), "   ",
        "LOOIC = ", round(x$modelfit[7],2),
        "\n")
    cat("\nNumber of clusters with the number of data points in each one\n")
    print(table(paste("group", x$localised_structure, sep="")))

  }else
  {
    #### Print out the model fitted
    cat("\n#################\n")
    cat("#### Model fitted\n")
    cat("#################\n")
    cat(x$model)
    if(class(x$formula)=="formula")
    {
      cat("Regression equation - ")
      print(x$formula)
    }else
    {
      cat("Regression equation - ")
      print(x$formula[[1]])
      cat("Zero probability equation - ")
      print(x$formula[[2]])
    }

    cat("Number of missing observations - ")
    cat(n_miss)
    cat("\n")

    #### Print out the results
    cat("\n############\n")
    cat("#### Results\n")
    cat("############\n")
    cat("Posterior quantities and DIC\n\n")
    print(x$summary_results[ ,-c(4,5)])
    cat("\nDIC = ", x$modelfit[1], "   ",
        "p.d = ", x$modelfit[2], "   ",
        "LMPL = ", round(x$modelfit[5],2), "   ",
        "LOOIC = ", round(x$modelfit[7],2),
        "\n")
  }
  return(invisible(x))
}
