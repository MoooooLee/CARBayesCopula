% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/F01_summary.CARBayes.R
\name{summary.CARBayes}
\alias{summary.CARBayes}
\title{Summarize a CARBayes Model}
\usage{
\method{summary}{CARBayes}(object, ...)
}
\arguments{
\item{object}{A CARBayes model object, typically the result of fitting a model with CARBayes.}

\item{...}{Additional arguments (currently not used).}
}
\value{
Invisibly returns the input \code{object}.
}
\description{
This function provides a summary of a fitted CARBayes model, including model details,
posterior quantities, Deviance Information Criterion (DIC), and Log Marginal Predictive Likelihood (LMPL).
}
\details{
The function prints the model formula, number of missing observations, and a summary of the posterior
results, including DIC, p.d, and LMPL. If the model includes a localized structure, it also prints the number of stepchanges
or clusters identified in the random effect surface.
}
