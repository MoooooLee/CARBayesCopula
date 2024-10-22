\name{CARBayesCopula-package}
\alias{CARBayesCopula-package}
\alias{CARBayesCopula}
\docType{package}
\title{
  Spatial Generalised Linear Mixed Copula Models for Areal Unit Data
}
\description{The CARBayesCopula package provides tools for fitting spatial generalized linear mixed models with conditional autoregressive (CAR) structures and copula-based bivariate models, currently supporting a binomial framework. It extends the framework of the CARBayes package (developed by Duncan Lee) to allow for the modeling of spatial dependencies in areal unit data through copula structures. This enables flexible Bayesian inference for spatial data analysis, accommodating spatial autocorrelation and joint modeling of two binomial response variables. The package is particularly useful in fields such as epidemiology, ecology, and social sciences, where spatial dependencies and multivariate outcomes need to be modeled together.
}
\details{
  The CARBayesCopula package is designed for fitting spatial generalized linear mixed models (GLMMs) with Conditional Autoregressive (CAR) priors and copula-based dependence structures. The primary focus of this package is to extend the univariate CAR models implemented in the CARBayes package (Duncan Lee) to accommodate bivariate spatial data using copula functions. The package currently supports binomial models, allowing joint modeling of two spatially dependent binomial outcomes.

  Key features include:
  \itemize{
    \item Flexible copula-based modeling of spatial dependencies in areal unit data.
    \item Supports several copula types, including Gaussian, Clayton, Gumbel, Joe, and Frank copulas.
    \item Ability to handle missing data and allow for Bayesian inference using MCMC.
    \item Extends the CARBayes framework to model bivariate binomial data in a spatial context.
  }

  This package is particularly useful in fields such as epidemiology, environmental science, social sciences, and any domain where spatial dependencies between areas and bivariate outcomes are essential to model. It enables flexible Bayesian inference and can be customized for other applications.
}

\author{
Mu Li

Maintainer: Mu Li <mu.li@anu.edu.au>
}

\references{
Besag, J. and  York, J and Mollie, A (1991). Bayesian image restoration with two
applications in spatial statistics. Annals of the Institute of Statistics and
Mathematics 43, 1-59.

Gelfand, A and Vounatsou, P (2003). Proper multivariate conditional autoregressive
models for spatial data analysis, Biostatistics, 4, 11-25.

Kavanagh, L., D. Lee, and G. Pryce (2016). Is Poverty Decentralising? Quantifying
Uncertainty in the Decentralisation of Urban Poverty, Annals of the American
Association of Geographers, 106, 1286-1298.

Lee, D. and Mitchell, R (2012). Boundary detection in disease mapping studies.
Biostatistics, 13, 415-426.

Lee, D and Sarran, C  (2015). Controlling for unmeasured confounding and spatial
misalignment in long-term air pollution and health studies, Environmetrics, 26,
477-487.

Leroux B, Lei X, Breslow N (2000). "Estimation of Disease Rates in SmallAreas: A
New Mixed Model for Spatial Dependence." In M Halloran, D Berry (eds.),
\emph{Statistical Models in Epidemiology, the Environment and Clinical Trials},
pp. 179-191. Springer-Verlag, New York.

Roberts, G and Rosenthal, J (1998). Optimal scaling of discrete approximations to
the Langevin diffusions, Journal of the Royal Statistical Society Series B 60,
255-268.

Sklar, A. (1973). Random variables, joint distribution functions, and copulas. Kybernetika, 9(6), 449-460.

Nelsen, R. B. (2006). An introduction to copulas. Springer.
}
\keyword{
  spatial, copula, Bayesian, GLMM, bivariate, CARBayes, binomial, epidemiology
}
\seealso{
  For more information and the latest updates, please visit the \href{https://github.com/MoooooLee/CARBayesCopula}{CARBayesCopula}.
}

