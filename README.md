# CARBayesCopula

## Overview

`CARBayesCopula` is an R package designed to fit a bivariate Bayesian hierarchical model with copula, which is developed by Mu Li and available on GitHub. The package implements bivariate spatial generalized linear mixed models for area unit data, with inference in a Bayesian setting using Markov chain Monte Carlo (MCMC) simulation, which is built on the `Rcpp` package and `CARBayes` package (Eddelbuettel and François 2011; Lee 2013).

As this is the first published version, the response variables can only be binomial. The package is designed to be user-friendly and computationally efficient. The package is still under development and more features will be added in the future.

## Installation

You can only install the pre-released version of `CARBayesCopula` from GitHub:

```r
# Install devtools if necessary
# install.packages("devtools")

# Install CARBayesCopula from GitHub
devtools::install_github("MoooooLee/CARBayesCopula")
```

As the package build by `Rcpp`, it may take a while to install.

## Features

- **Bivariate Bayesian hierarchical modeling** for spatial areal unit data with binomial responses.
- **Copula-based dependence modeling** between two response variables, allowing for flexible tail dependencies and correlations.
- **Spatial structure** modeled by CAR prior, accommodating spatial dependence between areal units.
- Available **copula functions**: Clayton, Frank, Gumbel, Joe, Gaussian, and their inverse counterparts for flexible dependence structures.
- Efficient implementation using **MCMC simulation** for Bayesian inference.
- Built on **Rcpp** and **CARBayes** for computational efficiency.

For more comprehensive details on model formulation, parameter selection, and advanced features, please refer to the package vignette, which provides in-depth explanations and examples. It is strongly recommended to review the [**vignette**](https://mooooolee.github.io/CARBayesCopula/vignettes/my-vignette.html) for a better understanding of the package's functionality.

## License

This project is open-source and available under the GPL-3.0 License.

## Contact

For any inquiries or feedback, please contact Mu Li at [mu.li@anu.edu.au](mu.li@anu.edu.au).

## References

Eddelbuettel, Dirk, and Romain François. 2011. “Rcpp: Seamless r and c++ Integration.” Journal of Statistical Software 40: 1–18.

Lee, Duncan. 2013. “CARBayes: An r Package for Bayesian Spatial Modeling with Conditional Autoregressive Priors.” Journal of Statistical Software 55 (13): 1–24.
