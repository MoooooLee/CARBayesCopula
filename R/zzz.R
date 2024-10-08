#' @useDynLib CARBayesCopula
#' @importFrom Rcpp sourceCpp
#' @importFrom stats cor cov dbinom glm model.frame model.matrix model.offset model.response na.pass quantile rbinom rnorm runif sd var
#' @importFrom utils setTxtProgressBar txtProgressBar
NULL

#' AEDC 2021 Data by SA3
#'
#' This dataset contains data from the 2021 Australian Early Development Census (AEDC) aggregated by Statistical Area Level 3 (SA3). As for modeling and confidentiality reasons, the data has been aggregated and modified. For areas that have less than 5 children, we will denote the data as 1 as it can be identified as a small number and support for modeling.
#'
#' @format A data frame with 1670 rows and 14 variables:
#' \describe{
#'   \item{SA3_NAME21}{A string value representing the name of the SA3 region under the 2021 Australian Statistical Geography Standard (ASGS2021).}
#'   \item{SA4_NAME21}{A string value representing the name of the SA4 region under the 2021 Australian Statistical Geography Standard (ASGS2021).}
#'   \item{GCC_NAME21}{A string value representing the name of the Greater Capital City under the 2021 Australian Statistical Geography Standard (ASGS2021).}
#'   \item{STE_NAME21}{A string value representing the name of the State or Territory under the 2021 Australian Statistical Geography Standard (ASGS2021).}
#'   \item{IRSD}{A numeric value representing the Index of Relative Socio-economic Disadvantage (IRSD) for the SA3 region, which is provided by the Australian Bureau of Statistics (ABS) and the value were calculated based on the 2021 Census data.}
#'   \item{IRSAD}{A numeric value representing the Index of Relative Socio-economic Advantage and Disadvantage (IRSAD) for the SA3 region, which is provided by the Australian Bureau of Statistics (ABS) and the value were calculated based on the 2021 Census data.}
#'   \item{IER}{A numeric value representing the Index of Education and Occupation (IER) for the SA3 region, which is provided by the Australian Bureau of Statistics (ABS) and the value were calculated based on the 2021 Census data.}
#'   \item{IEO}{A numeric value representing the Index of Economic Resources (IEO) for the SA3 region, which is provided by the Australian Bureau of Statistics (ABS) and the value were calculated based on the 2021 Census data.}
#'   \item{URP}{A numeric value representing the Usual Resident Population (URP) for the SA3 region, which is provided by the Australian Bureau of Statistics (ABS).}
#'   \item{domain}{A string value representing the domain of the AEDC data, which includes the following categories: "Physical Health and Wellbeing", "Social Competence", "Emotional Maturity", "Language and Cognitive Skills", and "Communication Skills and General Knowledge".}
#'   \item{vulnerable}{A numeric value representing the count number of children who are vulnerable in the domain in the SA3.}
#'   \item{N}{A numeric value representing the total number of children sampled in the survey in that domain in the SA3.}
#' }
#' @source Australian Early Development Census (AEDC), 2021, Australian Bureau of Statistics (ABS).
"aedc2021_sa3"


#' Neighborhood Matrix of SA3 Areas in Australia
#'
#' A matrix representing the neighboring relationships between Statistical Area Level 3 (SA3) regions.
#'
#' @format A square matrix with \(334\) rows and \(334\) columns, where \(334\) is the number of SA3 regions:
#' \describe{
#'   \item{Row/Column}{Are in the order of the SA3 regions that provided in the `aedc2021_sa3` dataset.}
#'   \item{Entries}{Binary values: 1 indicates that two regions are neighbors, 0 otherwise.}
#' }
#' @details The matrix is symmetric, and diagonal entries are typically 0 (regions are not neighbors with themselves).
#' @source Derived from spatial data for SA3 regions in Australia.
"aedc_W_sa3"





