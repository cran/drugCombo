#' Example checkerboard design drug combination data
#'
#' A dataset containing data from 10 different in-vitro drug combination 
#' experiments for an antiviral treatment (compound 1) using a checkerboard 
#' design. The variables are as follows:
#'
#' @name checkerboardData
#' @docType data
#' @format A data frame with 9360 observations of 5 variables:
#' \itemize{
#'  \item d1: dose of the first compound
#'  \item d2: dose of the second compound
#'  \item effect: observed effect (normalized cell counts)
#'  \item plate: plate ID (1 to 3 within experiment)
#'  \item exp: experiment ID (1 to 10)
#' }
NULL

#' Example ray design drug combination data
#'
#' A dataset containing data from an in-vitro drug combination experiment in 
#' oncology using a ray design. The variables are as follows:
#'
#' @name rayData
#' @docType data
#' @format A data frame with 378 observations of 4 variables:
#' \itemize{
#'  \item d1: dose of the first compound
#'  \item d2: dose of the second compound
#'  \item effect: observed effect (radioactivity level)
#'  \item ray: a character vector with values 0, 0.2, 0.35, 0.5, 
#' 0.65, 0.8 and 1, corresponding to the mixture factors. Each ray has 9 dose 
#' combinations with 6 replicates.
#' }
NULL
