#' simpledata: an artificial dataset with abstract trajectories
#'
#' An artificial dataset with abstract states and time scale. The
#' state space consists of two transient states (A,B) and one absorbing state (X).
#'
#' @format `simpledata`
#' A data frame with 10,529 rows and 3 columns:
#' \describe{
#'   \item{id}{Identifier of the units}
#'   \item{age}{Age}
#'   \item{state}{The state occupied by an unit at a given age}
#' }
"simpledata"

#' hrsdata: simulated working trajectories based on the Health and Retirement Study
#'
#' A simulated dataset of individuals' working trajectories during late working life and retirement
#' age. The state space consists of three transient states (Employed,Retired,Inactive)
#' and one absorbing state (Dead). The age range covers ages 50 to 99. The data is
#' simulated using transition probabilities published as part of Dudel & Myrskylä (2017).
#'
#' @format `hrsdata`
#' A data frame with 250,000 rows and 4 columns:
#' \describe{
#'   \item{ID}{Person identifier}
#'   \item{Gender}{Individuals' gender (0=men, 1=women)}
#'   \item{Age}{Age}
#'   \item{State}{The state occupied by an unit at a given age}
#' }
#' @source <https://doi.org/10.1007/s13524-017-0619-6>
"hrsdata"

