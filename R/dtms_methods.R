#' @export
## Function for data summary
summary.dtms_data <- function(object,...) {
  dtms_data_summary(data=object,...)
}

#' @export
## Function for transition probabilities
summary.dtms_probs <- function(object,...) {
  dtms_probs_summary(probs=object,...)
}
