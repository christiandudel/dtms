#' @export
## Summary function for data in transition format
summary.dtms_data <- function(object,...) {
  dtms_data_summary(data=object,...)
}

#' @export
## Summary function for transition probabilities
summary.dtms_probs <- function(object,...) {
  dtms_probs_summary(probs=object,...)
}

#' @export
## Summary function for distributional results
summary.dtms_distr <- function(object,...) {
  dtms_distr_summary(distr=object,...)
}

#' @export
## Plotting function for transition probabilities
plot.dtms_probs <- function(x,...) {
  dtms_plot(probs=x,...)
}
