#' Summary for distributional results
#'
#' @description
#' This function provides several summary measures based on results obtained
#' with dtms_visits, dtms_first, and dtms_last.
#'
#' @param distr An object created with dtms_visits, dtms_first, or dtms_last.
#'
#' @return A matrix with summary measures
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:20)
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' fit <- dtms_fit(data=estdata)
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' example <- dtms_visits(dtms=simple,
#'                        matrix=Tp,
#'                        risk="A")
#' summary(example)


dtms_distr_summary <- function(distr) {

  # Drop totals and similar
  allnames <- colnames(distr)
  allnames <- as.numeric(allnames)
  whichdrop <- is.na(allnames)
  distr <- distr[,!whichdrop]

  # Get values
  values <- colnames(distr)
  values <- as.numeric(values)

  # Mean
  MEAN <- apply(distr,1,function(x) sum(x*values))

  # Variance
  VARIANCE <- apply(distr,1,function(x) sum((values-sum(x*values))^2*x))

  # Standard deviation
  SD <- sqrt(VARIANCE)

  # Median
  MEDIAN <- apply(distr,1,function(x) values[min(which(cumsum(x)>0.5))])

  # Risk
  RISK0 <- apply(distr,1,function(x) x["0"])

  # Combine
  result <- cbind(MEAN,VARIANCE,SD,MEDIAN,RISK0)

  # Return
  return(result)

}
