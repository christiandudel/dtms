#'  Cleans data in transition format
#'
#' @description
#' Cleans data in transition format. It can handle issues regularly occurring
#' with such data: transitions starting from or ending in missing states,
#' observations not covered by the time range, transitions starting or ending in
#' a state which is not in the state space, and observations starting in
#' absorbing states.
#'
#' @details
#' Transitions starting or ending with a missing state often occur for three
#' reasons. First, the function \code{dtms_format} will always create a transition
#' with a missing receiving state for the last observation of a unit, whether
#' due to censoring or not. For instance, if t=20 is the last value of the
#' time scale, and a unit is in state A at that time, then there will be a
#' transition starting at time t=20 from state A, and with receiving
#' state missing. Such transitions can usually be safely ignored, in particular
#' if there is only one absorbing state. Second, if, say, for a unit the last
#' observation is at time t=10 and censored after, there will be a transition
#' starting at time t=10 with missing receiving state. Whether such transitions
#' can be ignored depends on the censoring mechanism. If it is uninformative
#' these transitions can be dropped. Third, there might be missing values in
#' the sequence of states. For instance, a unit might first be in state A, then
#' state B, then the state is missing, and then state is again A, giving the
#' sequence A, B, NA, A. This implies a transition from B to NA, and from NA to
#' A.
#'
#' Transitions which are out of the time range can occur, for instance, when
#' the researcher is interested in a shorter time frame than covered by data
#' collection. In a clinical trial, the time scale could capture follow-up
#' time since start of the trial in months and data might be available for 60
#' months. But perhaps the researcher is only interested in the first 36 months.
#'
#' Transitions which start or end in a state which is not in the state space
#' occur when the states in the transition data are not included in the `dtms`
#' object. This likely will apply to states which rarely occur and which the
#' researcher does not want to combine with other states.
#'
#' @param data Data frame, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable with starting state. Default is "from".
#' @param tovar Character (optional), name of variable with receiving state. Default is "to".
#' @param timevar Character (optional), name of variable with time scale. Default is "time".
#' @param dropState Logical (optional), drop transitions with states which are not part of the state space. Default is TRUE.
#' @param dropTime Logical (optional), drop transitions with values of time not covered by the model. Default is TRUE.
#' @param dropNA Logical (optional), drop transitions with gaps, last observations, and similar. Default is TRUE.
#' @param dropAbs Logical (optional), drop transitions starting from absorbing states. Default is TRUE.
#' @param verbose Logical (optional), print how many transitions were dropped. Default is TRUE
#'
#' @return Cleaned data frame in transition format.
#' @export
#'
#' @examples
#' # Define model
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:20)
#' # Transiton format, filling implicit missings with explicit missings
#' estdata <- dtms_format(data=simpledata,
#' dtms=simple,
#' idvar="id",
#' timevar="time",
#' statevar="state",
#' fill=TRUE)
#' # Clean data
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)

dtms_clean <- function(data,
                       dtms,
                       fromvar="from",
                       tovar="to",
                       timevar="time",
                       dropTime=TRUE,
                       dropState=TRUE,
                       dropNA=TRUE,
                       dropAbs=TRUE,
                       verbose=TRUE) {

  # Check
  dtms_proper(dtms)

  # Drop observations not in state space
  if(dropState) {
    allstates <- c(dtms$transient,dtms$absorbing,NA)
    whichrows <- unlist(data[,fromvar])%in%allstates & unlist(data[,tovar])%in%allstates
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows not in state space\n")
    }
  }

  # Drop observations not in time range
  if(dropTime) {
    whichrows <- unlist(data[,timevar])%in%dtms$timescale
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows not in time range\n")
    }
  }

  # Drop missing values
  if(dropNA) {
    whichrows <- !is.na(unlist(data[,fromvar])) & !is.na(unlist(data[,tovar]))
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows starting or ending in NA\n")
    }
  }

  # Drop transitions starting in absorbing states
  if(dropAbs) {
    whichrows <- !unlist(data[,fromvar])%in%dtms$absorbing
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows starting in absorbing state\n")
    }
  }

  # Return
  return(data)

}
