#'  Cleans data in transition format
#'
#' @description
#' Cleans data in transition format. It can handle issues regularly occurring
#' with such data: #' transitions starting from or ending in missing states,
#' observations not covered by the time range, and #' observations starting in
#' absorbing states.
#'
#' @details
#' Data in transition format often includes some transitions which do not
#' contribute to the estimation of a multistate model. This includes
#' transitions starting or ending with a missing state; transitions which
#' are out of the time range of the model; and transitions starting in
#' absorbing states.
#'
#' Transitions starting or ending with a missing state often occur for three
#' reasons. First, the function `dtms_format` will always create a transition
#' with a missing receiving state for the last observation of a unit, whether
#' due to #' censoring or not. For instance, if t=20 is the last value of the
#' time scale, and a unit is in state A at that #' time, then there will be a
#' transition starting at time t=20, starting from state A, and with receiving
#' state missing. These transitions can often be safely ignored, in particular
#' if there is only one absorbing state. Second, if, say, for a unit the last
#' observation is at time t=10 and censored after, there will be a transition
#' starting at time t=10 with missing receiving state. Whether such transitions
#' can be ignored depends on the censoring mechanism. If it is uninformative
#' these transitions can be dropped. Third, there might be missing values in
#' the sequence of states. For instance, a unit might first be in state A, then
#' state B, then the state is missing, and then state is again A, giving the
#' sequence A, B, NA, A. This implies a transition from B to NA, and from NA to
#' A. The function `dtms_format` by default will ignore missing observations
#' and will not create transitions based on them, but it will if the argument
#' `fill` is set to TRUE. Creating these transitions can lead to data sets
#' which are much larger than the original data. If they are created,
#' `dtms_clean` can handle them.
#'
#' Transitions which are out of the time range can occur, for instance, when
#' the researcher is interested in a shorter time frame than covered by data
#' collection. In a clinical trial, the time scale could capture follow-up
#' time since start of the trial in months and data might be available for 60
#' months. But perhaps the researcher is only interested in the first 36 months.
#'
#'
#' @param data Data frame (or similar) in transition format as created with `dtms_format`.
#' @param dtms DTMS object as created with `dtms`
#' @param fromvar Character string (optional), name of the variable with the starting state.
#' @param tovar Character string (optional), name of the variable with the receiving state.
#' @param timevar Character string (optional), name of the variable with the time scale.
#' @param timescale Numeric (optional), values of the time scale.
#' @param absorbing Character (optional), names of absorbing state.
#' @param dropTime Logical, drop transitions with values of time not covered by the model. Default is TRUE.
#' @param dropNA Logical, drop transitions with gaps, last observations, and similar. Default is TRUE.
#' @param dropAbs Logical, drop transitions starting from absorbing states. Default is TRUE.
#' @param verbose Logical, print how many transitions were dropped. Default is true.
#'
#' @return Returns a cleaned data frame in transition format.
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
                       dtms=NULL,
                       fromvar="from",
                       tovar="to",
                       timevar="time",
                       timecale=NULL,
                       absorbing=NULL,
                       dropTime=T,
                       dropNA=T,
                       dropAbs=T,
                       verbose=T) {

  # Use dtms if provided
  if(!is.null(dtms)) {

    # Check
    proper_dtms(dtms)

    # Use values
    timescale <- dtms$timescale
    absorbing <- dtms$absorbing
  }

  # Drop observations not in time range
  if(dropTime) {
    whichrows <- unlist(data[,timevar])%in%timescale
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows because not in time range\n")
    }
  }

  # Drop missing values
  if(dropNA) {
    whichrows <- !is.na(unlist(data[,fromvar])) & !is.na(unlist(data[,tovar]))
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows because gap, last obs, ...\n")
    }
  }

  # Drop transitions starting in absorbing states
  if(dropAbs) {
    whichrows <- !unlist(data[,fromvar])%in%absorbing
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows because starting in absorbing state\n")
    }
  }

}
