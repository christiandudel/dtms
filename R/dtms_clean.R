#'  Cleans data in transition format
#'
#' @description
#' Cleans data in transition format. It can handle gaps in sequences and final observations,
#' observatiosn not covered by the time range, and observations starting in absorbing states.
#'
#' @param data Data frame (or similar) in transition format as created with `dtms_format`.
#' @param dtms DTMS object as created with `dtms`
#' @param fromvar Character string (optional), name of the variable with the starting state.
#' @param tovar Character string (optional), name of the variable with the receiving state.
#' @param timevar Character string (optional), name of the variable with the time scale.
#' @param dropTime Logical, drop observations with values of time not covered by the model. Default is TRUE.
#' @param dropNA Logical, drop observations with gaps, last observations, and similar. Default is TRUE.
#' @param dropAbs Logical, drop observations starting from absorbing states. Default is TRUE.
#'
#' @return Returns a cleaned data frame in transition format.
#' @export
#'
#' @examples
#' # Define model
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' time=0:20)
#' # Transiton format
#' estdata <- dtms_format(data=simpledata,
#' dtms=simple,
#' idvar="id",
#' timevar="time",
#' statevar="state")
#' # Clean data
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
dtms_clean <- function(data, # Transition format data
                       dtms, # dtms object
                       fromvar="from", # Variable with starting state
                       tovar="to", # Variable with receiving state
                       timevar="time", # Variable name of time
                       dropTime=T, # Drop obs with time not in the model
                       dropNA=T, # Drop gaps, last obs, ...
                       dropAbs=T) { # Drop obs starting from absorbing state

  if(dropTime) {
    whichrows <- unlist(data[,timevar])%in%dtms$time
    count <- sum(!whichrows)
    cat("Dropping ",count," rows because not in time range\n")
    data <- data[whichrows,]
  }

  if(dropNA) {
    whichrows <- !is.na(unlist(data[,fromvar])) & !is.na(unlist(data[,tovar]))
    count <- sum(!whichrows)
    cat("Dropping ",count," rows because gap, last obs, ...\n")
    data <- data[whichrows,]
  }

  if(dropAbs) {
    whichrows <- !unlist(data[,fromvar])%in%dtms$absorbing
    count <- sum(!whichrows)
    cat("Dropping ",count," rows because starting in absorbing state\n")
    data <- data[whichrows,]
  }

}
