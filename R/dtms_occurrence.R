#' Generate variable with number of occurrence of state
#'
#' @description
#' This function creates a variable which measures the number of occurrence
#' of the states.
#'
#' @details
#' Counting starts with 1 and the first occurrence iof a state. For instance,
#' if for an unit the sequence of states A, A, A, B, B, A, C is observed,
#' the occurrence variable would include 1, 1, 1, 1, 1, 2, 1.
#'
#' The argument `ignoreleft` controls how left censoring is handled; i.e., what
#' happens when for a unit there are no observations at the beginning of the
#' time scale. If `TRUE`, left censoring is ignored, and counting starts at
#' the first observation for a unit. For instance, if the time scale starts
#' at t=0, but the first observation for a unit is at time t=2, and the
#' sequence of states is again A, A, A, B, B, A, C, then `ignoreleft=TRUE`
#' returns 1, 2, 3, 1, 2, 1, 1. If `ignoreleft=FALSE`, then the function
#' would return NA, NA, NA, NA, NA, NA, NA for this unit.
#'
#' The function handles gaps in the data by setting all further occurrences to
#' NA. For #' instance, if a unit is observed at times 1, 2, 4, 5, and 6, but
#' not at time #' 3, and the states are A, A, B, C, C, then the occurrence
#' variable will have the values 1, 1, NA, NA, NA. Note that in this case
#' it would be possible to return 1, 1, NA, 1, NA, NA, but the function
#' currently does not have this capability.
#'
#' @param data A data frame in long format.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param statevar Character (optional), name of the variable in the data frame with the states. Default is `state`.
#' @param idvar Character (optional), name of variable with unit ID. Default is `id`.
#' @param timevar Character (optional), name of variable with time scale. Default is `time`.
#' @param newname Character (optional), name of new variable if data set is returned. Default is "duration".
#' @param ignoreleft Logical (optional), ignore left censoring and start counting at the first observation? Default is TRUE.
#' @param vector Logical (optional), return vector (if TRUE) or data frame (if FALSE). Default is FALSE
#'
#' @return The data frame specified with `data` with an additional column (if `vector=FALSE`) or a vector (if `vector=TRUE`).
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:19)
#'
#' dtms_occurrence(data=simpledata,
#'                 dtms=simple)
dtms_occurrence <- function(data,
                            dtms,
                            newname="occurrence",
                            statevar="state",
                            idvar="id",
                            timevar="time",
                            ignoreleft=TRUE,
                            vector=FALSE) {

  # Check dtms only if specified
  if(!is.null(dtms)) dtms_proper(dtms)

  # Sort data
  dataorder <- order(data[,idvar],
                     data[,timevar])
  data <- data[dataorder,]

  # Apply helper
  tmp <- tapply(data[,c(statevar,timevar)],
                data[,idvar],
                function(x) dtms_occurrence_help(states=x[,1],
                                                 time=x[,2],
                                                 ignoreleft=ignoreleft,
                                                 dtms=dtms))

  # Return vector
  if(vector) return(unlist(tmp))

  # Assign new values
  data[,newname] <- unlist(tmp)

  # Return
  return(data)

}
