#' Generate variable with duration
#'
#' @param data A data set in transition format
#' @param dtms dtms object, as created with \code{dtms}.
#' @param newname Name of new variable
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param idvar Character (optional), name of variable in `data` with unit ID. Default is "id".
#' @param timevar Character (optional), name of variable in `data` with time scale. Default is "time".
#' @param ignoreleft Logical (optional), ignore left censoring and start counting at the first observation? Default is TRUE.
#'
#' @return The data frame specified with 'data' with an additional column.
#' @export
#'
#' @examples
dtms_duration <- function(data,
                          dtms,
                          newname="duration",
                          fromvar="from",
                          idvar="id",
                          timevar="time",
                          ignoreleft=TRUE) {

  # Check dtms only if specified
  if(!is.null(dtms)) dtms_proper(dtms)

  # Sort data
  dataorder <- order(data[,idvar],
                     data[,timevar])
  data <- data[dataorder,]

  # Apply helper ### ISSUE HERE
  tmp <- tapply(data[,c(fromvar,timevar)],
                data[,idvar],
                function(x) dtms_duration_help(states=as.data.frame(x)[,1],
                                                time=as.data.frame(x)[,2],
                                                ignoreleft=ignoreleft,
                                                dtms=dtms))

  # Assign new values
  data[,newname] <- unlist(tmp)

  # Return
  return(unlist(tmp))

}
