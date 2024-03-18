#' Summarize transition probabilities
#'
#' @description
#' Provides several summary statistics on transition probabilities.
#'
#' @param probs Object with transition probabilities as created with \code{dtms_transitions}.
#' @param fromvar Character (optional), name of variable with starting state. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state. Default is `to`.
#' @param timevar Character (optional), name of variable with time scale. Default is `time`.
#' @param Pvar Character (optional), name of variable with transition probabilities. Default is `P`.
#' @param digits Numeric (optional), number of digits to return, default is 6.
#' @param format Character (optional), show results in decimal format or percentage, either `decimal` or `percent`. Default is `decimal`.
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#'
#' @return A data frame
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:20)
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
#' summary(probs)

dtms_probs_summary <- function(probs,
                               fromvar="from",
                               tovar="to",
                               timevar="time",
                               Pvar="P",
                               sep="_",
                               digits=4,
                               format="decimal") {

  # Get short state names
  probs[,fromvar] <- dtms_getstate(probs[,fromvar],sep=sep)
  probs[,tovar] <- dtms_getstate(probs[,tovar],sep=sep)

  # Aggregate (starting from minimum)
  result <- stats::aggregate(probs[,Pvar]~probs[,fromvar]+
                                          probs[,tovar],
                             FUN=min)
  names(result) <- c(fromvar,tovar,"MIN")

  # Add time values
  result$MINtime <- probs[match(result$MIN,probs[,Pvar]),timevar]

  # Add max
  result$MAX <- stats::aggregate(probs[,Pvar]~probs[,fromvar]+
                                              probs[,tovar],
                                 FUN=max)[,3]
  result$MAXtime <- probs[match(result$MAX,probs[,Pvar]),timevar]

  # Add other statistics
  result$MEDIAN <- stats::aggregate(probs[,Pvar]~probs[,fromvar]+
                                                 probs[,tovar],
                                    FUN=stats::median)[,3]
  result$MEAN <- stats::aggregate(probs[,Pvar]~probs[,fromvar]+
                                               probs[,tovar],
                                  FUN=mean)[,3]

  # Order result
  ordering <- order(result[,fromvar],result[,tovar])
  result <- result[ordering,]

  # Rounding
  result[,c("MIN","MAX","MEDIAN","MEAN")] <-
    round(result[,c("MIN","MAX","MEDIAN","MEAN")],digits=digits)

  # For printing
  if(format=="percent") {
    result[,c("MIN","MAX","MEDIAN","MEAN")] <-
      result[,c("MIN","MAX","MEDIAN","MEAN")]*100
  }

  # Return
  return(result)

}
