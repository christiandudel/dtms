#' Summarize data in transition format
#'
#' @description
#' Returns a data frame with number of observed transitions (column COUNT),
#' relative proportion (column PROP), and raw transition probabilities (column
#' PROB).
#'
#' @param data Data frame, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable with starting state. Default is "from".
#' @param tovar Character (optional), name of variable with receiving state. Default is "to".
#' @param weights Character (optional), name of variable with weights. Default is NULL.
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
#' dtms_data_summary(estdata)

## Method
dtms_data_summary <- function(data,
                              dtms=NULL,
                              fromvar="from",
                              tovar="to",
                              weights=NULL) {

    # Weights per transition
    if(is.null(weights)) data$COUNT <- 1 else
      data <- dtms_rename(data,weights,"COUNT")

    # For handling of missing values
    data[is.na(data[,fromvar]),fromvar] <- "NA"
    data[is.na(data[,tovar]),tovar] <- "NA"

    # Aggregate
    formal <- paste0("COUNT~",fromvar,"+",tovar)
    formal <- stats::as.formula(formal)
    result <- stats::aggregate(formal,data,FUN=sum,drop=F)

    # If there are unused combinations COUNT could be NA
    result[is.na(result$COUNT),"COUNT"] <- 0

    # Order
    ordering <- order(result[,fromvar],result[,tovar])
    result <- result[ordering,]

    # Proportion
    N <- sum(result$COUNT)
    result$PROP <- result$COUNT/N

    # Raw transition probabilities
    probs <- tapply(result$COUNT,
                    result[,fromvar],
                    FUN=function(x) x/sum(x))
    probs <- unlist(probs)
    result$PROB <- probs

    # Row-numbers
    rownames(result) <- NULL

    # If dtms is provided
    if(!is.null(dtms))  {

      # Get data frame
      newframe <- expand.grid(from=dtms$transient,to=c(dtms$transient,dtms$absorbing))

      # Result
      newresult <- merge(newframe,result,all=TRUE)

      # Replace missings with zero
      newresult[is.na(newresult$COUNT),c("COUNT","PROP","PROB")] <- rep(0,3)

      # Warning
      if(any(!result[,fromvar]%in%c(dtms$transient,dtms$absorbing))) warning("Some fromvar values not in state space")
      if(any(!result[,tovar]%in%c(dtms$transient,dtms$absorbing))) warning("Some tovar values not in state space")

    }

    # Return
    return(result)
}
