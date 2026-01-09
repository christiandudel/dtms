#' Summary function for bootstrap results
#'
#' @description
#' Provides bootstrap percentiles for bootstrap replications created with
#' \code{dtms_boot()}.
#'
#' @details
#' Percentiles can be specified with the argument \code{probs}. This can be as
#' many percentiles as required by the user. If it is not
#' specified, the argument \code{alpha} is used instead. \code{alpha} is the
#' confidence level for the confidence intervals.
#'
#' The function passed to \code{dtms_boot()} needs to either return a numeric
#' vector, a matrix, or a data.frame, otherwise \code{dtms_boot_summary()}
#' returns an error #' message. A data.frame will be transformed into a matrix.
#'
#' @param boot Object created with \code{dtms_boot()}.
#' @param probs Numeric (optional), vector of percentiles. Default is NULL.
#' @param alpha Numeric (optional), confidence level. Default is 0.05.
#'
#' @return Either a vector or a matrix.
#' @export
#'
#' @examples
#'
#' ## Define model: Absorbing and transient states, time scale
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:20)
#' ## Reshape to transition format
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' ## Clean
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' # Bootstrap function
#' bootfun <- function(data,dtms) {
#'   fit <- dtms_fit(data=data)
#'   probs    <- dtms_transitions(dtms=dtms,
#'                                model = fit)
#'   S <- dtms_start(dtms=dtms,
#'                   data=data)
#'   dtms_expectancy(dtms=dtms,
#'                   probs=probs,
#'                   start_distr=S)
#' }
#' # Run bootstrap
#' bootstrap <- dtms_boot(data=estdata,
#'                        dtms=simple,
#'                        fun=bootfun,
#'                        rep=5)
#' summary(bootstrap,
#'         probs=c(0.025,0.5,0.975))

dtms_boot_summary <- function(boot,probs=NULL,alpha=0.05) {

  # Get class
  what <- class(boot[[1]] )[1]

  # Get probs
  if(is.null(probs)) probs <- c(alpha/2,1-alpha/2)
  probs <- round(probs,digits=3)

  # How many results
  nresults <- length(probs)

  # Check
  if(!what%in%c("matrix","numeric","data.frame")) stop("Format not supported")

  # If matrix
  if(what%in%c("matrix","data.frame")) {

    # Turn into matrix
    if(what=="data.frame") {
      boot <- lapply(boot, data.matrix)
    }

    # Format result
    dimwhat <- dim(boot[[1]])
    colwhat <- colnames(boot[[1]])
    rowwhat <- rownames(boot[[1]])
    result <- list()
    entry <- matrix(NA,nrow=dimwhat[1],ncol=dimwhat[2])
    colnames(entry) <- colwhat
    rownames(entry) <- rowwhat
    for(i in 1:nresults) result[[i]] <- entry
    names(result) <- paste0(probs*100,"%")

    # Go through entries
    for(row in 1:dimwhat[1]) {
      for(col in 1:dimwhat[2]) {
        values <- unlist(lapply(boot,function(x) x[row,col]))
        values <- stats::quantile(values,probs=probs)
        for(val in names(values)) result[[val]][row,col] <- values[val]
      }
    }

    return(result)

  }

  # If vector
  if(what=="numeric") {

    # Format result
    nwhat <- length(boot[[1]])
    nameswhat <- names(boot[[1]])
    result <- list()
    for(i in 1:nresults) {
      result[[i]] <- numeric(length=nwhat)
      names(result[[i]]) <- nameswhat
    }
    names(result) <- paste0(probs*100,"%")

    # Go through results
    for(i in 1:nwhat) {
      values <- unlist(lapply(boot,function(x) x[i]))
      values <- stats::quantile(values,probs=probs)
      for(val in names(values)) result[[val]][i] <- values[val]

    }

    return(result)

  }

}
