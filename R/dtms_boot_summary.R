#' Summary function for bootstrap results
#'
#' @description
#' Provides bootstrap percentiles for \code{dtms_boot()}.
#'
#' @details
#' Percentiles can be specified with the argument \code{probs}. If it is not
#' specified, the argument \code{alpha} is used instead. \code{alpha} is the
#' confidence level for
#'
#' @param boot Object created with \code{dtms_boot()}.
#' @param probs Numeric (optional), vector of percentiles. Default is NULL.
#' @param alpha Numeric (optional), confidence level. Default is 0.05.
#'
#' @return A list
#' @export

dtms_boot_summary <- function(boot,probs=NULL,alpha=0.05) {

  # Get class
  what <- class(boot[[1]] )[1]

  # Get probs
  if(is.null(probs)) probs <- c(alpha/2,1-alpha/2)
  probs <- round(probs,digits=3)

  # How many results
  nresults <- length(probs)

  # Check
  if(!what%in%c("matrix","numeric")) stop("Format not supported")

  # If matrix
  if(what=="matrix") {

    # Format result
    dimwhat <- dim(boot[[1]])
    result <- list()
    entry <- matrix(NA,nrow=dimwhat[1],ncol=dimwhat[2])
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
    result <- list()
    for(i in 1:nresults) result[[i]] <- numeric(length=nwhat)
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
