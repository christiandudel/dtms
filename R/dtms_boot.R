#' Bootstrap and block bootstrap
#'
#' @description
#' This function is a simple wrapper for bootstrapping and block-bootstrapping
#' data in transition format.
#'
#' @details
#' This is currently a placeholder. What goes here: fun(arguments -> 1st data, 2nd dtms), seed (if null), progress, parallel, results (list with rep entries)
#'
#' @param data Data frame in transition format as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fun Function to be repeatedly applied, see details.
#' @param rep Numeric, number of bootstrap replications.
#' @param method Character (optional), either "simple" for simple bootstrap or "block" for block bootstrap. Default is "simple".
#' @param idvar Character (optional), name of ID variable in `data' identifying units. Only required for block bootstrap. Default is "id".
#' @param seed Numeric (optional), seed for random numbers. See details.
#' @param verbose Logical (optional), print output which might be generated when running `fun`? Default is FALSE.
#' @param progress Logical (optional), indicate progress if simple bootstrap? Default is FALSE.
#' @param parallel Logical (optional), use parallel processing? Default is FALSE.
#' @param cores Numeric (optional), if parallel=TRUE, how many cores should be used? Default is 2.
#' @param .packages Character (optional), packages to be loaded when parallel processing. Default is `c("mclogit","VGAM","nnet","dtms")`
#' @param ... Arguments to be passed to `fun`, only works if `parallel=FALSE`.
#'
#' @return A list of results, see details
#' @export
#'
#' @examples
#' #Currently none

dtms_boot <- function(data,
                      dtms,
                      fun,
                      rep,
                      method="simple",
                      idvar="id",
                      seed=NULL,
                      verbose=FALSE,
                      progress=FALSE,
                      parallel=FALSE,
                      cores=2,
                      .packages=c("mclogit","VGAM","nnet","dtms"),
                      ...) {

  # Seed
  if(is.null(seed)) set.seed(as.numeric(Sys.time()))

  # For resampling
  if(method=="simple") n <- dim(data)[1]

  if(method=="block") {

    # Unique ids
    ids <- unique(data[,idvar])
    nids <- length(ids)

    # Rows belonging to each id
    rowsid <- vector("list",nids)
    names(rowsid) <- paste(ids)
    for(i in ids) {
      rowsid[[paste(i)]] <- which(data[,idvar]==i)
    }

  }

  # Parallel version
  if(parallel) {

    # Avoids import
    `%dopar%` <- foreach::`%dopar%`

    # Set up parallel processing
    cluster <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cluster)

    # Foreach "loop"
    result <- foreach::foreach(sample_no=1:rep,.packages=.packages,.verbose=verbose) %dopar% {

       # Simple resampling
      if(method=="simple") {
        samplerows <- sample(1:n,size=n,replace=TRUE)
        newdata <- data[samplerows,]
      }

      # Block bootstrap
      if(method=="block") {
        sampleids <- sample(ids,size=nids,replace=TRUE)
        samplerows <- unlist(rowsid[paste(sampleids)])
        newdata <- data[samplerows,]
        dim(newdata)
      }

      fun(newdata,dtms)

    }

    # Stop cluster
    parallel::stopCluster(cluster)

    # Return
    class(result) <- c("dtms_boot","list")
    return(result)

  }

  # Not parallel, list for results
  result <- list()

  # Repeat rep times
  for(i in 1:rep) {

    # Show progress
    if(progress & i%%10==0) cat(".")

    # Simple resampling
    if(method=="simple") {
      samplerows <- sample(1:n,size=n,replace=TRUE)
      newdata <- data[samplerows,]
    }

    # Block bootstrap
    if(method=="block") {
      sampleids <- sample(ids,size=nids,replace=TRUE)
      samplerows <- unlist(rowsid[paste(sampleids)])
      newdata <- data[samplerows,]
      dim(newdata)
    }

    if(verbose) result[[i]] <- fun(newdata,dtms,...) else
      utils::capture.output(result[[i]] <- fun(newdata,dtms,...),file=nullfile())

  }

  # Return
  class(result) <- c("dtms_boot","list")
  return(result)

}
