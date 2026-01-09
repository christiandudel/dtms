#' Bootstrap and block bootstrap
#'
#' @description
#' This function is a simple wrapper for bootstrapping and block-bootstrapping
#' data in transition format. Parallel processing is supported.
#'
#' @details
#' \code{dtms_boot()} takes a function specified with the argument `fun` and
#' applies it several times to resampled data, where the original data is
#' specified with the argument `data` and `rep` specifies the number of
#' replications. The argument `dtms` takes an object created with \code{dtms()}
#' and also passes it to `fun`. `data` is passed to `fun` as its first argument,
#' and `dtms` is passed as the second argument.
#'
#' The result of this function is a list with as many entries as there are
#' replications. Each entry is the result of calling `fun` for the respective
#' replication.
#'
#' Three methods are implemented and selected with the argument `method`. A simple
#' resampling bootstrap, which assumes that the rows in `data` are independent
#' of each other (`method=simple`). The block bootstrap which allows for
#' dependent observations; e.g., different units each contributing
#' several transitions (`method=block`). Moreover, a parametric bootstrap using
#' weights is also supported, assuming that observations are i.i.d. multinomial
#' (`method=weights`). If the block bootstrap is used the argument `idvar` sets
#' which variable in #' `data` contains information the unit/cluster identifier.
#' In case the parametric bootstrap is used the argument `weights` is used to
#' specify the name of the variable with the weights.
#'
#' For parallel computing the packages `foreach` and `doParallel` (and their
#' dependcies are used). See the documentation of these packages for details.
#'
#' @param data Data frame in transition format as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fun Function to be repeatedly applied, see details.
#' @param rep Numeric, number of bootstrap replications.
#' @param method Character (optional), either "simple" for simple bootstrap, "block" for block bootstrap, or "weights" for a weight-based parametric bootstrap. Default is "simple".
#' @param idvar Character (optional), name of ID variable in `data' identifying units. Only required for block bootstrap. Default is "id".
#' @param weights Character (optional), name of variable with weights. Only used if `method=weights`. Default is NULL.
#' @param slack Numeric (optional), used to in parametric resampling to replace 0. Default is 1.
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
#' @seealso
#' \code{\link{dtms_boot_summary}} to help with summarizing the results.
#'
#' @examples
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
#' # Simple resampling bootstrap
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
#' # Parametric bootstrap
#' aggdata <- dtms_aggregate(estdata)
#' # Bootstrap function
#' bootfun <- function(data,dtms) {
#'   fit <- dtms_fit(data=data,weights="count")
#'   probs <- dtms_transitions(dtms=dtms,
#'                             model = fit)
#'   S <- dtms_start(dtms=dtms,
#'                   data=data,
#'                   weights="count")
#'   dtms_expectancy(dtms=dtms,
#'                   probs=probs,
#'                   start_distr=S)
#' }
#' # Bootstrap
#' bootstrap <- dtms_boot(data=aggdata,
#'                        dtms=simple,
#'                        fun=bootfun,
#'                        rep=5,
#'                        weights="count",
#'                        method="weights")
#' # Results
#' summary(bootstrap,
#'         probs=c(0.025,0.5,0.975))


dtms_boot <- function(data,
                      dtms,
                      fun,
                      rep,
                      method="simple",
                      idvar="id",
                      weights=NULL,
                      slack=1,
                      verbose=FALSE,
                      progress=FALSE,
                      parallel=FALSE,
                      cores=2,
                      .packages=c("mclogit","VGAM","nnet","dtms"),
                      ...) {

  # Check if dtms correctly specified
  dtms_proper(dtms)

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

  if(method=="weights") n <- sum(data[,weights])

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
      }

      # Simple parametric
      if(method=="weights") {
        newweights <- stats::rmultinom(1,
                                       size=sum(data[,weights]),
                                       prob=data[,weights])
        newweights <- as.numeric(newweights)
        # Replace 0 values
        if(any(newweights==0)) {
          newweights[newweights==0] <- slack
          newweights <- newweights* (n/sum(newweights))
        }
        newweights <- as.numeric(newweights)
        newdata[,weights] <- newweights
      }

      # Run
      fun(newdata,dtms,...)

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
    }

    # Parametric/weight-based
    if(method=="weights") {
      # Resample
      newweights <- stats::rmultinom(1,
                                     size=sum(data[,weights]),
                                     prob=data[,weights])
      newweights <- as.numeric(newweights)
      # Replace 0 values
      if(any(newweights==0)) {
        newweights[newweights==0] <- slack
        newweights <- newweights* (n/sum(newweights))
      }
      newdata <- data
      newdata[,weights] <- newweights
    }

    if(verbose) result[[i]] <- fun(newdata,dtms,...)
    if(!verbose) utils::capture.output(result[[i]] <- fun(newdata,dtms,...),file=nullfile())

  }

  # Return
  class(result) <- c("dtms_boot","list")
  return(result)

}
