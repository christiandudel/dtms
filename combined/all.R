#' Create dtms object
#'
#' @description
#' This function creates an object of class 'dtms' to be passed to other
#' functions of the package.
#'
#' @details
#' \code{dtms} provides an abstract definition of a multistate model, including
#' the names of the transient states, the names of the absorbing states, the
#' values the time scale can take, and the step length of the time scale.
#'
#' The names of the absorbing and transient states should be provided as
#' character strings. However, numeric values also work. Factors are not
#' supported
#'
#' The step length of the time scale can be a vector with several values, which
#' allows for unevenly spaced observations. Note, however, that some functions
#' require one specific value for the step length; e.g.,
#' \code{dtms_transitions()}. For such functions, if several values are provided
#' the first value will be used.
#'
#' @param transient A character vector of names of the transient states in the state space.
#' @param absorbing A character vector of names of the absorbing states in the state space.
#' @param timescale A numeric vector with the time scale, including the starting time and the final time.
#' @param timestep Numeric (optional), step length of the time scale, will be guessed if NULL (default).
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#'
#' @return Returns an object of class 'dtms'
#' @export
#'
#' @examples
#' dtms(transient=c("A","B"),
#'      absorbing="X",
#'      timescale=1:10)

dtms <- function(transient,
                 absorbing,
                 timescale,
                 timestep=NULL,
                 sep="_") {

  # Guess time step?
  if(is.null(timestep)) {

    # Step lengths as specified by time scale
    step <- timescale |> diff() |> unique()

    # Number of different step lengths
    nstep <- length(step)

    # Guess timestep if possible, else error
    if(nstep==1) timestep <- step else
      stop("Not able to guess step length of time scale")

  }

  # Combine everything in a list
  result <- list(transient=transient,
                 absorbing=absorbing,
                 timescale=timescale,
                 timestep=timestep,
                 sep=sep)

  # Assign class
  class(result)[2] <- "dtms"

  # Return
  return(result)

}
#' Calculate the distribution of the time until entering an absorbing state
#'
#' @description
#' Calculates the distribution of the time until entering any of the absorbing
#' states.
#'
#' @details
#' In a discrete-time model, the time spent in a state depends on assumptions
#' about when transitions happen. Currently, this functions supports two
#' variants which can be specified with the argument `method`: mid-interval
#' transitions can be selected with the option `mid` and imply that transitions
#' happen at the middle of the time interval; and the option `end` assumes
#' that instead transitions happen at the end of the interval. The calculation
#' takes the step length of the time scale into account as specified by the
#' `dtms` object. If the #' step length is not one fixed value, the first entry
#' of `dtms$timestep` will be used.
#'
#' If a distribution of the starting states is provided with `start_distr` the
#' output table has an additional row, showing the waiting time distribution
#' unconditional on the starting state.
#'
#' @param matrix Matrix with transition probabilities, as generated with \code{dtms_matrix}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param start_distr Numeric (optional), distribution of starting states. If specified, average distribution over all starting states will be calculated.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#' @param method Character (optional), do transitions happen mid-interval (`mid`, default) or at the end of the interval (`end`), see details.
#'
#' @return A table with the distribution of time until entering any of the absorbing states.
#' @export
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
#' ## Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Get transition matrix
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' ## Get starting distribution
#' S <- dtms_start(dtms=simple,
#'                 data=estdata)
#' ## Distribution of visits
#' dtms_absorbed(dtms=simple,
#'               matrix=Tp,
#'               start_distr=S)


dtms_absorbed <- function(matrix,
                          dtms,
                          start_distr=NULL,
                          start_state=NULL,
                          start_time=NULL,
                          end_time=NULL,
                          method="mid") {

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=dtms$sep)

  # Matrix for results
  Umat <- dtms_absorbing(matrix=matrix)
  nstates <- dim(Umat)[1]
  distr <- matrix(data=NA,ncol=1,nrow=nstates)
  rownames(distr) <- rownames(Umat)

  # First step: Pr(Time to absorbtion=1)
  distr[,1] <- rowSums(diag(1,nstates)-Umat)
  step <- 1
  if(is.null(end_time)) max_steps <- length(dtms$timescale) else {
    max_steps <- length(1:which(dtms$timescale==end_time))
  }

  # Iterate
  while(!all(round(colSums(distr),digits=2)==1) & step<max_steps) {
    step <- step+1
    distr <- cbind(distr,rowSums(dtms_mtexp(Umat,step-1)-dtms_mtexp(Umat,step)))
  }

  # Name
  if(method=="mid") colnames(distr) <- paste((c(1:max_steps)-0.5)*dtms$timestep[1])
  if(method=="end") colnames(distr) <- paste(c(1:max_steps)*dtms$timestep[1])

  # Select starting states
  result <- distr[starting,]

  # Average
  if(!is.null(start_distr)) {

    # Overall average
    AVERAGE <- start_distr%*%result
    rownames(AVERAGE) <- "AVERAGE"
    result <- rbind(result,AVERAGE)

  }

  # Assign class
  class(result) <- c("dtms_distr","matrix")

  # Results
  return(result)

} # End of function
#' Aggregate data
#'
#' @description
#' This function takes any data set and returns a new data frame which only
#' includes the unique rows from the original data set and indicates how
#' often these rows appear in the original data.
#'
#' @details
#' Currently, missing values are not supported and will be dropped; consider
#' using factors if you want to keep them. The variable provided with the
#' argument `idvar` will be dropped from the aggregated data. If `weights` is
#' specified, the counts will be placed in a variable with the same name. If
#' `countvar` is used, any existing variable in the original data with the
#' same name will be replaced.
#'
#' @param data Data frame.
#' @param weights Character (optional). Name of variable with weights.
#' @param idvar Character (optional). Name of variable in `data` with unit ID. Default is "id".
#' @param countvar Character (optional). Name of new variable in data which provides the counts. Default is "count".
#'
#' @returns An aggregated data frame
#' @export
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
#' ## Aggregate
#' aggdata <- dtms_aggregate(estdata)
#' ## Fit model
#' fit <- dtms_fit(data=aggdata,
#'                 weights="count")

dtms_aggregate <- function(data,
                           weights=NULL,
                           idvar="id",
                           countvar="count") {

  # Transform to data frame, e.g., if tibble
  if(class(data)[1]!="data.frame") data <- as.data.frame(data)

  # Get variable names without weights, collapse for formula
  controls <- names(data)
  controls <- controls[which(controls!=idvar)]
  if(!is.null(weights)) controls <- controls[which(controls!=weights)]
  controls <- paste(controls,collapse="+")

  # If no weights
  if(is.null(weights)) {
    data[,countvar] <- 1
    aggformula <- paste0(countvar,"~",controls)

  # If weights
  } else {
    aggformula <- paste0(weights,"~",controls)
  }

  # Formula for aggregate
  aggformula <- stats::as.formula(aggformula)

  # Warning if missing values
  drops <- data |>
    stats::na.omit() |> dim()
  drops <- dim(data)[1]-drops[1]
  if(drops>0) warning(paste("Dropping",drops,"rows with missing values"))

  # Aggregate
  tmp <- stats::aggregate(by=aggformula,
                          x=data,
                          FUN=sum)

  # Return
  return(tmp)

}
#' Carry states backward
#'
#' @description
#' This function carries a state backward after its last occurrence.
#'
#' @details
#' This function carries a state backward after its first occurrence.
#' For instance, carrying the state "A" backward in the sequence `B, B, A, B, B`
#' will give the sequence `A, A, A, B, B`. The sequence `C, B, C, A, B, A, A, B`
#' will give `A, A, A, A, A, A, A, B`.
#'
#' This function works with data frames in transition format and in long format.
#' The default is transition format, using the arguments `fromvar` and `tovar`.
#' If, however, the argument `statevar` is specified, it is used instead.
#'
#' The argument `overwrite` is used to control what type of information is
#' replaced. If `overwrite==transient`, then only transient states are replaced
#' while missing values and absorbing states remain unchanged. For example,
#' carrying backward state "A" in the sequence `B, NA, B, B, X, A, X` with X
#' being an absorbing state will give `A, NA, A, A, X, A, X`. If
#' `overwrite==missing` then in addition to transient states also missing values
#' are replaced and for the example sequence `A, A, A, A, X, A, X` would be
#' returned. If `overwrite==absorbing` then in addition to transient states
#' absorbing states will be replaced; for the example sequence the result would
#' be `A, NA, A, A, A, A, X`. Finally, if `overwrite==all` then all values in
#' the sequence will be replaced: `A, A, A, A, A, A, X`.
#'
#' @seealso
#' \code{\link{dtms_forward}} to carry states forward.
#'
#' @param data A data frame in long format.
#' @param state Character, name of the state to be carried forward.
#' @param fromvar Character (optional), name of variable with starting state. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state. Default is `to`.
#' @param statevar Character (optional), name of the variable in the data frame in long format with the states. Default is NULL.
#' @param idvar Character (optional), name of variable with unit ID. Default is `id`.
#' @param timevar Character (optional), name of variable with time scale. Default is `time`.
#' @param dtms dtms object (optional), as created with \code{dtms}. Not required if `overwrite==transient`.
#' @param overwrite Character (optional), one of `transient`, `missing`, `absorbing`, and `all`, see details. Default is `transient`.
#' @param vector Logical (optional), return vector (if TRUE) or data frame (if FALSE). Default is FALSE. Argument is only used if argument `statevar` is specified.
#'
#' @return The data frame specified with `data` and the edited state variable (if `vector=FALSE`) or a vector (if `vector=TRUE`).
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:19)
#'
#' dtms_backward(data=simpledata,
#'               statevar="state",
#'               state="A",
#'               dtms=simple,
#'               overwrite="transient")

dtms_backward <- function(data,
                         state,
                         fromvar="from",
                         tovar="to",
                         statevar=NULL,
                         idvar="id",
                         timevar="time",
                         dtms=NULL,
                         overwrite="missing",
                         vector=FALSE) {

  # Check dtms only if specified
  if(!is.null(dtms)) dtms_proper(dtms)

  # Sort data
  dataorder <- order(data[,idvar],
                     data[,timevar])
  data <- data[dataorder,]

  # Apply helper to transition format
  if(is.null(statevar)) {

    # Move receiving state
    data[,fromvar][data[,tovar]==state] <- state

    # Carry backward starting state
    tmp1 <- tapply(data[,fromvar],
                  data[,idvar],
                  function(x) dtms_backward_help(x=x,
                                                 state=state,
                                                 overwrite=overwrite,
                                                 dtms=dtms))

    # Carry backward receiving state
    tmp2 <- tapply(data[,tovar],
                   data[,idvar],
                   function(x) dtms_backward_help(x=x,
                                                  state=state,
                                                  overwrite=overwrite,
                                                  dtms=dtms))

    # Assign new values
    data[,fromvar] <- unlist(tmp1)
    data[,tovar] <- unlist(tmp2)

  # Apply to long format
  } else {

    # Apply helper
    tmp <- tapply(data[,statevar],
                  data[,idvar],
                  function(x) dtms_backward_help(x=x,
                                                state=state,
                                                overwrite=overwrite,
                                                dtms=dtms))

    # Return vector
    if(vector) return(unlist(tmp))

    # Assign new values
    data[,statevar] <- unlist(tmp)

  }

  # Return
  return(data)

}
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
#' Two methods are implemented and selected with the argument `method`. A simple
#' resampling bootstrap, which assumes that the rows in `data` are independent
#' of each other. And the block bootstrap which allows for dependent
#' observations; e.g., different units each contributing several transitions.
#' If the block bootstrap is used the argument `idvar` sets which variable in
#' `data` contains information the unit/cluster identifier.
#'
#' For parallel computing the packages `foreach` and `doParallel` (and their
#' dependcies are used). See the documentation of these packages for details.
#'
#' @param data Data frame in transition format as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fun Function to be repeatedly applied, see details.
#' @param rep Numeric, number of bootstrap replications.
#' @param method Character (optional), either "simple" for simple bootstrap or "block" for block bootstrap. Default is "simple".
#' @param idvar Character (optional), name of ID variable in `data' identifying units. Only required for block bootstrap. Default is "id".
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
#' # Bootstrap function
#' bootfun <- function(data,dtms) {
#'   fit <- dtms_fit(data=data)
#'   probs    <- dtms_transitions(dtms=dtms,
#'                                model = fit)
#'   Tp <- dtms_matrix(dtms=dtms,
#'                     probs=probs)
#'   S <- dtms_start(dtms=dtms,
#'                   data=data)
#'   dtms_expectancy(dtms=dtms,
#'                   matrix=Tp,
#'                   start_distr=S)
#' }
#' # Run bootstrap
#' bootstrap <- dtms_boot(data=estdata,
#'                        dtms=simple,
#'                        fun=bootfun,
#'                        rep=5)
#' summary(bootstrap,
#'         probs=c(0.025,0.5,0.975))

dtms_boot <- function(data,
                      dtms,
                      fun,
                      rep,
                      method="simple",
                      idvar="id",
                      verbose=FALSE,
                      progress=FALSE,
                      parallel=FALSE,
                      cores=2,
                      .packages=c("mclogit","VGAM","nnet","dtms"),
                      ...) {

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
#'   Tp <- dtms_matrix(dtms=dtms,
#'                     probs=probs)
#'   S <- dtms_start(dtms=dtms,
#'                   data=data)
#'   dtms_expectancy(dtms=dtms,
#'                   matrix=Tp,
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
#' Left censoring, right censoring, and gaps in data
#'
#' @description
#' This function provides an overview of censoring and gaps in the data. It can
#' do so in several ways: by providing counts of units with left censoring,
#' right censoring, and gaps; by providing a cross-tabulation of the number of
#' units with left censoring and/or right censoring and/or gaps; and by
#' returning a data frame with added indicators on censoring and gaps.
#'
#' @details
#' Added variables can be at the unit level or at the observation level. This
#' is controlled by the argument "addtype". If it is set to "id" then the unit
#' level is used. In this case the added variables are the same for
#' each observation  of a unit. For instance, if a unit experiences any gap,
#' then the added variable has the value TRUE for all observations of that unit.
#' If "addtype" is set to "obs" the observation level is used and the indicators
#' are only set to TRUE if they apply to a specific observation. For instance,
#' if a unit experience right censoring, only the last observation will have
#' TRUE as the value for the right-censoring indicator; i.e., showing that after
#' this last observation there is right censoring. This can be helpful for
#' analyses to understand censoring better.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable in `data` with receiving state. Default is "to".
#' @param timevar Character (optional), name of variable in `data` with time scale. Default is "time".
#' @param idvar Character (optional), name of variable in `data` with unit ID. Default is "id".
#' @param print Logical (optional), print counts? Default is TRUE.
#' @param printlong Logical (optional), print cross-tabulation? Default is FALSE.
#' @param add Logical (optional), add indicators to data set? Default is FALSE. If TRUE the data frame specified with \code{data} is returned with added columns.
#' @param addtype Character (optional), what type of information should be added if add=TRUE. Either "id" or "obs", see details. Default is "id".
#' @param varnames Character vector (optional), names of added variables if add=T. Default is "c("LEFT","GAP","RIGHT")".
#'
#' @return Table or data frame.
#' @export
#'
#' @examples
#' ## Define model: Absorbing and transient states, time scale
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:19)
#' # Reshape to transition format
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' ## Clean
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' ## Censoring
#' dtms_censoring(data=estdata,
#'                dtms=simple)

dtms_censoring <- function(data,
                           dtms,
                           fromvar="from",
                           tovar="to",
                           timevar="time",
                           idvar="id",
                           print=TRUE,
                           printlong=FALSE,
                           add=FALSE,
                           addtype="id",
                           varnames=c("LEFT","GAP","RIGHT")) {

  # Check dtms
  dtms_proper(dtms)

  # Sort data
  dataorder <- order(data[,idvar],
                     data[,timevar])
  data <- data[dataorder,]

  # Check for left censoring
  left <- by(data[,timevar],
             data[,idvar],
             FUN=function(x) min(x)>min(dtms$timescale))

  # Check for gaps
  gap <- by(data[,timevar],
            data[,idvar],
            FUN=function(x) any(diff(x)!=dtms$timestep))

  # Check for right censoring
  right1 <- by(data[,tovar],
              data[,idvar],
              FUN=function(x) !x[length(x)]%in%dtms$absorbing)

  right2 <- by(data[,timevar],
               data[,idvar],
               FUN=function(x) !max(x)==max(dtms$timescale))

  # Vectorize and combine
  left <- as.logical(left)
  gap <- as.logical(gap)
  right <- as.logical(right1) & as.logical(right2)

  # Print
  if(print) {
    cat("Units with left censoring: ",sum(left),"\n")
    cat("Units with gaps: ",sum(gap),"\n")
    cat("Units with right censoring: ",sum(right),"\n")
  }
  if(printlong) {
    cat("Cross tabulation:","\n")
    print(table(left,right,gap))
  }

  # Add indicators to data by ID
  if(add & addtype=="id") {

    # Frame for merging
    tmp <- data.frame(id=unique(data[,idvar]),
                      left=left,
                      gap=gap,
                      right=right)

    names(tmp) <- c(idvar,varnames)

    # Merge with data
    data <- merge(data,tmp)

    # Return
    return(data)

  }

  # Add indicators to data by observation
  if(add & addtype=="obs") {

    # Left censoring
    left <- by(data[,timevar],
               data[,idvar],
               FUN=function(x) {
                 first <- x[1]>min(dtms$timescale)
                 rest <- rep(FALSE,length(x)-1)
                 return(c(first,rest))
                 })

    # GGaps
    gap <- by(data[,timevar],
              data[,idvar],
              FUN=function(x) c(diff(x)!=dtms$timestep,FALSE))

    # Right censoring: Last state is not absorbing
    right1 <- by(data[,tovar],
                 data[,idvar],
                 FUN=function(x) {
                   last <- !x[length(x)]%in%dtms$absorbing
                   rest <- rep(FALSE,length(x)-1)
                   return <- c(rest,last)
                   })

    # Right censoring: Last obs is not at end of timescale
    right2 <- by(data[,timevar],
                 data[,idvar],
                 FUN=function(x) {
                   last <- !max(x)==max(dtms$timescale)
                   rest <- rep(FALSE,length(x)-1)
                   return <- c(rest,last)
                   })

    # Vectorize and combine
    left <- unlist(left)
    gap <- unlist(gap)
    right <- unlist(right1) & unlist(right2)

    # Put into data
    data[,varnames[1]] <- left
    data[,varnames[2]] <- gap
    data[,varnames[3]] <- right

    # Return
    return(data)

  }

}
#'  Cleans data in transition format
#'
#' @description
#' Cleans data in transition format. It can handle issues regularly occurring
#' with such data: transitions starting from or ending in missing states,
#' observations not covered by the time range, transitions starting or ending in
#' a state which is not in the state space, and observations starting in
#' absorbing states.
#'
#' @details
#' Transitions starting or ending with a missing state often occur for three
#' reasons. First, the function \code{dtms_format} will always create a transition
#' with a missing receiving state for the last observation of a unit, whether
#' due to censoring or not. For instance, if t=20 is the last value of the
#' time scale, and a unit is in state A at that time, then there will be a
#' transition starting at time t=20 from state A, and with receiving
#' state missing. Such transitions can usually be safely ignored, in particular
#' if there is only one absorbing state. Second, if, say, for a unit the last
#' observation is at time t=10 and censored after, there will be a transition
#' starting at time t=10 with missing receiving state. Whether such transitions
#' can be ignored depends on the censoring mechanism. If it is uninformative
#' these transitions can be dropped. Third, there might be missing values in
#' the sequence of states. For instance, a unit might first be in state A, then
#' state B, then the state is missing, and then state is again A, giving the
#' sequence A, B, NA, A. This implies a transition from B to NA, and from NA to
#' A.
#'
#' Transitions which are out of the time range can occur, for instance, when
#' the researcher is interested in a shorter time frame than covered by data
#' collection. In a clinical trial, the time scale could capture follow-up
#' time since start of the trial in months and data might be available for 60
#' months. But perhaps the researcher is only interested in the first 36 months.
#'
#' Transitions which start or end in a state which is not in the state space
#' occur when the states in the transition data are not included in the `dtms`
#' object. This likely will apply to states which rarely occur and which the
#' researcher does not want to combine with other states.
#'
#' @param data Data frame, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable with starting state. Default is "from".
#' @param tovar Character (optional), name of variable with receiving state. Default is "to".
#' @param timevar Character (optional), name of variable with time scale. Default is "time".
#' @param dropState Logical (optional), drop transitions with states which are not part of the state space. Default is TRUE.
#' @param dropTime Logical (optional), drop transitions with values of time not covered by the model. Default is TRUE.
#' @param dropNA Logical (optional), drop transitions with gaps, last observations, and similar. Default is TRUE.
#' @param dropAbs Logical (optional), drop transitions starting from absorbing states. Default is TRUE.
#' @param verbose Logical (optional), print how many transitions were dropped. Default is TRUE
#'
#' @return Cleaned data frame in transition format.
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
                       dtms,
                       fromvar="from",
                       tovar="to",
                       timevar="time",
                       dropTime=TRUE,
                       dropState=TRUE,
                       dropNA=TRUE,
                       dropAbs=TRUE,
                       verbose=TRUE) {

  # Check
  dtms_proper(dtms)

  # Drop observations not in state space
  if(dropState) {
    allstates <- c(dtms$transient,dtms$absorbing,NA)
    whichrows <- unlist(data[,fromvar])%in%allstates & unlist(data[,tovar])%in%allstates
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows not in state space\n")
    }
  }

  # Drop observations not in time range
  if(dropTime) {
    whichrows <- unlist(data[,timevar])%in%dtms$timescale
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows not in time range\n")
    }
  }

  # Drop missing values
  if(dropNA) {
    whichrows <- !is.na(unlist(data[,fromvar])) & !is.na(unlist(data[,tovar]))
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows starting or ending in NA\n")
    }
  }

  # Drop transitions starting in absorbing states
  if(dropAbs) {
    whichrows <- !unlist(data[,fromvar])%in%dtms$absorbing
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows starting in absorbing state\n")
    }
  }

  # Return
  return(data)

}
#' simpledata: an artificial dataset with abstract trajectories
#'
#' An artificial dataset with abstract states and time scale. The
#' state space consists of two transient states (A,B) and one absorbing state (X).
#'
#' @format `simpledata`
#' A data frame with 12,179 rows and 3 columns:
#' \describe{
#'   \item{id}{Identifier of the units}
#'   \item{time}{Time scale}
#'   \item{state}{The state occupied by an unit at a given age}
#' }
"simpledata"

#' workdata: simulated working trajectories
#'
#' A simulated dataset of individuals' working trajectories during late working life and retirement
#' age. The state space consists of three transient states (working; retired; not working)
#' and one absorbing state (dead). The age range covers ages 50 to 99. The data is
#' simulated using transition probabilities published as part of Dudel & Myrskylä (2017).
#'
#' @format `workdata`
#' A data frame with 250,000 rows and 4 columns:
#' \describe{
#'   \item{ID}{Person identifier}
#'   \item{Gender}{Individuals' gender (0=men, 1=women)}
#'   \item{Age}{Age, the time scale of this example}
#'   \item{State}{The state occupied by an unit at a given age}
#' }
#' @source <https://doi.org/10.1007/s13524-017-0619-6>
"workdata"

#' Summarize data in transition format
#'
#' @description
#' Returns a data frame with number of observed transitions (column COUNT),
#' relative proportion of a transition relative to all transitions (column
#' PROP), and raw transition probabilities Pr(j|i) (column PROB).
#'
#' @param data Data frame, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable  in `data`with receiving state. Default is "to".
#' @param weights Character (optional), name of variablein `data` with weights. Default is NULL.
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
#' Calculate delta
#'
#' @description
#' Calculates delta, either to compare transition probabilities from two
#' different models, or to assess how including lags changes transition
#' probabilities.
#'
#' @details
#' Delta is the weighted average absolute difference between the predicted
#' transition probabilities from two multistate models. It can attain values
#' between 0 and 1, where 0 indicates perfect similarity and 1 indicates that
#' the two models always give predictions at the opposite extremes; i.e., for
#' all predicted probabilities, one model predicts a probability of 0 and the
#' other predicts a probability of 1.
#'
#' This function is designed to use delta to assess the impact of including
#' different lags of the state variable in the model.
#'
#' To compare two different models, the arguments `data`, `model1`, and `model2`
#' are needed. `data` specifies the data frame used for predicting transition
#' probabilities. It needs to have all variables required for predicting based
#' on both `model1` and `model2`. The latter two arguments are the names of
#' multistate models estimated with \code{dtms_fit}.
#'
#' To compare how the inclusion of different lags of the state variable affects
#' predictions, a model needs to be specified using `data` and `dtms`, as well
#' as potential covariates with `controls`. The argument `lags` sets which lags
#' are included These are always including lower lags; e.g., a model including
#' the state at t-3 also has the state at t-2, at t-1, and at t. All resulting
#' models are compared to a model which does not control for the current or any
#' past state. If `lags=NULL` the Markov model is comapred to this model not
#' accounting for the current state or any past states.
#'
#' The argument `keepNA` controls how missing values are handled. These will
#' often occur for lagged states. For instance, for the first transition
#' observed for an individual, the state at time t is known, but not at time
#' t-1. In this case, if a first-order lag is used, this observation could either
#' be dropped; or, a missing value of the state at time t-1 could be included
#' as a predictor. `keepNA=TRUE` will do the latter, while if `FALSE`, all
#' observations with missing states are dropped. This is done for all
#' models, irrespective of the lag, such that they are based on exactly the
#' same observations.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param model1 Name of object containing a model estimated with \code{dtms_fit}.
#' @param model2 Name of object containing a model estimated with \code{dtms_fit}.
#' @param lags Numeric (optional), vector containing the lags as positive integers.
#' @param keepNA Logical (optional), keep missing values of lags as predictor value? Default is TRUE.
#' @param controls Character (optional), names of control variables
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable in `data` with receiving state. Default is "to".
#' @param idvar Character (optional), name of variable in `data` with unit ID. Default is "id".
#' @param timevar Character (optional), name of variable in `data` with time scale. Default is "time".
#' @param full Logical (optional), estimate fully interacted model? Default is FALSE.
#' @param package Character, chooses package for multinomial logistic regression, currently `VGAM`, `nnet`, and `mclogit` are supported. Default is `VGAM`.
#' @param reference Numeric or character (optional). Reference level of multinomial logistic regression.
#' @param ... Further arguments passed to estimation functions.
#'
#' @return Vector of values of delta
#' @export
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
#' ## Fit models
#' fit1 <- dtms_fit(data=estdata,controls="time")
#' fit2 <- dtms_fullfit(data=estdata,controls="time")
#'
#' ## Compare
#' dtms_delta(data=estdata,model1=fit1,model2=fit2)

dtms_delta <- function(data,
                       dtms=NULL,
                       model1=NULL,
                       model2=NULL,
                       lags=1:5,
                       controls=NULL,
                       fromvar="from",
                       tovar="to",
                       timevar="time",
                       idvar="id",
                       reference=1,
                       package="VGAM",
                       full=FALSE,
                       keepNA=TRUE,
                       ...) {

  # Compare models
  if(!is.null(model1) & !is.null(model2)) {

    # Predict based on model class
    if(inherits(model1,c("vgam","mclogit"))) {
      pred1 <- stats::predict(model1,data,"response")
    }

    if(inherits(model2,c("vgam","mclogit"))) {
      pred2 <- stats::predict(model2,data,"response")
    }

    if(inherits(model1,"nnet")) {
      pred1 <- stats::predict(model1,data,"probs")
    }

    if(inherits(model2,"nnet")) {
      pred2 <- stats::predict(model2,data,"probs")
    }

    # Calculate delta
    delta <- mean(rowSums(0.5*abs(pred1-pred2)))

    # Return
    return(delta)

  }

  # List for formulas
  formulist <- list()

  # Model without history
  formulist[[1]] <- dtms_formula(controls=controls,
                                 fromvar=NULL,
                                 tovar=tovar,
                                 full=full)

  # Markov model
  formulist[[2]] <- dtms_formula(controls=controls,
                                 fromvar=fromvar,
                                 tovar=tovar,
                                 full=full)

  # Loop over lags
  for(addlag in lags) {

    # Name of lag variable
    varname <- paste0(fromvar,"l",addlag)

    # Add to data
    data[,varname] <- dtms_lag(data=data,
                               dtms=dtms,
                               lag=addlag,
                               fromvar=fromvar,
                               idvar=idvar,
                               timevar=timevar)

    # Keep NA?
    if(keepNA) data[is.na(data[,varname]),varname] <- "NA"

    # As factor
    data[,varname] <- as.factor(data[,varname])

    # Add lag to controls
    controls <- c(controls,varname)

    # Generate formula
    nlag <- which(lags==addlag)
    formulist[[2+nlag]] <- dtms_formula(controls=controls,
                                      fromvar=fromvar,
                                      tovar=tovar,
                                      full=full)

  }

  # Keep NA?
  if(!keepNA) data <- stats::na.omit(data)

  # Number of models
  nmodels <- length(formulist)

  # Factors (needed by some packages)
  data[,fromvar] <- as.factor(data[,fromvar])
  data[,tovar] <- as.factor(data[,tovar])
  data[,tovar] <- stats::relevel(data[,tovar],ref=reference)

  # List for regression results
  fitlist <- list()

  # Making list assignment below simple
  count <- 1

  # Loop over models
  for(model in formulist) {

    # VGAM
    if(package=="VGAM") {

      # Estimate
      fitlist[[count]] <- VGAM::vgam(formula=model,
                                     family=VGAM::multinomial(refLevel=reference),
                                     data=data,
                                     ...)
    }

    #nnet
    if(package=="nnet") {

      # Estimate
      fitlist[[count]] <- nnet::multinom(formula=model,
                                         data=data,
                                         ...)
    }

    #mclogit
    if(package=="mclogit") {

      # Estimate
      fitlist[[count]] <- mclogit::mblogit(formula=model,
                                           data=data,
                                           ...)
    }

    count <- count + 1

  }

  # List for predicted results
  predictlist <- list()
  count <- 1

  # Predict
  for(fit in fitlist) {

    if(inherits(fit,c("vgam","mclogit"))) {
      predictlist[[count]] <- stats::predict(fit,data,"response")
    }

    if(inherits(fit,"nnet")) {
      predictlist[[count]] <- stats::predict(fit,data,"probs")
    }

    count <- count + 1

  }

  # Get differences
  resultlength <- nmodels-1
  result <- numeric(resultlength)

  for(i in 1:resultlength) {

    result[i] <- mean(rowSums(abs(predictlist[[1]] - predictlist[[i+1]])))

  }

  # Names
  names(result) <- c("Markov",paste("Lag",lags))

  # Return results
  return(result)

}
#' Summary for distributional results
#'
#' @description
#' This function provides several summary measures for results obtained
#' with dtms_visits, dtms_first, and dtms_last.
#'
#' @param distr An object of class `dtms_distr` created with \code{dtms_visits}, \code{dtms_first}, or \code{dtms_last}.
#'
#' @return A matrix with summary measures
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:19)
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
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' example <- dtms_visits(dtms=simple,
#'                        matrix=Tp,
#'                        risk="A")
#' summary(example)


dtms_distr_summary <- function(distr) {

  # Drop totals and similar
  allnames <- colnames(distr)
  allnames <- as.numeric(allnames)
  whichdrop <- is.na(allnames)
  distr <- distr[,!whichdrop]

  # Get values
  values <- colnames(distr)
  values <- as.numeric(values)

  # Mean
  MEAN <- apply(distr,1,function(x) sum(x*values))

  # Variance
  VARIANCE <- apply(distr,1,function(x) sum((values-sum(x*values))^2*x))

  # Standard deviation
  SD <- sqrt(VARIANCE)

  # Median
  MEDIAN <- apply(distr,1,function(x) values[min(which(cumsum(x)>0.5))])

  # Risk
  RISK0 <- apply(distr,1,function(x) x["0"])

  # Combine
  result <- cbind(MEAN,VARIANCE,SD,MEDIAN,RISK0)

  # Return
  return(result)

}
#' Generate variable with duration
#'
#' @description
#' This function creates a variable which measures the duration in a state.
#'
#' @details
#' Counting starts with 1 and the first occurence in a state. For instance,
#' if for an unit the sequence of states A, A, A, B, B, A, C is observed,
#' the duration variable would include 1, 2, 3, 1, 2, 1, 1.
#'
#' The argument `ignoreleft` controls how left censoring is handled; i.e., what
#' happens when for a unit there are no observations at the beginning of the
#' time scale. If `TRUE`, left censoring is ignored, and counting starts at
#' the first observation for a unit. For instance, if the time scale starts
#' at t=0, but the first observation for a unit is at time t=2, and the
#' sequence of states is again A, A, A, B, B, A, C, then `ignoreleft=TRUE`
#' returns 1, 2, 3, 1, 2, 1, 1. If `ignoreleft=FALSE`, then the function
#' would return NA, NA, NA, 1, 2, 1, 1.
#'
#' The function handles gaps in the data by setting the duration to NA. For
#' instance, if a unit is observed at times 1, 2, 4, 5, and 6, but not at time
#' 3, and the states are A, A, B, C, C, then the duration variable will
#' have the values 1, 2, NA, 1, 2.
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
#' dtms_duration(data=simpledata,
#'               dtms=simple)
dtms_duration <- function(data,
                          dtms,
                          newname="duration",
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
                function(x) dtms_duration_help(states=x[,1],
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
#' Calculate state expectancy
#'
#' @description
#' This function calculates the expected time spent in the transient states
#' (state expectancy).
#'
#' @details
#' If the argument `start_distr` is specified, the average of the state
#' expectancies over all starting states is calculated. The names and length
#' of `start_distr` need to match the starting states generated by this
#' function which are based on the `dtms` object.
#'
#' The partial expectancy for the time spent in the transient states can be
#' calculated using the arguments `start_time` and `end_time`.
#'
#' IF the argument `risk` is specified, then only the remaining life expectancy
#' for the state specified with this argument is shown, but for all time units
#' of the time scale.
#'
#' Two corrections to the results will be applied per default. Both corrections
#' are required as the underlying formulas do actually not provide the
#' expected time spent in a state, but the number of visits to a state. Time
#' and visits are only equal under certain conditions; in particular, only if
#' transitions between states happen mid-interval and the step length of the
#' time scale is equal to one. The first correction will remove a certain amount
#' of time spent in a certain state if its equal to the starting state. This is
#' controlled with the argument `correction` which is applied multiplicative.
#' For instance, its default value 0.5 means that the state expectancy in some
#' state X starting from state X is reduced by 0.5 time steps. The second
#' correction uses the entry `timestep` of `dtms`, and multiplies results with its value.
#'
#' @param matrix Matrix with transition probabilities, as generated with \code{dtms_matrix}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param risk Character (otpional), name of one transient state. If specified expectancies are only shown for this state but by values of the time scale.
#' @param start_distr Numeric (optional), distribution of starting states. If specified, average expectancy over all starting states will be calculated. Only applied if risk=NULL.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#' @param correction Numeric (optional), correction for expectancy when starting state and state under consideration match, see details. Defaults to 0.5.
#' @param total Logical (optional), calculate total expectancy. Default is TRUE. Only applied if risk=NULL.
#' @param verbose Logical (optional), print some information on what is computed. Default is FALSE.
#' @param fundamental Logical (optional), return fundamental matrix? Default is FALSE.
#'
#' @return A matrix with state expectancies.
#' @export
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
#' # Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Get transition matrix
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' ## Get starting distribution
#' S <- dtms_start(dtms=simple,
#'                 data=estdata)
#' ## State expectancies
#' dtms_expectancy(dtms=simple,
#'                 matrix=Tp,
#'                 start_distr=S)

dtms_expectancy <- function(matrix,
                            dtms,
                            risk=NULL,
                            start_distr=NULL,
                            start_time=NULL,
                            start_state=NULL,
                            end_time=NULL,
                            correction=0.5,
                            total=TRUE,
                            fundamental=FALSE,
                            verbose=FALSE) {

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=dtms$sep)

  # Number of starting and receiving states
  nstart <- length(starting)
  ntransient <- length(dtms$transient)

  # Remove absorbing states
  matrix <- dtms_absorbing(matrix)

  # All transient states
  allstates <- rownames(matrix)

  # Fundamental matrix
  nstates <- dim(matrix)[1]
  Nmat <- solve(diag(1,nstates)-matrix)

  # Correction
  if(is.numeric(correction)) {

    # Adjust
    diag(Nmat) <- diag(Nmat) - correction

    # Output
    if(verbose) cat("(Applying correction)","\n\n")
  }

  # Only return fundamental matrix?
  if(fundamental) {
    return(Nmat)
  }

  # Variant 1: Expectation of all transient states
  if(is.null(risk)) {

    # Matrix for results
    result <- matrix(data=NA,ncol=ntransient,nrow=nstart)
    rownames(result) <- paste0("start:",starting)
    colnames(result) <- dtms$transient

    for(i in 1:ntransient) {

      # Get states
      selector <- dtms_in(allstates,dtms$transient[i],dtms$sep)

      # Use end_time if specified
      if(!is.null(end_time)) {
        times <- dtms_gettime(allstates,dtms$sep)
        times <- times<=end_time
        times[!is.logical(times)] <- F
        selector <- selector & times
      }

      # Calculate results and place
      if(nstart>1) tmp <- rowSums(Nmat[starting,selector]) else tmp <- sum(Nmat[starting,selector])

      # Place
      result[,dtms$transient[i]] <- tmp
    }

  }

  # Variant 2: Expectation in one state by time scale
  if(!is.null(risk)) {

    # Check
    if(length(risk)!=1) stop("Only one state allowed for 'risk'")

    # Get time right
    first <- which(dtms$timescale==start_time)
    if(is.null(end_time)) last <- length(dtms$timescale) else
      last <- which(dtms$timescale==end_time)
    times <- dtms$timescale[first:last]
    ntimes <- length(times)

    # Get right columns from fundamental matrix
    selector1 <- dtms_in(allstates,risk,dtms$ep)
    selector2 <- dtms_gettime(allstates,dtms$sep)%in%times
    selector <- selector1 & selector2

    # Get result
    tmp <- rowSums(Nmat[,selector])

    # Matrix with results
    result <- matrix(data=tmp,ncol=ntimes,nrow=nstart,byrow=T)
    rownames(result) <- paste0("start:",starting)
    colnames(result) <- paste(times)

  }

  # Calculate average if starting distribution is provided
  if(!is.null(start_distr) & is.null(risk)) {

    # Check if matching
    if(length(start_distr)!=dim(result)[1]) stop("Starting distribution too long or short")

    # Match to starting/row ordering of result
    start_distr <- start_distr[match(names(start_distr),starting)]

    # Calculate
    AVERAGE <- colSums(result*start_distr)

    # Put into matrix for results
    result <- rbind(result,AVERAGE)
  }

  # Add row totals
  if(total & is.null(risk)) {
    TOTAL <- rowSums(result)
    result <- cbind(result,TOTAL)
  }

  # Adjust for time step
  if(dtms$timestep!=1) {
    result <- result*dtms$timestep
    if(verbose) cat("Adjusting for step length","\n\n")
  }

  # Return result
  return(result)

}
#' Time needed to reach a subset of states for the first time
#'
#' @description
#' This function calculates the distribution of the time needed to reach a
#' subset of states for the first time.
#'
#' @details
#' The resulting distribution is conditional on ever reaching the subset of
#' states, as it is not defined if the set is never reached. If the
#' argument `rescale` is set to FALSE, the distribution will not sum to one but
#' to the lifetime risk of ever reaching the subset.
#'
#' The state(s) which count to the time are specified with the argument `risk`.
#' If several states are specified, the resulting distribution refers to the
#' lifetime spent in any of the specified states.
#'
#' In a discrete-time model, the time spent in a state depends on assumptions
#' about when transitions happen. Currently, this functions supports two
#' variants which can be specified with the argument `method`: mid-interval
#' transitions can be selected with the option `mid` and imply that transitions
#' happen at the middle of the time interval; and the option `end` assumes
#' that instead transitions happen at the end of the interval. In this latter
#' case the distribution of the time spent in a state is equivalent to the
#' number of visits to that state. The calculation takes the step length of
#' the time scale into account as specified by the `dtms` object. If the
#' step length is not one fixed value, the first entry of `dtms$timestep` will
#' be used.
#'
#' If a distribution of the starting states is provided with `start_distr` the
#' output table has two additional rows. One shows the distribution
#' unconditional on the starting state. The other shows the distribution
#' conditional on not starting in any state of the risk set.
#'
#' @param matrix Matrix with transition probabilities, as generated with \code{dtms_matrix}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param risk Character, name of state(s) for which risk is of interest.
#' @param start_distr Numeric (optional), distribution of starting states. If specified, average distribution over all starting states will be calculated.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#' @param method Character (optional), do transitions happen mid-interval (`mid`, default) or at the end of the interval (`end`), see details.
#' @param total Logical (optional), should total of distribution be shown? See details. Default is FALSE.
#' @param rescale Logical (optional), should distribution be rescaled to sum to 1? See details. Default is TRUE.
#'
#' @return A table of the distribution of the time needed to reach the subset of states
#' @export
#'
#' @seealso
#' \code{\link{dtms_distr_summary}} to help with summarizing the resulting distribution.
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
#' # Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Get transition matrix
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' ## Get starting distribution
#' S <- dtms_start(dtms=simple,
#'                 data=estdata)
#' ## First visit
#' dtms_first(dtms=simple,
#'            matrix=Tp,
#'            risk="A",
#'            start_distr=S)

dtms_first  <- function(matrix,
                        dtms,
                        risk,
                        start_time=NULL,
                        start_state=NULL,
                        start_distr=NULL,
                        end_time=NULL,
                        method="mid",
                        total=TRUE,
                        rescale=TRUE) {

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=dtms$sep)

  # Time scale: Only transitions starting up to T-1 relevant
  timescale_reduced <- dtms$timescale[-length(dtms$timescale)]

  # States of the transition matrix
  allstates <- rownames(matrix)
  nstates <- length(allstates)

  # Select subset
  selectorU <- dtms_in(allstates,risk,dtms$sep)

  # Generate matrix
  P_E <- matrix(data=0,nrow=nstates,ncol=nstates)
  P_E[!selectorU,!selectorU] <- matrix[!selectorU,!selectorU]

  # Generate t
  if(is.null(start_time)) t <- 0 else t <- which(start_time==timescale_reduced)-1

  # Generate max
  if(is.null(end_time)) maxtime <- length(timescale_reduced)-1 else maxtime <- which(end_time==timescale_reduced)-1

  # Generate W_t_0 and W_t_0.5, initial conditions
  results <- vector("list",2)
  names(results) <- c("W_0","W_0.5")
  results[["W_0"]] <- dtms_mtexp(matrix,t)
  if(t==0) {
    results[["W_0.5"]] <- matrix(data=0,ncol=nstates,nrow=nstates)
    diag(results[["W_0.5"]][!selectorU,!selectorU]) <- 1
  } else results[["W_0.5"]] <- dtms_mtexp(P_E,t)

  # Variables
  upcoming <- 1.5
  past_steps <- c(0,0.5)
  steps <- 1
  end <- 0

  # Loop to generate results
  totalsteps <- length(t:maxtime)
  for(i in 1:totalsteps) {
    past_steps <- c(past_steps,upcoming)
    tmp <- vector("list",1)
    tmp[[1]] <- t(dtms_mtexp(P_E,upcoming-0.5)%*%
                    matrix(data=1,nrow=nstates,ncol=nstates))*results[["W_0.5"]]
    names(tmp) <- paste("W",upcoming,sep="_")
    results <- c(results,tmp)
    upcoming <- upcoming+1
  }

  # Generate V
  results_V <- vector("list",length(past_steps)-1)
  for(i in 1:length(results_V)) {
    names(results_V)[i] <- paste("V",past_steps[i],sep="_")
    results_V[[paste("V",past_steps[i],sep="_")]] <- results[[paste("W",past_steps[i],sep="_")]]-
                                                     results[[paste("W",past_steps[i+1],sep="_")]]
    colnames(results_V[[paste("V",past_steps[i],sep="_")]]) <- allstates
  }

  # Distribution
  tmp <- unlist(lapply(results_V, function(y) colSums(y)[starting]))

  # Result object
  result <- matrix(data=tmp,ncol=length(past_steps)-1,nrow=length(starting))
  rownames(result) <- paste0("start:",starting)

  # Get times right (interval length, mid interval vs end of interval)
  steps <- past_steps*dtms$timestep[1]
  if(method=="end") steps[-1] <- steps[-1]+0.5*dtms$timestep[1]
  colnames(result) <- paste(steps)[-length(steps)]

  # Average
  if(!is.null(start_distr)) {

    # Overall average
    AVERAGE <- start_distr%*%result
    rownames(AVERAGE) <- "AVERAGE"
    result <- rbind(result,AVERAGE)

    # Conditional on not starting in state in risk set
    whererisk <- !dtms$transient%in%risk
    tmp_distr <- start_distr[whererisk]/sum(start_distr[whererisk])
    tmp_names <- paste0("start:",names(tmp_distr))
    tmp <- tmp_distr%*%result[tmp_names,]
    result <- rbind(result,tmp)
    rownames(result)[dim(result)[1]] <- "AVERAGE(COND.)"

  }

  # Rescale
  if(rescale) {
    result <- t(apply(result,1,function(x) x/sum(x)))
  }

  # Total
  if(total) {
    TOTAL <- rowSums(result)
    result <- cbind(result,TOTAL)
    if(rescale) colnames(result)[dim(result)[2]] <- "TOTAL(RESCALED)"
  }

  # Assign class
  class(result) <- c("dtms_distr","matrix")

  # Output
  return(result)

} # End of function
#' Estimate (un)constrained discrete-time multistate model
#'
#' @description
#' This function estimates a (un)constrained discrete-time multistate model
#' using multinomial logistic regression.
#'
#' @details
#' The argument `data` takes a data set in transition format. The model formula
#' can either be specified by using the argument `formula`. Alternatively, it
#' can be specified with the arguments `fromvar`, `tovar`, and `controls`. These
#' are used if `formula` is not specified. `fromvar` takes the name of the
#' variable with the starting state as a character string, `tovar` the same for
#' the receiving state, and `controls` is an optional vector of control
#' variables. `fromvar` and `tovar` have default values which match other
#' functions of this package, making them a convenient alternative to `formula`
#' (see example).
#'
#' If `full=TRUE` a fully interacted model will be estimated in which each
#' control variable is interacted with all starting states. This is equivalent
#' to a full or unconstrained multistate model in which each transition is a
#' regression equation.
#'
#' The argument `package` is used choose the package used for estimation.
#' Currently, `VGAM` (default), `nnet`, and `mclogit` are supported.
#' The functions #' used for estimation are, respectively, `vgam`, `multinom`,
#' and `mblogit`. Arguments for these functions are passed via `...`.
#'
#' The argument `reference` sets the reference category for the multinomial
#' logistic regression. Weights for the regression can be passed via the
#' arguments `weights`. See the documentation of the package and function
#' used for estimation for details.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param controls Character (optional), names of control variables
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable in `data` with receiving state. Default is "to".
#' @param formula Formula (optional). If no formula is specified, it will be build from the information specified with controls, fromvar, tovar, and timevar.
#' @param full Logical (optional), estimate fully interacted model? Default is FALSE.
#' @param package Character, chooses package for multinomial logistic regression, currently `VGAM`, `nnet`, and `mclogit` are supported. Default is `VGAM`.
#' @param reference Numeric or character (optional). Reference level of multinomial logistic regression.
#' @param weights Character (optional). Name of variable with survey weights.
#' @param ... Further arguments passed to estimation functions.
#'
#' @return Returns an object with class depending on the package used.
#' @export
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
#' ## Fit model
#' fit <- dtms_fit(data=estdata)

dtms_fit <- function(data,
                     controls=NULL,
                     formula=NULL,
                     weights=NULL,
                     fromvar="from",
                     tovar="to",
                     reference=1,
                     package="VGAM",
                     full=FALSE,
                     ...) {

  # Require package used for estimation (requireNamespace does not help here)
  require(package,character.only=TRUE,quietly=TRUE)

  # Build formula if not specified
  if(is.null(formula)) formula <- dtms_formula(controls=controls,
                                               fromvar=fromvar,
                                               tovar=tovar,
                                               full=full)

  # Make sure environment for formula is correct (ugh)
  environment(formula) <- environment()

  # Get weights if specified
  if(!is.null(weights)) weights <- data[,weights]

  # Factors (needed by most packages)
  data[,fromvar] <- as.factor(data[,fromvar])
  data[,tovar] <- as.factor(data[,tovar])
  data[,tovar] <- stats::relevel(data[,tovar],ref=reference)

  # VGAM
  if(package=="VGAM") {

    # Estimate
    fitted <- VGAM::vgam(formula=formula,
                         family=VGAM::multinomial(refLevel=reference),
                         data=data,
                         weights=weights,
                         ...)
  }

  #nnet
  if(package=="nnet") {

    # Estimate
    fitted <- nnet::multinom(formula=formula,
                             data=data,
                             weights=weights,
                             ...)

  }

  #mclogit
  if(package=="mclogit") {

    # Estimate
    fitted <- mclogit::mblogit(formula=formula,
                               data=data,
                               weights=weights,
                               ...)

  }


  # Return results
  return(fitted)

}
#' Reshape data to transition format
#'
#' @description
#' Takes a data frame in long format and reshapes it into transition format.
#'
#' @details
#' The data frame supplied with the `data` argument has to be in long format,
#' where X is a time-constant covariate and Z(t) is a time-dependent covariate:
#'
#' \tabular{lllll}{
#' idvar \tab timevar \tab statevar \tab X \tab Z(t) \cr
#' 1 \tab 0 \tab A \tab x_1 \tab z_1(0) \cr
#' 1 \tab 1 \tab A \tab x_1 \tab z_1(1) \cr
#' 1 \tab 2 \tab B \tab x_1 \tab z_1(2) \cr
#' 1 \tab 3 \tab A \tab x_1 \tab z_1(3) \cr
#' 2 \tab 0 \tab B \tab x_2 \tab z_2(0) \cr
#' 2 \tab 1 \tab A \tab x_2 \tab z_2(1) \cr
#' ... \tab ... \tab ... \tab ... \tab ...
#' }
#'
#' If it is not in long format it has to be reshaped. The state variable
#' should provide the states as character strings or numbers; factors are not
#' supported.
#'
#' `dtms_format` turns the data set above into a data frame in transition
#' format:
#'
#' \tabular{llllll}{
#' id \tab time \tab fromvar \tab tovar \tab X \tab Z(t) \cr
#' 1 \tab 0 \tab A \tab A \tab x_1 \tab z_1(0) \cr
#' 1 \tab 1 \tab A \tab B \tab x_1 \tab z_1(1)  \cr
#' 1 \tab 2 \tab B \tab A \tab x_1 \tab z_1(2)  \cr
#' 2 \tab 0 \tab B \tab A \tab x_2 \tab z_2(0)  \cr
#' ... \tab ... \tab ... \tab ... \tab ... \tab ...
#' }
#'
#' Covariates do not need to be specified and are handled implicitly. The
#' transition from time t to t+1 takes covariate values from time t. By default
#' the variable names of the ID variable and the time variable are changed to
#' `id` and `time`, as the other functions of the package use these as default
#' names. If renaming the variables is not possible because these variable
#' names already appear in the data then the original names are used. If
#' `keepnames=TRUE` the original names for `id` and `time` are kept.
#'
#' `dtms_format` by default drops gaps in the data, as no transitions are
#' observed. For instance, in the following example there is no observation at
#' time 4, and thus no transition is observed from t=3 to t=4; and no
#' transition from t=4 to t=5:
#'
#' \tabular{lll}{
#' idvar \tab timevar \tab statevar \cr
#' 1 \tab 0 \tab A \cr
#' 1 \tab 1 \tab A \cr
#' 1 \tab 2 \tab B \cr
#' 1 \tab 3 \tab A \cr
#' 1 \tab 5 \tab A \cr
#' }
#'
#' In this example, `dtms_format` will return the following:
#'
#'  \tabular{llll}{
#' id \tab time \tab from \tab to \cr
#' 1 \tab 0 \tab A \tab A \cr
#' 1 \tab 1 \tab A \tab B \cr
#' 1 \tab 2 \tab B \tab A \cr
#' }
#'
#' If `fill=T`, then `dtms_format` will return the following:
#'
#'  \tabular{llll}{
#' id \tab time \tab from \tab to \cr
#' 1 \tab 0 \tab A \tab A \cr
#' 1 \tab 1 \tab A \tab B \cr
#' 1 \tab 2 \tab B \tab A \cr
#' 1 \tab 3 \tab A \tab NA \cr
#' 1 \tab 4 \tab NA \tab A \cr
#' }
#'
#' The argument \code{absorbing} controls if the first observed absorbing state
#' is carried over to later observations. For instance, if a unit is first
#' observed to be in transient state 'A' at time 1, then in absorbing state 'X'
#' at time 2, and then in transient state 'A' at time 3, \code{absorbing=TRUE}
#' will lead to replacement of the state at time 3 with 'X'.
#'
#' @param data Data frame in long format.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param idvar Character (optional), name of variable in `data` with unit ID. Default is "id".
#' @param timevar Character (optional), name of variable in `data` with time scale. Default is "time".
#' @param statevar Character (optional), name of variable in `data` with state, default is 'state'.
#' @param fromvar Character (optional), name of variable`with starting state in reshaped data. Default is "from".
#' @param tovar Character (optional), name of variable with receiving state in reshaped data. Default is "to".
#' @param absorbing Logical (optional), use first observed absorbing state consistently? See details. Default is TRUE.
#' @param keepnames Logical (optional), keep original names for id and time variable? Default is FALSE; i.e., not keeping original names.
#' @param fill Logical (optional), fill implicit missing values with explicit NA? Default is FALSE.
#' @param verbose Logical (optional), create output to console if changing variable names is not possible? Default is TRUE.
#' @param steplength Logical (optional), if true, the time to the next state is returned as a variable. Default is FALSE.
#' @param stepvar Character (optional), if \code{steplength==TRUE}, this specifies the name of the variable with the step length. Default is `steplength`.
#'
#' @return A data set reshaped to transition format
#' @export
#'
#' @seealso
#' \code{\link{dtms_data_summary}} to summarize data in transition format.
#' \code{\link{dtms_censoring}} for descriptives on censoring.
#' \code{\link{dtms_clean}} for fast data cleaning.
#'
#' @examples
#' # Define model
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:20)
#' # Transiton format
#' estdata <- dtms_format(data=simpledata,
#' dtms=simple,
#' idvar="id",
#' timevar="time",
#' statevar="state")

dtms_format <- function(data,
                        dtms,
                        idvar="id",
                        timevar="time",
                        statevar="state",
                        fromvar="from",
                        tovar="to",
                        absorbing=TRUE,
                        keepnames=FALSE,
                        fill=FALSE,
                        verbose=TRUE,
                        steplength=FALSE,
                        stepvar="steplength") {

  # Transform to data frame, e.g., if tibble
  if(class(data)[1]!="data.frame") data <- as.data.frame(data)

  # Check
  dtms_proper(dtms)

  # Fill data
  if(fill) {

    # Get ID values
    idvalues <- data[,idvar] |> unique()

    # Full data
    fulldata <- expand.grid(dtms$timescale,idvalues,
                            stringsAsFactors=FALSE)
    names(fulldata) <- c(timevar,idvar)

    # Merge with data
    data <- merge(fulldata,data,
                  by=c(idvar,timevar),
                  all=T)

    # Drop temporary data
    rm(fulldata)

  }

  # Sort data
  dataorder <- order(data[,idvar],
                     data[,timevar])
  data <- data[dataorder,]

  # Check if next value is valid
  consecutive <- dtms_consecutive(data=data,
                                  idvar=idvar,
                                  timevar=timevar,
                                  timestep=dtms$timestep)

  # Absorbing states carry forward
  if(absorbing) {
    tmp <- tapply(data[,statevar],data[,idvar],function(x) dtms_carry(x=x,dtms=dtms))
    data[,statevar] <- unlist(tmp)
  }

  # Get next state
  data[,tovar] <- NA
  tovalues <- c(data[-1,statevar],NA)
  data[consecutive$true,tovar] <- tovalues[consecutive$true]
  if(steplength) data[consecutive$true,stepvar] <- consecutive$numeric[consecutive$true]

  # Rename from variable
  data <- dtms_rename(data,statevar,fromvar)

  # Change names of id and time if possible
  if(!keepnames) {
    if(timevar!='time' & !'time'%in%names(data))
      data <- dtms_rename(data,timevar,"time") else
      if(verbose) cat("Kept original name for time \n")
    if(idvar!='id' & !'id'%in%names(data))
      data <- dtms_rename(data,idvar,"id") else
      if(verbose) cat("Kept original name for id \n")
  }

  # Class
  class(data) <- c("dtms_data","data.frame")

  # Return result
  return(data)

}
#' Carry states forward
#'
#' @description
#' This function carries a state forward after its first occurrence.
#'
#' @details
#' This function carries a state forward after its first occurrence.
#' For instance, carrying the state "A" forward in the sequence `B, B, A, B, B`
#' will give the sequence `B, B, A, A, A`. The sequence `C, B, C, A, B, A, A, B`
#' will give `C, B, C, A, A, A, A, A`.
#'
#' This function works with data frames in transition format and in long format.
#' The default is transition format, using the arguments `fromvar` and `tovar`.
#' If, however, the argument `statevar` is specified, it is used instead.
#'
#' The argument `overwrite` is used to control what type of information is
#' replaced. If `overwrite==transient`, then only transient states are replaced
#' while missing values and absorbing states remain unchanged. For example,
#' carrying forward state "A" in the sequence `B, B, A, B, NA, X, X` with X
#' being an absorbing state will give `B, B, A, A, NA, X, X`. If
#' `overwrite==missing` then in addition to transient states also missing values
#' are replaced and for the example sequence `B, B, A, A, A, X, X` would be
#' returned. If `overwrite==absorbing` then in addition to transient states
#' absorbing states will be replaced; for the example sequence the result would
#' be `B, B, A, A, NA, A, A`. Finally, if `overwrite==all` then all values in
#' the sequence will be replaced: `B, B, A, A, A, A, A`.
#'
#' @seealso
#' \code{\link{dtms_backward}} to carry states backward.
#'
#' @param data A data frame in long format.
#' @param state Character, name of the state to be carried forward.
#' @param fromvar Character (optional), name of variable with starting state. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state. Default is `to`.
#' @param statevar Character (optional), name of the variable in the data frame in long format with the states. Default is NULL.
#' @param idvar Character (optional), name of variable with unit ID. Default is `id`.
#' @param timevar Character (optional), name of variable with time scale. Default is `time`.
#' @param dtms dtms object (optional), as created with \code{dtms}. Not required if `overwrite==transient`.
#' @param overwrite Character (optional), one of `transient`, `missing`, `absorbing`, and `all`, see details. Default is `transient`.
#' @param vector Logical (optional), return vector (if TRUE) or data frame (if FALSE). Default is FALSE. Argument is only used if argument `statevar` is specified.
#'
#' @return The data frame specified with `data` and the edited state variable (if `vector=FALSE`) or a vector (if `vector=TRUE`).
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:19)
#'
#' dtms_forward(data=simpledata,
#'              statevar="state",
#'              state="A",
#'              dtms=simple,
#'              overwrite="transient")

dtms_forward <- function(data,
                         state,
                         fromvar="from",
                         tovar="to",
                         statevar=NULL,
                         idvar="id",
                         timevar="time",
                         dtms=NULL,
                         overwrite="missing",
                         vector=FALSE) {

  # Check dtms only if specified
  if(!is.null(dtms)) dtms_proper(dtms)

  # Sort data
  dataorder <- order(data[,idvar],
                     data[,timevar])
  data <- data[dataorder,]

  # Apply helper for transition format
  if(is.null(statevar)) {

    # Move starting state
    data[,tovar][data[,fromvar]==state] <- state

    # Carry forward starting state
    tmp1 <- tapply(data[,fromvar],
                   data[,idvar],
                   function(x) dtms_forward_help(x=x,
                                                 state=state,
                                                 overwrite=overwrite,
                                                 dtms=dtms))
    # Carry forward receiving state
    tmp2 <- tapply(data[,tovar],
                   data[,idvar],
                   function(x) dtms_forward_help(x=x,
                                                 state=state,
                                                 overwrite=overwrite,
                                                 dtms=dtms))

    # Assign new values
    data[,fromvar] <- unlist(tmp1)
    data[,tovar] <- unlist(tmp2)


  } else {

    # Apply helper for long format
    tmp <- tapply(data[,statevar],
                  data[,idvar],
                  function(x) dtms_forward_help(x=x,
                                                state=state,
                                                overwrite=overwrite,
                                                dtms=dtms))

    # Return vector
    if(vector) return(unlist(tmp))

    # Assign new values
    data[,statevar] <- unlist(tmp)

  }

  # Return
  return(data)

}
#' Estimate unconstrained discrete-time multistate model
#'
#' @description
#' This function estimates an unconstrained discrete-time multistate model
#' using multinomial logistic regression. This is achieved by interacting
#' the starting state with all predictors in the model. It is a wrapper for
#' \code{dtms_fit()} with `full=TRUE` and otherwise slightly less arguments.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param controls Character (optional), names of control variables
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable in `data` with receiving state. Default is "to".
#' @param formula Formula (optional). If no formula is specified, it will be build from the information specified with controls, fromvar, tovar, and timevar.
#' @param weights Character (optional), name of variable in `data` with survey weights.
#' @param reference Numeric or character (optional). Reference level of multinomial logistic regression.
#' @param package Character, chooses package for multinomial logistic regression, currently `VGAM`, `nnet`, and `mclogit` are supported. Default is `VGAM`.
#' @param ... Further arguments passed to estimation functions.
#'
#' @return Returns an object with class depending on the package used.
#' @export
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
#' ## Fit model
#' fit <- dtms_fullfit(data=estdata)

dtms_fullfit <- function(data,
                         controls=NULL,
                         formula=NULL,
                         weights=NULL,
                         fromvar="from",
                         tovar="to",
                         reference=1,
                         package="VGAM",
                         ...) {

  dtms_fit(data=data,
           controls=controls,
           formula=formula,
           weights=weights,
           fromvar=fromvar,
           tovar=tovar,
           reference=reference,
           package=package,
           full=TRUE)

}
### Create transition matrix without absorbing states
dtms_absorbing <- function(matrix) { # matrix=full transition matrix

  # Get states which are absorbing
  to_remove <- which(diag(matrix)==1)

  # Remove from matrix
  removed <- matrix[-to_remove,-to_remove]

  # Output reduced matrix
  return(removed)

}

### Check if object is proper dtms object
dtms_proper <- function(dtms) { # dtms=object to be checked

  # Error message
  message <- "Not a proper dtms object."

  # Check class
  if(!class(dtms)[2]=="dtms") stop(message)

  # Check names
  listnames <- c("transient" ,"absorbing" ,"timescale" ,"timestep" ,"sep")
  if(!all(names(dtms)==listnames)) stop(message)

  # Check types (and length)
  if(!is.character(dtms$transient)&!is.numeric(dtms$transient)) stop(message)
  if(!is.character(dtms$absorbing)&!is.numeric(dtms$absorbing)) stop(message)
  if(!is.character(dtms$sep)) stop(message)
  if(!is.numeric(dtms$timescale) | length(dtms$timescale)<2) stop(message)
  if(!is.numeric(dtms$timestep)) stop(message)

}

### Check if values of time scale are consecutive
dtms_consecutive <- function(data,idvar,timevar,timestep) {

  # Make sure no missing times and ids
  if(any(is.na(data[,timevar]))) stop("Missing values in time variable not allowed")
  if(any(is.na(data[,idvar]))) stop("Missing values in ID variable not allowed")

  # Get diff to next time step
  consecutive <- by(data[,timevar],data[,idvar],FUN=diff)

  # Add last obs
  consecutive <- lapply(consecutive,function(x) c(x,-1))

  # Unlist
  consecutive <- unlist(consecutive)

  # TRUE if equal to timestep, FALSE otherwise
  result <- data.frame(true=consecutive%in%timestep,
                       numeric=consecutive)

  # Return
  return(result)

}

### Rename variables
dtms_rename <- function(data,oldnames,newnames) {

  # Which names to change?
  changenames <- match(oldnames,names(data))

  # Change names
  names(data)[changenames] <- newnames

  # Return
  return(data)

}

### Combines values; e.g., to combine state names with time scale values
dtms_combine <- function(values1,values2,sep) {

  # Generate output vector
  output <- character(0)

  # Get values
  for(value in values1) {
    output <- c(output,paste(value,values2,sep=sep))
  }

  # Return
  return(output)

}

### Check if a short state name is in long name
dtms_in <- function(vector,name,sep) {
  res <- lapply(strsplit(vector,split=sep),function(y) any(name%in%y[[1]]) )
  res <- unlist(res)
  return(res)
}

### Get time from long name
dtms_gettime <- function(vector,sep) {
  res <- lapply(strsplit(vector,split=sep),function(y) y[2] )
  res <- unlist(res) |> as.numeric()
  return(res)
}

### Get state from long name
dtms_getstate <- function(vector,sep) {
  res <-lapply(strsplit(vector,split=sep),function(x) x[1])
  res <- unlist(res)
  return(res)
}

### Calculate power of matrix
dtms_mtexp <- function(matrix,n) {

  res <- diag(nrow=nrow(matrix))
  rep <- 0

  while(rep<n) {
    res <- res %*% matrix
    rep <- rep+1
  }

  return(res)

}

### Generate formula for estimation
dtms_formula <- function(controls, # Arguments the same as for dtms_fit
                         fromvar,
                         tovar,
                         full) {

  # If fromvar is NULL (not default, needs explicit call)
  if(is.null(fromvar)) fromvar <- "1"

  # Constrained model
  if(!full) {
    formula <- paste0(tovar,"~",fromvar)
    if(!is.null(controls)) {
      controls <- paste(controls,collapse="+")
      formula <- paste(formula,controls,sep="+")
    }
    formula <- stats::as.formula(formula)
    # Unconstrained/fully interacted model
  } else {
    formula <- paste0(tovar,"~",fromvar)
    if(!is.null(controls)) {
      varlist <- paste(controls,fromvar,sep="*")
      varlist <- paste(varlist,collapse="+")
      formula <- paste(formula,varlist,sep="+")
    }
    formula <- stats::as.formula(formula)
  }

  # Return
  return(formula)

}

### Get lagged state variable, accounting for gaps in the data
dtms_lag <- function(data,
                     dtms,
                     lag,
                     fromvar="from",
                     idvar="id",
                     timevar="time") {

  # Make data smaller
  data <- data[,c(idvar,timevar,fromvar)]

  # Get ID values
  idvalues <- data[,idvar] |> unique()

  # Full data
  fulldata <- expand.grid(dtms$timescale,idvalues,
                          stringsAsFactors=FALSE)
  names(fulldata) <- c(timevar,idvar)

  # Merge with data
  fulldata <- merge(fulldata,data,
                    by=c(idvar,timevar),
                    all=T)

  # shift state
  stateshift <- by(fulldata[,fromvar],
                   fulldata[,idvar],function(x) c(rep(NA,min(lag,length(x))),
                                                  x[-( max(0,(length(x)-(lag-1))): length(x))]))

  # shift time
  timeshift <- by(fulldata[,timevar],
                  fulldata[,idvar],function(x) c(rep(NA,min(lag,length(x))),
                                                 diff(x,lag=lag) ))

  # Unlist
  stateshift <- unlist(stateshift)
  timeshift <- unlist(timeshift)

  # drop wrong spacing
  stateshift[!timeshift%in%c(NA,dtms$timestep*lag)] <- NA

  # Merge back
  fulldata$stateshift <- stateshift
  data <- merge(data,fulldata,by=c(idvar,timevar,fromvar),sort=FALSE)

  # Return
  return(data$stateshift)

}

### Carry over absorbing values
dtms_carry <- function(x,
                       dtms) {

  # Check if action necessary
  if(any(x%in%dtms$absorbing)) {

    # Length
    n <- length(x)

    # From where to carry over
    whichfirst <- which(x%in%dtms$absorbing)[1]

    # What to carry over
    whichvalue <- x[whichfirst]

    # Carry over
    x[whichfirst:n] <- whichvalue

    # Return
    return(x)

  } else return(x)

}

### Move states forward
dtms_forward_help <- function(x, # Vector of states
                              state, # State name as character string
                              overwrite="missing", # transient, missing, absorbing, all
                              dtms=NULL) {

  # Length
  nx <- length(x)

  # Find appearances
  whichfirst <- which(x==state)

  # Stop if no appearance
  if(length(whichfirst)==0) return(x)

  # Get first
  whichfirst <- min(whichfirst)

  # Replace: both missing and absorbing
  if(overwrite=="all") {
    dochange <- whichfirst:nx
    x[dochange] <- state
  }

  # Replace: only missing
  if(overwrite=="missing") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    dochange <- whichfirst:nx
    dochange <- setdiff(dochange,whichabsorbing)
    x[dochange] <- state
  }

  # Replace: only absorbing
  if(overwrite=="absorbing") {
    whichmissing <- which(is.na(x))
    dochange <- whichfirst:nx
    dochange <- setdiff(dochange,whichmissing)
    x[dochange] <- state
  }

  # Do not replace absorbing or missing
  if(overwrite=="transient") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    whichmissing <- which(is.na(x))
    whichboth <- union(whichabsorbing,whichmissing)
    dochange <- whichfirst:nx
    dochange <- setdiff(dochange,whichboth)
    x[dochange] <- state
  }

  # Return
  return(x)

}

### Carry states backward
dtms_backward_help <- function(x, # Vector of states
                               state, # State name as character string
                               overwrite="missing", # transient, missing, absorbing, all
                               dtms=NULL) {

  # Find appearances
  whichfirst <- which(x==state)

  # Stop if no appearance
  if(length(whichfirst)==0) return(x)

  # Get last
  whichfirst <- max(whichfirst)

  # Replace: both missing and absorbing
  if(overwrite=="all") {
    dochange <- 1:whichfirst
    x[dochange] <- state
  }

  # Replace: only missing
  if(overwrite=="missing") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    dochange <- 1:whichfirst
    dochange <- setdiff(dochange,whichabsorbing)
    x[dochange] <- state
  }

  # Replace: only absorbing
  if(overwrite=="absorbing") {
    whichmissing <- which(is.na(x))
    dochange <- 1:whichfirst
    dochange <- setdiff(dochange,whichmissing)
    x[dochange] <- state
  }

  # Do not replace absorbing or missing
  if(overwrite=="transient") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    whichmissing <- which(is.na(x))
    whichboth <- union(whichabsorbing,whichmissing)
    dochange <- 1:whichfirst
    dochange <- setdiff(dochange,whichboth)
    x[dochange] <- state
  }

  # Return
  return(x)

}

### Duration helper
dtms_duration_help <- function(states, # Vector of states of one unit
                               time, # Vector of time scale values
                               dtms, # dtms object
                               ignoreleft) { # TRUE or FALSE, as per dtms_duration

  # Get lengths
  lengths <- rle(states)
  duration <- unlist(lapply(lengths$length, function(x) 1:x))

  # Change if left censoring
  if(!ignoreleft & time[1]!=dtms$timescale[1]) duration[1:lengths[1]] <- NA

  # Remove gaps
  timediffs <- diff(time)
  timediffs <- c(FALSE,!timediffs%in%dtms$timestep)

  # Go through all gaps
  if(any(timediffs)) {

    # Find entries
    dropwhich <- which(timediffs)
    cumlengths <- cumsum(lengths$length)
    cumwhere <- c(1,cumlengths[-length(cumlengths)]+1)

    # Loop and replace
    for(drop in dropwhich) {

      # Which entries need to be replaced
      whichlengths <- which(cumlengths>=drop)[c(1,2)]
      from <- drop
      to <- cumwhere[whichlengths[2]]-1
      if(is.na(to)) to <- cumwhere[whichlengths[1]]
      # Replace
      duration[from:to] <- NA
    }
  }

  # Return
  return(duration)

}


## Occurrence helper
dtms_occurrence_help <- function(states, # Vector of states of one unit
                                time, # Vector of time scale values
                                dtms, # dtms object
                                ignoreleft) { # TRUE or FALSE, as per dtms_occurence

  # Change if left censoring
  if(!ignoreleft & time[1]!=dtms$timescale[1]) {
    result <- rep(NA,length(states))
    return(result)
  }

  # Get states of spells
  spells <- rle(states)$values

  # Count occurence of spells
  if(any(is.na(spells))) occurences <- table(spells,useNA="always") else occurences <- table(spells)

  # Start counting from 1
  expanded_occurences <- lapply(occurences,function(x) 1:x)

  # Get names, vectorize
  names_expanded <- names(expanded_occurences)
  occurences <- unlist(expanded_occurences)

  # Match ordering
  ordering <- unlist(lapply(names_expanded,function(x) {
    if(is.na(x)) which(is.na(spells)) else which(x==spells)
  }))

  # Get counts in right order
  result <- numeric(length(occurences))
  result[ordering] <- occurences

  # Repeat counts for each observation in spell
  result <- inverse.rle(list(values=result,lengths=rle(states)$lengths))

  # Handle gaps: find if any
  timediffs <- diff(time)
  timediffs <- c(!timediffs%in%dtms$timestep)

  # If gap fill with NA starting from first gap
  if(any(timediffs)) {
      dropwhich <- min(which(timediffs))
      nresult <- length(result)
      result[(dropwhich+1):nresult] <- NA
    }

  # Return result
  return(result)

}
#' Calculate the distribution of the time until a subset of states is left for
#' the last time.
#'
#' @description
#' Calculates the distribution of the until a subset of states is left for the
#' very last time.
#'
#' @details
#' The resulting distribution is conditional on ever experiencing the final
#' exit, as the waiting time otherwise is not a finite number. The argument
#' `rescale` can be used to control whether the distribution is rescaled to
#' sum to 1; it usually will do without rescaling.
#'
#' The state(s) which count to the time are specified with the argument `risk`.
#' If several states are specified, the resulting distribution refers to the
#' lifetime spent in any of the specified states. The optional argument
#' `risk_to` can be used to restrict results to exits from the set `risk` to
#' another specific subset defined by `risk_to`; i.e., this way, not all
#' transitions out of `risk` count for the final exit, but only those to
#' specific states.
#'
#' In a discrete-time model, the time spent in a state depends on assumptions
#' about when transitions happen. Currently, this functions supports two
#' variants which can be specified with the argument `method`: mid-interval
#' transitions can be selected with the option `mid` and imply that transitions
#' happen at the middle of the time interval; and the option `end` assumes
#' that instead transitions happen at the end of the interval. In this latter
#' case the distribution of the time spent in a state is equivalent to the
#' number of visits to that state. The calculation takes the step length of
#' the time scale into account as specified by the `dtms` object. If the
#' step length is not one fixed value, the first entry of `dtms$timestep` will
#' be used.
#'
#' If a distribution of the starting states is provided with `start_distr` the
#' output table has two additional rows. One shows the distribution
#' unconditional on the starting state. The other shows the distribution
#' conditional on not starting in any state of the risk set.
#'
#' The distribution of partial waiting times can be generated using the arguments
#' `start_state` and `start_time` in combination with `end_time`.
#'
#' @param matrix Matrix with transition probabilities, as generated with \code{dtms_matrix}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param risk Character, name of state(s) for which risk is of interest.
#' @param risk_to Character (optional), names of one or several states to which the states specified in `risk` are left. See details.
#' @param start_distr Numeric (optional), distribution of starting states. If specified, average distribution over all starting states will be calculated.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#' @param method Character (optional), do transitions happen mid-interval (`mid`, default) or at the end of the interval (`end`), see details.
#' @param rescale Logical (optional), should distribution be rescaled to sum to 1? See details. Default is TRUE.
#' @param total Logical, should total of distribution be shown (always sums to 1)? Default is FALSE.
#'
#' @return Matrix with the distribution(s) of the waiting time.
#' @export
#'
#' @seealso
#' \code{\link{dtms_distr_summary}} to help with summarizing the resulting distribution.
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
#' # Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Get transition matrix
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' ## Get starting distribution
#' S <- dtms_start(dtms=simple,
#'                 data=estdata)
#' ## First visit
#' dtms_last(dtms=simple,
#'            matrix=Tp,
#'            risk="A",
#'            start_distr=S)

dtms_last <- function(matrix,
                      dtms,
                      risk,
                      risk_to=NULL,
                      start_time=NULL,
                      start_state=NULL,
                      start_distr=NULL,
                      end_time=NULL,
                      method="mid",
                      total=TRUE,
                      rescale=TRUE){

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=dtms$sep)

  # Time scale: Only transitions starting up to T-1 relevant
  timescale_reduced <- dtms$timescale[-length(dtms$timescale)]

  # States of the transition matrix
  allstates <- rownames(matrix)
  nstates <- length(allstates)

  # Select subset
  selectorD <- dtms_in(allstates,risk,dtms$sep)
  if(is.null(risk_to)) selectorU <- !selectorD else
    selectorU <- dtms_in(allstates,risk_to,dtms$sep)

  # Get maxtime
  if(is.null(end_time)) maxtime <- length(timescale_reduced)-1 else
    maxtime <- which(end_time==timescale_reduced)-1

  # Generate t
  if(is.null(start_time)) t <- 0 else
    t <- which(start_time==timescale_reduced)-1

  # Partition transition matrix
  P_E <- matrix(data=0,nrow=nstates,ncol=nstates)
  P_S <- matrix(data=0,nrow=nstates,ncol=nstates)
  P_E[selectorD,selectorU] <- matrix[selectorD,selectorU]
  P_S[!selectorD,] <- matrix[!selectorD,]
  P_S[selectorD,!selectorU] <- matrix[selectorD,!selectorU]

  # Special case E_0
  ones <- matrix(data=1,nrow=nstates,ncol=nstates)
  results <- vector("list",1)
  names(results) <- "E_0"
  results[["E_0"]] <- t(dtms_mtexp(P_S,maxtime-t+1)%*%ones) *
                          dtms_mtexp(matrix,t)
  colnames(results[["E_0"]]) <- allstates

  # Loop for other E_x
  step <- 0.5
  past.steps <- c(0)

  while(!step>=(maxtime-t+1)) {
    past.steps <- c(past.steps,step)
    e <- step-0.5
    tmp <- vector("list",1)
    names(tmp) <- paste("E",step,sep="_")
    tmp[[paste("E",step,sep="_")]] <- t(dtms_mtexp(matrix,e)%*%P_E%*%
                                        dtms_mtexp(P_S,maxtime-t-e)%*%ones) *
                                        dtms_mtexp(matrix,t)
    colnames(tmp[[paste("E",step,sep="_")]]) <- allstates
    results <- c(results,tmp)
    step <- step+1
  }

  # Get distribution
  tmp <- unlist(lapply(results, function(z) colSums(z)[starting]))

  # Result object
  result <- matrix(data=tmp,ncol=length(past.steps),nrow=length(starting))
  rownames(result) <- paste0("start:",starting)

  # Get times right (interval length, mid interval vs end of interval)
  steps <- past.steps*dtms$timestep[1]
  if(method=="end") steps[-1] <- steps[-1]+0.5*dtms$timestep[1]
  colnames(result) <- paste(steps)

  # Average
  if(!is.null(start_distr)) {

    # Overall average
    AVERAGE <- start_distr%*%result
    rownames(AVERAGE) <- "AVERAGE"
    result <- rbind(result,AVERAGE)

    # Conditional on not starting in state in risk set
    whererisk <- !dtms$transient%in%risk
    tmp_distr <- start_distr[whererisk]/sum(start_distr[whererisk])
    tmp_names <- paste0("start:",names(tmp_distr))
    tmp <- tmp_distr%*%result[tmp_names,]
    result <- rbind(result,tmp)
    rownames(result)[dim(result)[1]] <- "AVERAGE(COND.)"

  }

  # Rescale
  if(rescale) {
    result <- result[,-1]
    result <- t(apply(result,1,function(x) x/sum(x)))
  }

  # Total
  if(total) {
    TOTAL <- rowSums(result)
    result <- cbind(result,TOTAL)
    if(rescale) colnames(result)[dim(result)[2]] <- "TOTAL(RESCALED)"
  }

  # Assign class
  class(result) <- c("dtms_distr","matrix")

  # Output
  return(result)

}
#' Creates a transition matrix from transition probabilities
#'
#' @description
#' This function creates a transiton matrix based on transition probabilities
#' predicted using the function `dtms_transitions`.
#'
#' @param probs Data frame with transition probabilities, as created with \code{dtms_transitions}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable in `probs` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable in `probs` with receiving state. Default is "to".
#' @param Pvar Character (optional), name of variable in `probs` with transition probabilities. Default is `P`.
#' @param enforcedeath Logical (optional), make sure that every unit moves to absorbing state after last value of time scale? Default is TRUE.
#' @param rescale Logical (optional), rescale transition probabilities to sum to 1? Default is TRUE.
#' @param reshapesep Character (optional), used in re-arranging the transition probabilities; should not appear in any state name. Default is `:`.
#'
#' @return Returns a transition matrix.
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:20)
#' ## Reshape to transition format
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' ## Clean
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' ## Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Get transition matrix
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)

dtms_matrix <- function(probs,
                        dtms=NULL,
                        fromvar="from",
                        tovar="to",
                        Pvar="P",
                        enforcedeath=T,
                        rescale=T,
                        reshapesep=":") {

  # Check
  dtms_proper(dtms)

  # Combine states and time
  transient_states <- dtms_combine(dtms$transient,dtms$timescale,sep=dtms$sep)
  absorbing <- paste(dtms$absorbing)
  all_states <- c(transient_states,absorbing)

  # Get variable names in probs right
  probs <- dtms_rename(probs,c(fromvar,tovar,Pvar),c("from","to","P"))

  # Subset
  getthem <- probs$from%in%transient_states & probs$to%in%all_states
  probs <- subset(probs,subset=getthem)

  # Total number of transient and absorbing states
  s_states <- length(transient_states)
  a_states <- length(absorbing)
  n_states <- length(all_states)

  # Reshape
  Tm <- stats::reshape(probs[,c("from","to","P")],
                  timevar="to",
                  idvar="from",
                  direction="wide",
                  sep=reshapesep)

  # Edit a bit
  Tm[is.na(Tm)] <- 0
  keepnames <- Tm$from
  Tm <- Tm[,-1]

  # Generate matrix
  Tm <- as.matrix(Tm)
  rownames(Tm) <- keepnames

  # Column names
  oldnames <- strsplit(colnames(Tm),split=paste0("[",reshapesep,"]"))
  oldnames <- lapply(oldnames,function(x) x[2])
  colnames(Tm) <- unlist(oldnames)

  # Add "missing" starting states, if any
  addnames <- rownames(Tm)[!rownames(Tm)%in%colnames(Tm)]
  nadd <- length(addnames)
  if(nadd>0) {
    add <- matrix(data=0,ncol=nadd,nrow=dim(Tm)[1])
    colnames(add) <- addnames
    rownames(add) <- rownames(Tm)
    Tm <- cbind(Tm,add)
  }

  # Add potentially missing final states
  addnames <- colnames(Tm)[!colnames(Tm)%in%rownames(Tm)]
  nadd <- length(addnames)
  if(nadd>0) {
    add <- matrix(data=0,nrow=nadd,ncol=dim(Tm)[2])
    rownames(add) <- addnames
    colnames(add) <- colnames(Tm)
    Tm <- rbind(Tm,add)
  }

  # Add death (the column should already be there)
  Tm <- rbind(Tm,rep(0,n_states))
  rownames(Tm)[(s_states+1):n_states] <- absorbing

  # The dead stay dead (hopefully)
  if(a_states==1) Tm[absorbing,absorbing] <- 1
  if(a_states>1) diag(Tm[absorbing,absorbing]) <- 1

  # Sort a little
  Tm <- Tm[all_states,all_states]

  # Numbers please
  class(Tm) <- "numeric"

  # Make sure everyone dies at the end
  if(enforcedeath==T) {
    last_states <- paste(dtms$transient,max(dtms$timescale),sep=dtms$sep)
    if(length(absorbing)==1) {
      Tm[last_states,] <- 0
      Tm[last_states,absorbing] <- 1
    }
    if(length(absorbing)>1) {
      Tm[last_states,] <- 0
      Tm[last_states,absorbing[1]] <- 1
    }
  }

  # Rescale
  if(rescale) Tm <- t(apply(Tm,1,function(x) x/sum(x)))

  # Class
  class(Tm) <- c("dtms_matrix","matrix")

  # Return
  return(Tm)

}
#' @export
## Summary function for data in transition format
summary.dtms_data <- function(object,...) {
  dtms_data_summary(data=object,...)
}

#' @export
## Summary function for transition probabilities
summary.dtms_probs <- function(object,...) {
  dtms_probs_summary(probs=object,...)
}

#' @export
## Summary function for distributional results
summary.dtms_distr <- function(object,...) {
  dtms_distr_summary(distr=object,...)
}

#' @export
## Summary function for bootstrap
summary.dtms_boot <- function(object,...) {
  dtms_boot_summary(object,...)
}

#' @export
## Plotting function for transition probabilities
plot.dtms_probs <- function(x,...) {
  dtms_plot(probs=x,...)
}
#' Nonparametric estimates of transition probabilities
#'
#' @description
#' This function calculates nonparametric estimates of transition probabilities.
#' Standard errors assume that all observations are independent.
#'
#' @details
#' The argument `data` takes a data set in transition format. Predicted
#' transition probabilities are returned as a data frame, and not
#' as a transition matrix. While the latter is required for applying Markov
#' chain methods, the data frame is more convenient for viewing and
#' analyzing the transition probabilities themselves.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable with starting state in `data`. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state in `data`. Default is `to`.
#' @param timevar Character (optional), name of variable with time scale in `data`. Default is `time`.
#' @param Pvar Character (optional), name of variable with transition probabilities in the returned data frame. Default is `P`.
#' @param weights Character (optional). Name of variable with survey weights.
#' @param se Logical (optional), return standard errors of predicted probabilites. Default is `TRUE`.
#' @param vcov Logical (optional), return variance-covariance matrix of predicted probabilites. Default is `FALSE`.
#' @param CI Logical (optional), return confidence intervals? See details. Default is FALSE.
#' @param alpha Numeric (optional), if CI=TRUE, what confidence level is used? Default is 0.05.
#'
#' @returns A data frame with transition probabilities.
#' @export
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
#' ## Nonparametric transition probabilities
#' probs <- dtms_nonparametric(data=estdata,
#'                             dtms=simple)

dtms_nonparametric <- function(data,
                               dtms,
                               fromvar="from",
                               tovar="to",
                               timevar="time",
                               Pvar="P",
                               weights=NULL,
                               se=TRUE,
                               vcov=FALSE,
                               CI=FALSE,
                               alpha=0.05) {

  # Check
  dtms_proper(dtms)

  # Adjust time scale (transitions in the model)
  timescale_reduced <- dtms$timescale[-length(dtms$timescale)]
  ntime <- length(timescale_reduced)

  # Get full state space
  all_states <- paste(c(dtms$transient,dtms$absorbing))

  # Create empty frame
  model_frame <- expand.grid(from=dtms$transient,
                             to=all_states,
                             time=timescale_reduced,
                             stringsAsFactors=FALSE)

  # Get names right
  names(model_frame) <- c(fromvar,tovar,timevar)

  # Weights per transition
  if(is.null(weights)) data$COUNT <- 1 else
    data <- dtms_rename(data,weights,"COUNT")

  # Warning if missing values
  if(any(is.na(data[,c(fromvar,tovar,timevar)]))) warning("Missing values dropped")

  # Aggregate (denominators)
  formal1 <- paste0("COUNT~",fromvar,"+",timevar)
  formal1 <- stats::as.formula(formal1)
  denominators <- stats::aggregate(formal1,data,FUN=sum,drop=F)

  # Aggregate (numerators)
  formal2 <- paste0("COUNT~",fromvar,"+",tovar,"+",timevar)
  formal2 <- stats::as.formula(formal2)
  numerators <- stats::aggregate(formal2,data,FUN=sum,drop=F)

  # Merge
  probs <- merge(numerators,denominators,by=c(fromvar,timevar))
  model_frame <- merge(model_frame,probs,by=c(fromvar,tovar,timevar))

  # Replace missing with 0
  model_frame$COUNT.x[is.na(model_frame$COUNT.x)] <- 0

  # Calculate
  model_frame[,Pvar] <- model_frame$COUNT.x/model_frame$COUNT.y

  # Standard error/confidence interval/vcov?
  if(se|CI) {

  }

  # Warning if empty cells etc cause missing values
  if(any(is.na(model_frame[,Pvar]))) warning("Some probabilities are missing")

  # Values of starting state
  model_frame[,fromvar] <- paste(model_frame[,fromvar],model_frame[,timevar],sep=dtms$sep)

  # Values of receiving state (state name + time)
  rightrows <- model_frame[,tovar]%in%dtms$transient
  oldvalues <- model_frame[rightrows,tovar]
  timevalues <- model_frame[rightrows,timevar]+dtms$timestep
  model_frame[rightrows,tovar] <- paste(oldvalues,
                                        timevalues,
                                        sep=dtms$sep)

  # Only keep relevant variables
  model_frame <- model_frame[,names(model_frame)%in%c(fromvar,tovar,timevar,Pvar,"se","CIlow","CIup")]

  # Class
  class(model_frame) <- c("dtms_probs","data.frame")

  # Return
  return(model_frame)

}
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
#' Plotting transition probabilities
#'
#' @description
#' A simple function for plotting transition probabilities with base R. This is
#' fast, but it is much easier to produce nicer looking results with
#' dtms_simplify.
#'
#' @param probs Object with transition probabilities as created with \code{dtms_transitions}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable in `probs` with starting state. Default is `from`.
#' @param tovar Character (optional), name of variable in `probs` with receiving state. Default is `to`.
#' @param timevar Character (optional), name of variable in `probs` with time scale. Default is `time`.
#' @param Pvar Character (optional), name of variable in `probs` with transition probabilities. Default is `P`.
#' @param ... Further arguments passed to plot().
#'
#' @export
#'
#' @examples
#' ## Model setup
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
#' ## Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Plot
#' dtms_plot(probs=probs,
#'           dtms=simple)

dtms_plot <- function(probs,
                      dtms,
                      fromvar="from",
                      tovar="to",
                      timevar="time",
                      Pvar="P",
                      ...) {

  # Check
  if(!is.null(dtms)) dtms_proper(dtms)

  # Simplify data
  probs <- dtms_simplify(probs,
                         fromvar=fromvar,
                         tovar=tovar,
                         sep=dtms$sep)

  # All states
  if(!is.null(dtms)) {
    transient <- dtms$transient
    allstates <- c(dtms$transient,dtms$absorbing)
  } else {
    transient <- unique(probs[,fromvar])
    allstates <- unique(probs[,tovar])
  }

  # Number of rows and columns
  nrows <- length(transient)
  ncols <- length(allstates)

  # Number of panels for layout
  npanels <- ncols*nrows

  # Layout
  graphics::layout(matrix(1:npanels,nrow=nrows))

  # Max probs
  maxprob <- max(probs[,Pvar])
  maxprob <- round(maxprob,digits=1)

  # Plot
  for(fromstate in transient) {

    for(tostate in allstates) {

      tmp <- subset(probs,
                    probs[,fromvar]%in%fromstate &
                      probs[,tovar]%in%tostate)

      plot(x=tmp[,timevar],
           y=tmp[,Pvar],
           type="l",
           ylim=c(0,maxprob),
           main=paste(fromstate,"to",tostate),
           ...)

    }

  }

  # Back to normal layout
  graphics::layout(1)

}
#' Summarize transition probabilities
#'
#' @description
#' Provides several summary statistics on transition probabilities.
#'
#' @param probs Object with transition probabilities as created with \code{dtms_transitions}.
#' @param fromvar Character (optional), name of variable with starting state in `probs`. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state in `probs`. Default is `to`.
#' @param timevar Character (optional), name of variable with time scale in `probs`. Default is `time`.
#' @param Pvar Character (optional), name of variable with transition probabilities in `probs`. Default is `P`.
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
                               digits=4,
                               format="decimal",
                               sep="_") {

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
#' Markov chain with rewards
#'
#' @description
#' This function calculates the expected rewards by starting state in a
#' Markov chain with rewards.
#'
#' @param matrix Matrix with transition probabilities, as generated with \code{dtms_matrix}.
#' @param reward Matrix with rewards, has to be of same dimensions as `matrix`.
#' @param dtms dtms object, as created with \code{dtms}.
#'
#' @return A matrix with expected rewards.
#' @export
#'
#' @seealso [dtms_rewardmatrix()]
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
#' # Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Get transition matrix
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' ## Reward matrix
#' Rw <- diag(1,dim(Tp)[1])
#' ## State expectancies
#' dtms_reward(dtms=simple,
#'             matrix=Tp,
#'             reward=Rw)

dtms_reward <- function(matrix,
                        reward,
                        dtms) {

  # Check
  dtms_proper(dtms)

  # Starting states, long names
  starting <- dtms_combine(dtms$transient,
                           min(dtms$timescale),
                           sep=dtms$sep)

  # Number of starting and receiving states
  nstart <- length(starting)
  ntransient <- length(dtms$transient)
  nabsorbing <- length(dtms$absorbing)

  # Remove absorbing states
  nmatrix <- dtms_absorbing(matrix)

  # Fundamental matrix
  nstates <- dim(nmatrix)[1]
  Nmat <- solve(diag(1,nstates)-nmatrix)

  # Additional matrices
  Z <- cbind(diag(1,nstates),matrix(data=0,ncol=nabsorbing,nrow=nstates))
  one <- rep(1,nabsorbing+nstates)

  # e
  result <- Nmat %*% Z %*%(matrix %*% reward) %*% one
  colnames(result) <- "Reward"

  # Return result
  return(result)

}
#' Generate the reward matrix for a Markov chain with rewards
#'
#' @description
#' This function generates a reward matrix which can be used with
#' \code{dtms_reward}.
#'
#' @param dtms dtms object, as created with \code{dtms}.
#' @param starting Character (optional), name or names of starting states. If NULL (default) any transition to the state or states specififed with \code{receiving} will get the reward.
#' @param receiving Character, name or names of states to which transitioning generates the reward. Can be both transient or absorbing states.
#' @param reward Numeric, reward value to be placed in matrix.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#'
#' @return A matrix with rewards.
#' @export
#'
#' @seealso [dtms_reward()]
#'
#' @examples
#' ## Define model: Absorbing and transient states, time scale
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:20)
#' dtms_rewardmatrix(dtms=simple,receiving="B",reward=0.3)

dtms_rewardmatrix <- function(dtms,
                              starting=NULL,
                              receiving,
                              reward,
                              start_time=NULL,
                              end_time=NULL) {

  # Check
  dtms_proper(dtms)

  # Combine states and time
  transient_states <- dtms_combine(dtms$transient,dtms$timescale,sep=dtms$sep)
  absorbing <- paste(dtms$absorbing)
  all_states <- c(transient_states,absorbing)

  # Starting states, long names
  if(is.null(starting)) starting <- dtms$transient
  starting_time <- dtms$timescale
  ntime <- length(starting_time)
  if(!is.null(start_time)) starting_time <- starting_time[which(dtms$timescale==start_time):ntime]
  if(!is.null(end_time)) starting_time <- starting_time[1:which(dtms$timescale==end_time)]
  starting_states <- dtms_combine(starting,
                                  starting_time,
                                  sep=dtms$sep)

  # Receiving states, time
  receiving_time <- dtms$timescale
  ntime <- length(receiving_time)
  if(!is.null(start_time)) receiving_time <- receiving_time[which(dtms$timescale==start_time):ntime]
  if(!is.null(end_time)) receiving_time <- receiving_time[1:which(dtms$timescale==end_time)]

  # Receiving states, split between transient and absorbing
  rec_tra <- receiving[receiving%in%dtms$transient]
  rec_abs <- receiving[receiving%in%dtms$absorbing]

  if(length(rec_tra)>0) {
    rec_tra <- dtms_combine(receiving,
                            receiving_time,
                            sep=dtms$sep)
  }

  receiving_states <- c(rec_tra,rec_abs)

  # Number of states
  nall <- length(all_states)

  # Build empty matrix
  result <- matrix(data=0,ncol=nall,nrow=nall)
  rownames(result) <- all_states
  colnames(result) <- all_states

  # Fill
  result[starting_states,receiving_states] <- reward

  # Return result
  return(result)

}
#' Calculate the lifetime risk of ever reaching a state
#'
#' @description
#' The function `dtms_risk` calculates the (partial) lifetime risk of ever
#' reaching a state specified with the argument `risk`.
#'
#' @param matrix Matrix with transition probabilities, as generated with \code{dtms_matrix}.
#' @param risk Character, name of state(s) for which risk is of interest.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param start_distr Numeric (optional), distribution of starting states. If specified, average distribution over all starting states will be calculated.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#'
#' @return Probability of ever reaching state `risk`.
#' @export
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
#' # Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Get transition matrix
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' ## Get starting distribution
#' S <- dtms_start(dtms=simple,
#'                 data=estdata)
#' ## Lifetime risk
#' dtms_risk(dtms=simple,
#'           matrix=Tp,
#'           risk="A")

dtms_risk <- function(matrix,
                      risk,
                      dtms,
                      start_distr=NULL,
                      start_state=NULL,
                      start_time=NULL,
                      end_time=NULL) {

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=dtms$sep)

  # States of the transition matrix
  allstates <- rownames(matrix)

  # Partition of states: all states which belong to risk
  selectorU <- dtms_in(allstates,risk,dtms$sep)

  # Use end_time if specified
  if(!is.null(end_time)) {
    times <- dtms_gettime(allstates,dtms$sep)
    times <- times<=end_time
    times[is.na(times)] <- FALSE
    selectorU <- selectorU & times
  }

  # Invert selection
  selectorD <- !selectorU

  # New transition matrix
  newmatrix <- matrix[selectorD,selectorD]
  newstates <- rownames(newmatrix)
  nnewstates <- length(newstates)

  # Probability of moving to risk
  probrisk <- 1-rowSums(newmatrix)

  # Add probability
  newmatrix <- cbind(newmatrix,probrisk)
  newmatrix <- rbind(newmatrix,
                     c(rep(0,nnewstates),1))
  colnames(newmatrix)[nnewstates+1] <- "Risk"
  rownames(newmatrix)[nnewstates+1] <- "Risk"

  # Get N
  Umat <- dtms_absorbing(newmatrix)
  nstates <- dim(Umat)[1]
  Nmat <- solve(diag(1,nstates)-Umat)

  # Get R
  whichabsorbing <- which(diag(newmatrix)==1)
  R <- newmatrix[-whichabsorbing,whichabsorbing]

  # Get results
  results <- Nmat%*%R

  # Get results in shape
  result <- rep(1,length(starting))
  whererisk <- !start_state%in%risk
  result[whererisk] <- results[starting[whererisk],"Risk"]
  names(result) <- starting

  # Add average
  if(!is.null(start_distr)) {

    # Overall average
    tmp <- sum(start_distr*result)
    result <- c(result,tmp)
    names(result)[length(result)] <- "AVERAGE"

    # Average conditional on not starting in state
    tmp_distr <- start_distr[whererisk]/sum(start_distr[whererisk])
    tmp <- sum(tmp_distr*result[names(tmp_distr)])
    result <- c(result,tmp)
    names(result)[length(result)] <- "AVERAGE(COND.)"

  }

  # Return
  return(result)

}
#' Simplify state names
#'
#' @description
#' This function turns long state names into short state names. It is
#' particularly useful for plotting and when used in pipes, see the example.
#'
#' @param probs Object with transition probabilities as created with \code{dtms_transitions}.
#' @param fromvar Character (optional), name of variable in `probs` with starting state. Default is `from`.
#' @param tovar Character (optional), name of variable in `probs` with receiving state. Default is `to`.
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#'
#' @return Data frame
#' @export
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
#' ## Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Simplify
#' probs |>  dtms_simplify()
#' ## NOT RUN: requires ggplot2
#' # library(ggplot2)
#' # probs |>  dtms_simplify() |>
#' #   ggplot(aes(x=time,y=P,color=to)) +
#' #   geom_line() +
#' #   facet_wrap(~from)

dtms_simplify <- function(probs,
                          fromvar="from",
                          tovar="to",
                          sep="_") {

  # Simplify names
  probs[,fromvar] <- dtms_getstate(probs[,fromvar],sep=sep)
  probs[,tovar] <- dtms_getstate(probs[,tovar],sep=sep)

  # Return
  return(probs)

}
#' Simulation of Markov chains
#'
#' @description
#' This function simulates trajectories based on a Markov chain using the
#' `markovchain` package.
#'
#' @param matrix Matrix, a matrix of transition probabilities as created with \code{dtms_matrix()},
#' @param dtms dtms object, as created with \code{dtms}.
#' @param size Numeric, number of trajectories which will be simulated. Default is 100.
#' @param start_distr Numeric (optional), distribution of starting states. If NULL, starting states will be assumed to be equally distributed.
#' @param droplast Logical (optional), drop final time step after the time scale in which every unit is absorbed? Default is TRUE.
#' @param varnames Character (optional), suffix for variable names in simulated data. Will be pasted with values of the timescale. Default is "T_".
#'
#' @return A data frame with simulated trajectories in wide format.
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:19)
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' fit <- dtms_fit(data=estdata,package="mclogit")
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' S <- dtms_start(dtms=simple,
#'                 data=estdata)
#' dtms_simulate(matrix=Tp,dtms=simple)

dtms_simulate <- function(matrix,
                          dtms,
                          size=100,
                          start_distr=NULL,
                          droplast=TRUE,
                          varnames="T_") {

  # Load markovchain package, because of "new" below
  requireNamespace("markovchain")

  # Check
  dtms_proper(dtms)

  # Time
  start_time <- min(dtms$timescale)
  ntime <- length(dtms$timescale)

  # Starting states, short and long names
  start_state <- dtms$transient
  starting <- dtms_combine(start_state,start_time,sep=dtms$sep)
  nstarting <- length(starting)

  # Starting distribution
  if(is.null(start_distr)) start_distr <- rep(1/nstarting,nstarting)

  # All states
  all_states <- colnames(matrix)

  # Set class of matrix
  class(matrix) <- "matrix"

  # Setup
  sim <- methods::new("markovchain",
                      states = all_states,
                      transitionMatrix = as.matrix(matrix),
                      name = "sim")

  # Data frame
  simdata <- data.frame(matrix(nrow=0,ncol=ntime))

  # Simulate
  for(i in 1:size) {

    # Draw initial state
    initial_state <- sample(starting,
                            size=1,
                            prob=start_distr)

    # Generate rest of sequence
    simseq <- markovchain::rmarkovchain(n = ntime-as.numeric(droplast),
                                        object = sim,
                                        t0 = initial_state,
                                        include.t0=T)

    # Put in data frame
    simdata <- rbind(simdata,simseq)

  }

  # Nicer names
  if(droplast) names(simdata) <- paste0(varnames,dtms$timescale) else
    names(simdata) <- paste0(varnames,c(dtms$timescale,max(dtms$timescale)+dtms$timestep))

  # Return
  return(simdata)

}
#' Tabulate starting distribution
#'
#' @description
#' Tabulates the starting distribution.
#'
#' @details
#' Per default, the starting distribution is the distribution of transient
#' states at the first value of the time scale in the data. This can be
#' changed to any value of the time scale, and any set of states. The
#' distribution can also be conditional on further covariate values which can be
#' specified with the argument `variables`.
#'
#' `variables` takes a named list where each entry of the list is named like
#' the corresponding variable and with the values to be selected.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param variables List (optional), a named list with covariate values which are used to restrict the data.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If several values are specified, the average distribution over all these values is calculated. In this case the first value specified with this argument is used to construct the long state name. If NULL (default) first value of time scale will be used.
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is `from`.
#' @param timevar Character (optional), name of variable in `data` with time scale. Default is `time`.
#' @param weights Character (optional). Name of variable with survey weights.
#'
#' @return Returns a table of the starting distribution.
#' @export
#'
#' @examples
## Define model: Absorbing and transient states, time scale
#' work <- dtms(transient=c("Working","Non-working","Retired"),
#'              absorbing="Dead",
#'              timescale=50:99)
#' ## Reshape
#' estdata <- dtms_format(data=workdata,
#'                        dtms=work,
#'                        idvar="ID",
#'                        timevar="Age",
#'                        statevar="State")
#' ## Drop dead-to-dead transitions etc
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=work)
#' ## Starting distributions
#' # Men
#' Sm <- dtms_start(dtms=work,
#'                  data=estdata,
#'                  variables=list(Gender=0))
#' # Women
#' Sw <- dtms_start(dtms=work,
#'                  data=estdata,
#'                  variables=list(Gender=1))

dtms_start <- function(data,
                       dtms,
                       variables=NULL,
                       start_state=NULL,
                       start_time=NULL,
                       fromvar="from",
                       timevar="time",
                       weights=NULL) {

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time[1],sep=dtms$sep)

  # Restrict data: States and time
  data <- data[data[,fromvar]%in%start_state,]
  data <- data[data[,timevar]%in%start_time,]

  # Apply restrictions, if any
  varnames <- names(variables)
  for(var in varnames) {
    data <- data[data[,var]%in%variables[[var]],]
  }

  # Tabulate
  if(is.null(weights)) {
    tab <- data[,fromvar] |> table() |> prop.table()
  } else {
    # https://stackoverflow.com/questions/18585977/frequency-tables-with-weighted-data-in-r
    tmp <- stats::aggregate(x = data[,weights], by = list(data[,fromvar]), FUN = sum)
    tab <- tmp[,2]
    names(tab) <- tmp[,1]
    tab <- prop.table(tab)
  }

  # Match with starting names
  result <- numeric(length(start_state))
  names(result) <- start_state
  matchnames <- names(tab)[names(tab)%in%start_state]
  result[matchnames] <- tab[matchnames]

  # Fix if issues and procude warning
  wrong <- which(result%in%c(NA,NaN,Inf,-Inf))
  if(length(wrong)>0) {
    warning("Something might have gone wrong with the starting distribution")
    result[wrong] <- 0
    result <- result/sum(result)
  }

  # Proper names
  names(result) <- starting

  # Return
  return(result)

}
#' Calculate the survivorship function
#'
#' @description
#' Calculates the proportion of units surviving up to certain values of the
#' time scale.
#'
#' @details
#' If a distribution of the starting states is provided with `start_distr` the
#' output table has an additional row, showing the waiting time distribution
#' unconditional on the starting state.
#'
#' @param matrix Matrix with transition probabilities, as generated with \code{dtms_matrix}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param start_distr Numeric (optional), distribution of starting states. If specified, average distribution over all starting states will be calculated.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#'
#' @return A table with the survivorship function.
#' @export
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
#' ## Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Get transition matrix
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' ## Get starting distribution
#' S <- dtms_start(dtms=simple,
#'                 data=estdata)
#' ## Distribution of visits
#' dtms_survivor(dtms=simple,
#'               matrix=Tp,
#'               start_distr=S)


dtms_survivor <- function(matrix,
                          dtms,
                          start_distr=NULL,
                          start_state=NULL,
                          start_time=NULL,
                          end_time=NULL) {

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=dtms$sep)

  # Matrix for results
  Umat <- dtms_absorbing(matrix=matrix)
  nstates <- dim(Umat)[1]
  distr <- matrix(data=NA,ncol=1,nrow=nstates)
  rownames(distr) <- rownames(Umat)

  # First step: Pr(Time to absorbtion=1)
  distr[,1] <- rowSums(diag(1,nstates)-Umat)
  step <- 1
  if(is.null(end_time)) max_steps <- length(dtms$timescale) else {
    max_steps <- length(1:which(dtms$timescale==end_time))
  }

  # Iterate
  while(!all(round(colSums(distr),digits=2)==1) & step<max_steps) {
    step <- step+1
    distr <- cbind(distr,rowSums(dtms_mtexp(Umat,step-1)-dtms_mtexp(Umat,step)))
  }

  # Edit to survivor function
  distr <- apply(distr,1,function(x) 1-cumsum(x))
  distr <- t(distr)
  distr <- cbind(1,distr)
  colnames(distr) <- 0:(dim(distr)[2]-1)

  # Select starting states
  result <- distr[starting,]

  # Shift states starting after t=0
  statenames <- rownames(result)
  for(state in statenames) {
    shorttime <- dtms_gettime(state,sep=dtms$sep)
    wherestart <- which(dtms$timescale==shorttime)
    whatlength <- length(wherestart:max_steps)
    movevalues <- result[state,1:whatlength]
    result[state,wherestart:max_steps] <- movevalues
    result[state,1:(wherestart-1)] <- 1
  }

  # Average
  if(!is.null(start_distr)) {

    # Overall average
    AVERAGE <- start_distr%*%result
    rownames(AVERAGE) <- "AVERAGE"
    result <- rbind(result,AVERAGE)

  }

  # Assign class
  class(result) <- c("dtms_distr","matrix")

  # Results
  return(result)

} # End of function
#' Predict transition probabilities
#'
#' @description
#' `dtms_transitions()` predicts transition probabilities based on a model
#' estimated with `dtms_fit()`.
#'
#' @details
#' Predicted transition probabilities are returned as a data frame, and not
#' as a transition matrix. While the latter is required for applying Markov
#' chain methods, the data frame is more convenient for viewing and
#' analyzing the transition probabilities themselves.
#'
#' Depending on the model specification, the prediction of transition
#' probabilities will require values for predictor variables which can be
#' specified with the argument `controls`. This is done using a named list
#' where each entry name must correspond to a variable name in the model.
#' For time-constant variables, each list entry is of length one and provides
#' a value for the corresponding time-constant variable. For time-varying
#' variables, each entry must have the length of the time scale minus one, and
#' provide a value for each (potential) transition in the model; i.e., starting
#' from time t=0, starting from time t=1, etc., until time t=T-1. Alternatively,
#' it can be of the same length as the time scale; in this case, the last value
#' is dismissed.
#'
#' If `vcov=TRUE` the full variance-covariance matrix of the transition
#' probabilities will be returned instead of the transition probabilities. If
#' `CI=TRUE`, confidence intervals will be returned. Note that the calculation
#' uses a normal approximation and results below 0 or above 1 are possible.
#'
#' The argument `dropvar` controls whether the covariate values used for
#' prediction are dropped. If `FALSE` each row of the resulting data frame will
#' have the covariate values which were used to predict the corresponding
#' probability.
#'
#' @param model Model estimated with \code{dtms_fit}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param controls List (optional) with values for predictors (see details).
#' @param se Logical (optional), return standard errors of predicted probabilites. Default is `TRUE`.
#' @param vcov Logical (optional), return variance-covariance matrix of predicted probabilites. Default is `FALSE`.
#' @param CI Logical (optional), return confidence intervals? See details. Default is FALSE.
#' @param alpha Numeric (optional), if CI=TRUE, what confidence level is used? Default is 0.05.
#' @param dropvar Logical (optional), should covariate values used for prediction be returned (see details). Default is `TRUE`.
#' @param fromvar Character (optional), name of variable with starting state in the returned data frame. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state in the returned data frame. Default is `to`.
#' @param timevar Character (optional), name of variable with time scale in the returned data frame. Default is `time`.
#' @param Pvar Character (optional), name of variable with transition probabilities in the returned data frame. Default is `P`.
#'
#' @return A data frame with transition probabilities.
#' @export
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
#' ## Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)

dtms_transitions <- function(model,
                             dtms,
                             controls=NULL,
                             dropvar=TRUE,
                             timevar="time",
                             fromvar="from",
                             tovar="to",
                             Pvar="P",
                             se=TRUE,
                             vcov=FALSE,
                             CI=FALSE,
                             alpha=0.05) {

  # Check
  dtms_proper(dtms)

  # Adjust time scale (transitions in the model)
  timescale_reduced <- dtms$timescale[-length(dtms$timescale)]
  ntime <- length(timescale_reduced)

  # Get full state space
  all_states <- paste(c(dtms$transient,dtms$absorbing))

  # Create empty frame
  model_frame <- expand.grid(from=dtms$transient,
                             time=timescale_reduced,
                             stringsAsFactors=FALSE)

  # Get names right
  names(model_frame) <- c(fromvar,timevar)

  # Deal with controls
  varnames <- names(controls)
  for(var in varnames) {

    # Get values
    value <- controls[[var]]
    nvalue <- length(value)

    # Check
    if(!nvalue%in%c(1,ntime,ntime+1)) stop("Wrong number of time-varying values")

    # Act depending on how many values
    if(nvalue==1) model_frame[var] <- value

    if(nvalue==ntime) {
      assign_values <- match(model_frame[,timevar],timescale_reduced)
      model_frame[var] <- value[assign_values]
    }

    if(nvalue==ntime+1) {
      value <- value[-nvalue]
      assign_values <- match(model_frame[,timevar],timescale_reduced)
      model_frame[var] <- value[assign_values]
    }


  }

  # Predict
  if(inherits(model,c("vgam","mclogit"))) {
    model_frame[,all_states] <- stats::predict(model,model_frame,"response")[,all_states]
  }

  if(inherits(model,"nnet")) {
    model_frame[,all_states] <- stats::predict(model,model_frame,"probs")[,all_states]
  }

  # Values of starting state
  model_frame[,fromvar] <- paste(model_frame[,fromvar],model_frame[,timevar],sep=dtms$sep)

  # Reshape
  model_frame <- stats::reshape(model_frame,
                                varying=all_states,
                                idvar=fromvar,
                                timevar=tovar,
                                times=all_states,
                                direction="long",
                                v.names=Pvar)

  # SE/CI/vcov
  if(se|vcov|CI) {

    # Simplify starting state (needed for model.matrix below)
    model_frame[,fromvar] <- dtms_simplify(model_frame)$from

    # Coefficients
    if(inherits(model,"mclogit")) {
      C <- stats::coef(model)
      Cstates <- dtms_getstate(names(C),sep="~")
      Cstates <- unique(Cstates)
      C <- matrix(data=C,
                  ncol=length(all_states)-1,
                  byrow=T)
    }

    if(inherits(model,"vgam")) {
      C <- stats::coef(model)
      C <- matrix(C,
                  ncol=length(all_states)-1,
                  byrow=T)
      Cstates <- model@extra$colnames.y[-model@extra$use.refLevel]
    }

    if(inherits(model,"nnet")) {
      C <- stats::coef(model)
      Cstates <- rownames(C)
      C <- t(C)
    }

    # vcov of coefficients
    Vml <- stats::vcov(model)

    if(inherits(model,"mclogit")) {
      ordernames <- sort(rownames(Vml))
      Vml <- Vml[ordernames,ordernames]
      Vnames <- dtms_getstate(rownames(Vml),sep=c("~"))
    }

    if(inherits(model,"vgam")) {

      # Get nice names (assigned below)
      nicenames <- dtms_getstate(rownames(Vml),sep=c(":"))
      nicenames <- unique(nicenames)
      nicenames <- sort(dtms_combine(Cstates,nicenames,sep=":"))

      # Reorder
      newnames <- unlist(lapply(strsplit(colnames(Vml),split=":"),function(x) paste0(x[2],":",x[1])))
      colnames(Vml) <- rownames(Vml) <- newnames
      ordernames <- sort(rownames(Vml))
      Vml <- Vml[ordernames,ordernames]

      # Assign nice names
      colnames(Vml) <- rownames(Vml) <- nicenames

      # State names
      Vnames <- dtms_getstate(rownames(Vml),sep=c(":"))
    }

    if(inherits(model,"nnet")) Vnames <- dtms_getstate(rownames(Vml),sep=c(":"))

    # Number of probabilities, coefs, states
    nprobs <- dim(model_frame)[1]
    ncoef <- dim(Vml)[1]
    nstates <- length(all_states)

    # Model matrix
    form <- stats::formula(model)
    mm <- stats::model.matrix(object=form,data=model_frame)

    # Scores (denominator for predicted prob)
    dscores <- matrix(data=1,
                      ncol=nstates,
                      nrow=nprobs)
    colnames(dscores) <- sort(all_states)
    dscores[,Cstates] <- exp(mm%*%C)

    # Full score (numerator for predicted prob)
    fullscore <- rowSums(dscores)

    # Parts of full derivative (n'*z-n*z')/z^2
    Z <- matrix(data=fullscore,
                nrow=nprobs,
                ncol=ncoef)

    N <- stats::model.matrix(object=~to,data=model_frame)
    N[,1] <- N[,1]-rowSums(N[,-1])
    N <- rowSums(N*dscores)
    N <- matrix(data=N,
                nrow=nprobs,
                ncol=ncoef)

    varvalues <- do.call("cbind",replicate(nstates-1,mm,simplify=FALSE))
    Zdash <- varvalues*dscores[,Vnames]

    Ndash <- outer(model_frame$to,Vnames,FUN=`==`)
    Ndash <- Ndash*N*varvalues

    # Matrix of derivatives
    G <- (Ndash*Z-N*Zdash)/(Z^2)

    # Re-order
    if(inherits(model,"mclogit")) colnames(G) <- dtms_combine(Cstates,colnames(mm),sep="~")
    if(inherits(model,"vgam")) colnames(G) <- dtms_combine(Cstates,colnames(mm),sep=":")
    if(inherits(model,"nnet")) colnames(G) <- dtms_combine(Cstates,colnames(mm),sep=":")

    G <- G[,rownames(Vml)]

    # Vcov matrix
    Vp <- G%*%Vml%*%t(G)

    # Return vcov matrix
    if(vcov) return(Vp)

    # Full starting state
    model_frame[,fromvar] <- paste(model_frame[,fromvar],model_frame[,timevar],sep=dtms$sep)

    # SE
    if(se) model_frame$se <- sqrt(diag(Vp))

    # CI
    if (CI) {
      z <- (1-alpha/2)
      z <- stats::qnorm(z)
      se <- sqrt(diag(Vp))
      model_frame$CIlow <- model_frame[,Pvar]-z*se
      model_frame$CIup <- model_frame[,Pvar]+z*se
    }

  }

  # Values of receiving state (state name + time)
  rightrows <- model_frame[,tovar]%in%dtms$transient
  oldvalues <- model_frame[rightrows,tovar]
  timevalues <- model_frame[rightrows,timevar]+dtms$timestep
  model_frame[rightrows,tovar] <- paste(oldvalues,
                                        timevalues,
                                        sep=dtms$sep)

  # Drop row names
  rownames(model_frame) <- NULL

  # Drop covariate values for prediction
  if(dropvar) {
    model_frame <- model_frame[,names(model_frame)%in%c(fromvar,tovar,timevar,Pvar,"se","CIlow","CIup")]
  }

  # Class
  class(model_frame) <- c("dtms_probs","data.frame")

  # Return
  return(model_frame)

}
#' Calculate the distribution of the time spent in a subset of states
#'
#' @description
#' Calculates the distribution of the time spent in a state or a subset of states.
#'
#' @details
#' The state(s) which count to the time are specified with the argument `risk`.
#' If several states are specified, the resulting distribution refers to the
#' lifetime spent in any of the specified states.
#'
#' In a discrete-time model, the time spent in a state depends on assumptions
#' about when transitions happen. Currently, this functions supports two
#' variants which can be specified with the argument `method`: mid-interval
#' transitions can be selected with the option `mid` and imply that transitions
#' happen at the middle of the time interval; and the option `end` assumes
#' that instead transitions happen at the end of the interval. In this latter
#' case the distribution of the time spent in a state is equivalent to the
#' number of visits to that state. The calculation takes the step length of
#' the time scale into account as specified by the `dtms` object. If the
#' step length is not one fixed value, the first entry of `dtms$timestep` will
#' be used.
#'
#' If a distribution of the starting states is provided with `start_distr` the
#' output table has two additional rows. One shows the distribution
#' unconditional on the starting state. The other shows the distribution
#' conditional on not starting in any state of the risk set.
#'
#' @param matrix Matrix with transition probabilities, as generated with \code{dtms_matrix}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param risk Character, name of state(s) for which risk is of interest.
#' @param start_distr Numeric (optional), distribution of starting states. If specified, average distribution over all starting states will be calculated.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#' @param method Character (optional), do transitions happen mid-interval (`mid`, default) or at the end of the interval (`end`), see details.
#' @param total Logical, should total of distribution be shown (always sums to 1)? Default is FALSE.
#'
#' @return A table with the distribution of time spent in a subset of states.
#' @export
#'
#' @seealso
#' \code{\link{dtms_distr_summary}} to help with summarizing the resulting distribution.
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
#' ## Fit model
#' fit <- dtms_fit(data=estdata)
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Get transition matrix
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)
#' ## Get starting distribution
#' S <- dtms_start(dtms=simple,
#'                 data=estdata)
#' ## Distribution of visits
#' dtms_visits(dtms=simple,
#'             matrix=Tp,
#'             risk="A",
#'             start_distr=S,
#'             total=TRUE)

dtms_visits <- function(matrix,
                        dtms,
                        risk,
                        start_time=NULL,
                        start_state=NULL,
                        start_distr=NULL,
                        end_time=NULL,
                        method="mid",
                        total=F) {

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=dtms$sep)

  # States of the transition matrix
  allstates <- rownames(matrix)
  nstates <- length(allstates)

  # Partition of states: all states which belong to risk
  selectorU <- dtms_in(allstates,risk,dtms$sep)

  # Use end_time if specified
  if(!is.null(end_time)) {
    times <- dtms_gettime(allstates,dtms$sep)
    times <- times<=end_time
    times[is.na(times)] <- FALSE
    selectorU <- selectorU & times
  }

  # Invert
  selectorD <- !selectorU

  # Empty partitions of transition matrix
  U_U <- matrix(data=0,nrow=nstates,ncol=nstates)
  U_D <- matrix(data=0,nrow=nstates,ncol=nstates)
  U_UD <- matrix(data=0,nrow=nstates,ncol=nstates)

  # Partition for mid-interval transitions
  if(method=="mid") {
    U_U[selectorU,selectorU] <- matrix[selectorU,selectorU]
    U_D[selectorD,selectorD] <- matrix[selectorD,selectorD]
    U_UD[selectorD,selectorU] <- matrix[selectorD,selectorU]
    U_UD[selectorU,selectorD] <- matrix[selectorU,selectorD]
  }

  # Partition for transitions at end of interval
  if(method=="end") {
    U_U[,selectorU] <- matrix[,selectorU]
    U_D[,selectorD] <- matrix[,selectorD]
  }

  # Get n
  if(is.null(end_time)) n <- length(dtms$timescale) else
    n <- which(end_time==dtms$timescale)

  # Time steps
  t_series <- 0:n
  if(method=="mid") d_series <- seq(0,n,by=0.5) else d_series <- t_series
  t_transitions <- length(t_series)
  t_steps <- 1:t_transitions

  # Generate initial conditions
  initial_conditions <- vector("list",t_transitions)
  for(i in t_steps) {
    initial_conditions[[i]] <- matrix(data=0,nrow=nstates,ncol=nstates)
    initial_conditions[[i]][selectorD,selectorD] <- dtms_mtexp(matrix[selectorD,selectorD],(i-1))
    rownames(initial_conditions[[i]]) <- allstates
    colnames(initial_conditions[[i]]) <- allstates
  }
  names(initial_conditions) <- paste("F",t_series,"0",sep="_")

  # Generate special case
  initial_conditions[["F_0_0"]] <- diag(1,nstates)
  rownames(initial_conditions[["F_0_0"]]) <- allstates
  colnames(initial_conditions[["F_0_0"]]) <- allstates

  # Generate matrices for which k>=n+1
  kgeq <- vector("list",t_transitions)
  for(i in t_steps) {
    kgeq[[i]] <- dtms_mtexp(matrix,i-1)
    rownames(kgeq[[i]]) <- allstates
    colnames(kgeq[[i]]) <- allstates
  }
  names(kgeq) <- paste("F",t_series,t_steps,sep="_")
  if(method=="mid") {
    kgeq_tmp <- kgeq
    names(kgeq_tmp) <- paste("F",t_series,t_steps-0.5,sep="_")
    kgeq <- c(kgeq,kgeq_tmp)
  }

  # Set 'negative' matrices to zero
  neg_matrices <- vector("list",t_transitions)
  for(i in t_steps) {
    neg_matrices[[i]] <- matrix(data=0,nrow=nstates,ncol=nstates)
    rownames(neg_matrices[[i]]) <- allstates
    colnames(neg_matrices[[i]]) <- allstates
  }
  names(neg_matrices) <- paste("F",t_series,"-0.5",sep="_")

  # Combine initial conditions and k>=n+1
  results <- c(initial_conditions,kgeq,neg_matrices)

  # Generate results step by step
  for(i in 1:n) {
    which_d <-  d_series[-1][d_series[-1]<=i]
    n_d <- length(which_d)
    tmp <- vector("list",n_d)
    names(tmp) <- paste("F",i,which_d,sep="_")

    for(k in 1:n_d) { # Equation (7)
      if(method=="mid") tmp[[paste("F",i,which_d[k],sep="_")]] <-
          U_U%*%results[[paste("F",i-1,which_d[k]-1,sep="_")]] +
          U_D%*%results[[paste("F",i-1,which_d[k],sep="_")]]   +
          U_UD%*%results[[paste("F",i-1,which_d[k]-0.5,sep="_")]]
      if(method=="end") tmp[[paste("F",i,which_d[k],sep="_")]] <-
          U_U%*%results[[paste("F",i-1,which_d[k]-1,sep="_")]] +
          U_D%*%results[[paste("F",i-1,which_d[k],sep="_")]]
      colnames(tmp[[paste("F",i,which_d[k],sep="_")]]) <- allstates
      rownames(tmp[[paste("F",i,which_d[k],sep="_")]]) <- allstates
    }

    results <- c(results,tmp)
  }

  # Get results of interest
  if(method=="mid") results <- results[paste("F",n,c(d_series,n+0.5,n+1),sep="_")]
  if(method=="end") results <- results[paste("F",n,c(d_series,n+1),sep="_")]

  # Distribution of k<=x
  tmp <- unlist(lapply(results, function(y) rowSums(y)[starting]))

  # Result object
  result <- matrix(data=tmp,nrow=length(starting))
  rownames(result) <- paste0("start:",starting)

  # Distribution of k=x
  result <- cbind(result[,1],t(apply(result,1,diff)))

  # Column names
  if(method=="mid") colnames(result) <- paste(c(d_series,n+0.5,n+1)*dtms$timestep[1])
  if(method=="end") colnames(result) <- paste(c(d_series,n+1)*dtms$timestep[1])

  # Drop last col if mid-transitions
  if(method=="mid") result <- result[,which(colnames(result)!=paste0(n+1))]

  # Average
  if(!is.null(start_distr)) {

    # Overall average
    AVERAGE <- start_distr%*%result
    rownames(AVERAGE) <- "AVERAGE"
    result <- rbind(result,AVERAGE)

    # Conditional on not starting in state in risk set
    whererisk <- !dtms$transient%in%risk
    tmp_distr <- start_distr[whererisk]/sum(start_distr[whererisk])
    tmp_names <- paste0("start:",names(tmp_distr))
    tmp <- tmp_distr%*%result[tmp_names,]
    result <- rbind(result,tmp)
    rownames(result)[dim(result)[1]] <- "AVERAGE(COND.)"

  }

  # Total
  if(total) {
    TOTAL <- rowSums(result)
    result <- cbind(result,TOTAL)
  }

  # Assign class
  class(result) <- c("dtms_distr","matrix")

  # Return
  return(result)

}
