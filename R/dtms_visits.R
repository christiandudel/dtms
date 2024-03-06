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
#' number of visits to that state.
#'
#' If a distribution of the starting states is provided with `start_distr` the
#' output table has two additional rows. One shows the distribution
#' unconditional on the starting state. The other shows the distribution
#' conditional on not starting in any state of the risk set.
#'
#' @param matrix Matrix with transition probabilities, as generated with `dtms_matrix`.
#' @param dtms DTMS object as created with `dtms`.
#' @param risk Character (required), names of one or several states for which the time spent should be calculated.
#' @param start_state Character (optional), names of one or several starting states. If NULL (default), all transient states will be considered.
#' @param start_time Numeric (optional), value of time scale at start. If NULL (default), first value of time scale is used.
#' @param end_time Numeric (optional), value of time scale at end. If NULL (default), last value of time scale is used.
#' @param start_distr Numeric (optional), distribution of starting states. Needs to be consistent with starting states. If provided, average distribution is provided; see details.
#' @param method Character, do transitions happen mid-interval (`mid`, default) or at the end of the interval (`end`), see details.
#' @param total Logical, should total of distribution be shown? Default is FALSE, as the total always is 1.
#' @param transient Character (optional), short names of transient states. If NULL (default) transient states are taken from `dtms` object.
#' @param timescale Numeric (optional), values of time scale. If NULL (default) obtained from `dtms` object.
#' @param timestep Numeric (optional), step length of time scale. If NULL (default) obtained from `dtms` object.
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#'
#' @return A table with the distribution of time spent in a subset of states.
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
#' dtms_visits(dtms=simple,
#'             matrix=Tp,
#'             risk="A",
#'             start_distr=S,
#'             total=TRUE)

dtms_visits <- function(matrix,
                        dtms=NULL,
                        risk,
                        start_time=NULL,
                        start_state=NULL,
                        start_distr=NULL,
                        end_time=NULL,
                        method="mid",
                        transient=NULL,
                        timescale=NULL,
                        timestep=NULL,
                        sep="_",
                        total=F) {

  # Use dtms if provided
  if(!is.null(dtms)) {

    # Check
    dtms_proper(dtms)

    # Use values
    transient <- dtms$transient
    timescale <- dtms$timescale
    timestep <- dtms$timestep
  }

  # Starting state and time
  if(is.null(start_state)) start_state <- transient
  if(is.null(start_time)) start_time <- min(timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=sep)

  # States of the transition matrix
  allstates <- rownames(matrix)
  nstates <- length(allstates)

  # Partition of states: all states which belong to risk
  selectorU <- dtms_in(allstates,risk,sep)

  # Use end_time if specified
  if(!is.null(end_time)) {
    times <- dtms_gettime(allstates,sep)
    times <- times<=end_time
    times[!is.logical(times)] <- F
    selector <- selector & times
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
  if(is.null(end_time)) n <- length(timescale) else n <- which(end_time==timescale)

  # Time steps
  t_series <- 0:n
  if(method=="mid") d_series <- seq(0,n,by=0.5) else d_series <- t_series
  t_transitions <- length(t_series)
  t_steps <- 1:t_transitions

  # Generate initial conditions
  initial_conditions <- vector("list",t_transitions)
  for(i in t_steps) {
    initial_conditions[[i]] <- matrix(data=0,nrow=nstates,ncol=nstates)
    initial_conditions[[i]][selectorD,selectorD] <- Biodem::mtx.exp(matrix[selectorD,selectorD],(i-1))
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
    kgeq[[i]] <- Biodem::mtx.exp(matrix,i-1)
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
      if(method=="mid") tmp[[paste("F",i,which_d[k],sep="_")]] <- U_U%*%results[[paste("F",i-1,which_d[k]-1,sep="_")]] +
          U_D%*%results[[paste("F",i-1,which_d[k],sep="_")]]   +
          U_UD%*%results[[paste("F",i-1,which_d[k]-0.5,sep="_")]]
      if(method=="end") tmp[[paste("F",i,which_d[k],sep="_")]] <- U_U%*%results[[paste("F",i-1,which_d[k]-1,sep="_")]] +
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
  if(method=="mid") colnames(result) <- paste(c(d_series,n+0.5,n+1)*timestep)
  if(method=="end") colnames(result) <- paste(c(d_series,n+1)*timestep)

  # Drop last col if mid-transitions
  if(method=="mid") result <- result[,which(colnames(result)!=paste0(n+1))]

  # Average
  if(!is.null(start_distr)) {

    # Overall average
    AVERAGE <- start_distr%*%result
    rownames(AVERAGE) <- "AVERAGE"
    result <- rbind(result,AVERAGE)

    # Conditional on not starting in state in risk set
    whererisk <- !transient%in%risk
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
