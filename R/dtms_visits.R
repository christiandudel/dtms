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
