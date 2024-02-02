#' Calculate the distribution of the time spent in a subset of states
#'
#' @description
#' A short description...
#'
#' @details
#' Mid-interval transitions and transitions at end of interval. Single state
#' versus subset of states. Partial distribution.
#'
#' @param matrix Matrix with transition probabilities, as generated with `dtms_matrix`.
#' @param dtms DTMS object as created with `dtms`.
#' @param risk Character (required), names of one or several states for which the time spent should be calculated.
#' @param start_state Character (optional), names of one or several starting states. If NULL (default), all transient states will be considered.
#' @param start_time Numeric (optional), value of time scale at start. If NULL (default), first value of time scale is used.
#' @param end_time Numeric (optional), value of time scale at end. If NULL (default), last value of time scale is used.
#' @param start_distr Numeric (optional), distribution of starting states. Needs to be consistent with starting states.
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

dtms_visits <- function(matrix,# Matrix with transition probabilities generated with dtms_matrix
                        transient=NULL, # Names of transient states
                        timescale=NULL, # Time scale
                        timestep=NULL, # Step length
                        dtms=NULL,# DTMS model
                        risk, # name of state(s) for which risk is of interest
                        start_time=NULL,# Starting time, will be lowest time in 'dtms' if NULL
                        start_state=NULL,# Names of starting states, will be all transient states in 'dtms' if NULL
                        start_distr=NULL,# Distribution of starting states for average
                        end_time=NULL, # Time up to which lifetime risks are calculated, making them partial lifetime risks
                        method="mid", # Mid-interval transitions or end of interval transitions, "mid" or "end"
                        sep="_",# Separator
                        total=F) { # Should total risk be added?

  # Use dtms if provided
  if(!is.null(dtms) & class(dtms)[2]=="dtms") {
    if(is.null(transient)) transient <- dtms$transient
    if(is.null(timescale)) {
      timescale <- dtms$timescale
      timescale <- timescale[-length(timescale)]
    }
    if(is.null(timestep)) timestep <- dtms$timestep
  }

  # Starting states
  if(is.null(start_state)) {
    starting <- levels(interaction(transient,min(timescale),sep=sep))
  } else {
    starting <- levels(interaction(start_state,start_time,sep=sep))
  }

  # States of the transition matrix
  states <- rownames(matrix)
  nstates <- length(states)

  # Partition of states: all states which belong to risk
  selectorU <- unlist(lapply(strsplit(states,sep),function(y) any(risk%in%y) ))

  # Use end time
  if(!is.null(end_time)) {
    times <- as.numeric(unlist(lapply(strsplit(states,split=sep),function(y) y[2] )))
    times <- times<=end_time
    times[!is.logical(times)] <- F
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
    rownames(initial_conditions[[i]]) <- states
    colnames(initial_conditions[[i]]) <- states
  }
  names(initial_conditions) <- paste("F",t_series,"0",sep="_")

  # Generate special case
  initial_conditions[["F_0_0"]] <- diag(1,nstates)
  rownames(initial_conditions[["F_0_0"]]) <- states
  colnames(initial_conditions[["F_0_0"]]) <- states

  # Generate matrices for which k>=n+1
  kgeq <- vector("list",t_transitions)
  for(i in t_steps) {
    kgeq[[i]] <- Biodem::mtx.exp(matrix,i-1)
    rownames(kgeq[[i]]) <- states
    colnames(kgeq[[i]]) <- states
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
    rownames(neg_matrices[[i]]) <- states
    colnames(neg_matrices[[i]]) <- states
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
      colnames(tmp[[paste("F",i,which_d[k],sep="_")]]) <- states
      rownames(tmp[[paste("F",i,which_d[k],sep="_")]]) <- states
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

  # Return
  return(result)

}  # End of function
