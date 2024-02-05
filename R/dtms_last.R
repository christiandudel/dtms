#' Calculate the distribution of the time until a subset of states is left for
#' the last time.
#'
#' @description
#' Calculates the distribution of the until a subset of states is left for the
#' very last time.
#'
#' @details
#' Details go here.
#'
#' @param matrix Matrix with transition probabilities, as generated with `dtms_matrix`.
#' @param dtms DTMS object as created with `dtms`.
#' @param risk Character (required), names of one or several states for which the waiting time should be calculated.
#' @param risk_to Character (optional), names of one or several states to which the states specified in `risk` are left. See details.
#' @param start_state Character (optional), names of one or several starting states. If NULL (default), all transient states will be considered.
#' @param start_time Numeric (optional), value of time scale at start. If NULL (default), first value of time scale is used.
#' @param end_time Numeric (optional), value of time scale at end. If NULL (default), last value of time scale is used.
#' @param start_distr Numeric (optional), distribution of starting states. Needs to be consistent with starting states. If provided, average distribution is provided; see details.
#' @param method Character, do transitions happen mid-interval (`mid`, default) or at the end of the interval (`end`), see details.
#' @param rescale Logical, should distribution be rescaled to sum to 1? See details. Default is TRUE.
#' @param total Logical, should total of distribution be shown? Default is FALSE, as the total always is 1.
#' @param transient Character (optional), short names of transient states. If NULL (default) transient states are taken from `dtms` object.
#' @param timescale Numeric (optional), values of time scale. If NULL (default) obtained from `dtms` object.
#' @param timestep Numeric (optional), step length of time scale. If NULL (default) obtained from `dtms` object.
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#'
#' @return The distribution of the waiting time until a subset of states is left for the last time.
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
#' ## First visit
#' dtms_last(dtms=simple,
#'            matrix=Tp,
#'            risk="A",
#'            start_distr=S)

dtms_last <- function(matrix,
                      dtms=NULL,
                      risk,
                      risk_to=NULL,
                      start_time=NULL,
                      start_state=NULL,
                      start_distr=NULL,
                      end_time=NULL,
                      transient=NULL,
                      timescale=NULL,
                      timestep=NULL,
                      method="mid",
                      sep="_",
                      total=TRUE,
                      rescale=TRUE){

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

  # Time scale: Only transitions starting up to T-1 relevant
  timescale <- timescale[-length(timescale)]

  # States of the transition matrix
  allstates <- rownames(matrix)
  nstates <- length(allstates)

  # Select subset
  selectorD <- dtms_in(allstates,risk,sep)
  if(is.null(risk_to)) selectorU <- !selectorD else selectorU <- dtms_in(allstates,risk_to,sep)

  # Get maxtime
  if(is.null(end_time)) maxtime <- length(timescale)-1 else maxtime <- which(end_time==timescale)

  # Generate t
  if(is.null(start_time)) t <- 0 else t <- which(start_time==timescale)-1

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
  results[["E_0"]] <- t(Biodem::mtx.exp(P_S,maxtime-t+1)%*%ones) *
                        Biodem::mtx.exp(matrix,t)
  colnames(results[["E_0"]]) <- allstates

  # Loop for other E_x
  step <- 0.5
  past.steps <- c(0)

  while(!step>=(maxtime-t+1)) {
    past.steps <- c(past.steps,step)
    e <- step-0.5
    tmp <- vector("list",1)
    names(tmp) <- paste("E",step,sep="_")
    tmp[[paste("E",step,sep="_")]] <- t(Biodem::mtx.exp(matrix,e)%*%P_E%*%
                                        Biodem::mtx.exp(P_S,maxtime-t-e)%*%ones) *
                                        Biodem::mtx.exp(matrix,t)
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
  steps <- past.steps*timestep
  if(method=="end") steps[-1] <- steps[-1]+0.5*timestep
  colnames(result) <- paste(steps)

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

  # Output
  return(result)

}
