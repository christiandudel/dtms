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
#' number of visits to that state.
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
#' @param start_distr Numeric (optional), distribution of starting states. If specified, average expectancy over all starting states will be calculated. Only applied if risk=NULL.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#' @param method Character (optional), do transitions happen mid-interval (`mid`, default) or at the end of the interval (`end`), see details.
#' @param rescale Logical (optional), should distribution be rescaled to sum to 1? See details. Default is TRUE.
#' @param total Logical, should total of distribution be shown? Default is FALSE, as the total always is 1.
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#'
#' @return Matrix with the distribution(s) of the waiting time.
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
                      dtms,
                      risk,
                      risk_to=NULL,
                      start_time=NULL,
                      start_state=NULL,
                      start_distr=NULL,
                      end_time=NULL,
                      method="mid",
                      sep="_",
                      total=TRUE,
                      rescale=TRUE){

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=sep)

  # Time scale: Only transitions starting up to T-1 relevant
  timescale_reduced <- dtms$timescale[-length(dtms$timescale)]

  # States of the transition matrix
  allstates <- rownames(matrix)
  nstates <- length(allstates)

  # Select subset
  selectorD <- dtms_in(allstates,risk,sep)
  if(is.null(risk_to)) selectorU <- !selectorD else
    selectorU <- dtms_in(allstates,risk_to,sep)

  # Get maxtime
  if(is.null(end_time)) maxtime <- length(timescale_reduced)-1 else
    maxtime <- which(end_time==timescale_reduced)

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
  steps <- past.steps*dtms$timestep
  if(method=="end") steps[-1] <- steps[-1]+0.5*dtms$timestep
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
