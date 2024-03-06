#' Time needed to reach a subset of states for the first time
#'
#' @description
#' This function calculates the distribution of the time needed to reach a
#' subset of states for the first time.
#'
#' @details
#' The resulting distribution is conditional on ever reaching the subset of
#' states, as it is not a finite number if the set is never reached. If the
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
#' number of visits to that state.
#'
#' If a distribution of the starting states is provided with `start_distr` the
#' output table has two additional rows. One shows the distribution
#' unconditional on the starting state. The other shows the distribution
#' conditional on not starting in any state of the risk set.
#'
#' @param matrix Matrix with transition probabilities generated with dtms_matrix.
#' @param dtms DTMS object, as created with `dtms`.
#' @param risk Character, name of state(s) for which risk is of interest.
#' @param start_state Character (optional), names of one or several starting states. If NULL (default), all transient states will be considered.
#' @param start_time Numeric (optional), value of time scale at start. If NULL (default), first value of time scale is used.
#' @param end_time Numeric (optional), value of time scale at end. If NULL (default), last value of time scale is used.
#' @param start_distr Numeric (optional), distribution of starting states. Needs to be consistent with starting states. If provided, average distribution is provided; see details.
#' @param method Character, do transitions happen mid-interval (`mid`, default) or at the end of the interval (`end`), see details.
#' @param total Logical, should total of distribution be shown? Default is FALSE, as the total always is 1.
#' @param rescale Logical, should distribution be rescaled to sum to 1? See details. Default is TRUE.
#' @param transient Character (optional), short names of transient states. If NULL (default) transient states are taken from `dtms` object.
#' @param timescale Numeric (optional), values of time scale. If NULL (default) obtained from `dtms` object.
#' @param timestep Numeric (optional), step length of time scale. If NULL (default) obtained from `dtms` object.
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#'
#' @return A table of the distribution of the time needed to reach the subset of states
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
#' dtms_first(dtms=simple,
#'            matrix=Tp,
#'            risk="A",
#'            start_distr=S)

dtms_first  <- function(matrix,
                        dtms=NULL,
                        risk,
                        start_time=NULL,
                        start_state=NULL,
                        start_distr=NULL,
                        end_time=NULL,
                        transient=NULL,
                        timescale=NULL,
                        timestep=NULL,
                        method="mid",
                        total=TRUE,
                        rescale=TRUE,
                        sep="_") {

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
  selectorU <- dtms_in(allstates,risk,sep)

  # Generate matrix
  P_E <- matrix(data=0,nrow=nstates,ncol=nstates)
  P_E[!selectorU,!selectorU] <- matrix[!selectorU,!selectorU]

  # Generate t
  if(is.null(start_time)) t <- 0 else t <- which(start_time==timescale)-1

  # Generate max
  if(is.null(end_time)) maxtime <- length(timescale)-1 else maxtime <- which(end_time==timescale)-1

  # Generate W_t_0 and W_t_0.5, initial conditions
  results <- vector("list",2)
  names(results) <- c("W_0","W_0.5")
  results[["W_0"]] <- Biodem::mtx.exp(matrix,t)
  if(t==0) {
    results[["W_0.5"]] <- matrix(data=0,ncol=nstates,nrow=nstates)
    diag(results[["W_0.5"]][!selectorU,!selectorU]) <- 1
  } else results[["W_0.5"]] <- Biodem::mtx.exp(P_E,t)

  # Variables
  upcoming <- 1.5
  past_steps <- c(0,0.5)
  steps <- 1
  end <- 0

  # Loop to generate results
  for(i in 1:maxtime) {
    past_steps <- c(past_steps,upcoming)
    tmp <- vector("list",1)
    tmp[[1]] <- t(Biodem::mtx.exp(P_E,upcoming-0.5)%*%
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
  steps <- past_steps*timestep
  if(method=="end") steps[-1] <- steps[-1]+0.5*timestep
  colnames(result) <- paste(steps)[-length(steps)]

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
