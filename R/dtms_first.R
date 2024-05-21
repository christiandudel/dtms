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
#' number of visits to that state.
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
  for(i in 1:maxtime) {
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
  steps <- past_steps*dtms$timestep
  if(method=="end") steps[-1] <- steps[-1]+0.5*dtms$timestep
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
