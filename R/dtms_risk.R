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
    times[!is.logical(times)] <- F
    selector <- selector & times
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
  whererisk <- !dtms$transient%in%risk
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
