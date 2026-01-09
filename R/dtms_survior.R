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
#' @param probs Data frame with transition probabilities, as created with \code{dtms_transitions}.
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
#' ## Get starting distribution
#' S <- dtms_start(dtms=simple,
#'                 data=estdata)
#' ## Distribution of visits
#' dtms_survivor(dtms=simple,
#'               probs=probs,
#'               start_distr=S)


dtms_survivor <- function(probs=NULL,
                          matrix=NULL,
                          dtms,
                          start_distr=NULL,
                          start_state=NULL,
                          start_time=NULL,
                          end_time=NULL) {

  # Check
  dtms_proper(dtms)

  # Get matrix if not specified
  if(is.null(matrix)) matrix <- dtms_matrix(probs=probs,
                                            dtms=dtms)

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
