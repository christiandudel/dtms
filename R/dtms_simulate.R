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

  initial_state <- sample(starting,
                          size=1,
                          prob=start_distr)

  simseq <- markovchain::rmarkovchain(n = ntime-1,
                                      object = sim,
                                      t0 = initial_state,
                                      include.t0=T)

  simdata <- rbind(simdata,simseq)

  }

  # Nicer names
  names(simdata) <- paste0(varnames,dtms$timescale)

  # Return
  return(simdata)

}
