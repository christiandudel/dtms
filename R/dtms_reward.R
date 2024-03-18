#' Markov chain with rewards
#'
#' @description
#' This function calculates the expected rewards by starting state in a
#' Markov chain with rewards.
#'
#' @param matrix Matrix with transition probabilities, as generated with \code{dtms_matrix}.
#' @param reward Matrix with rewards, has to be of same dimensions as `matrix`.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#'
#' @return A matrix with expected rewards.
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
#' ## Reward matrix
#' Rw <- diag(1,dim(Tp)[1])
#' ## State expectancies
#' dtms_reward(dtms=simple,
#'             matrix=Tp,
#'             reward=Rw)

dtms_reward <- function(matrix,
                        reward,
                        dtms,
                        sep="_") {

  # Check
  dtms_proper(dtms)

  # Starting states, long names
  starting <- dtms_combine(dtms$transient,
                           min(dtms$timescale),
                           sep=sep)

  # Number of starting and receiving states
  nstart <- length(starting)
  ntransient <- length(dtms$transient)
  nabsorbing <- length(dtms$absorbing)

  # Remove absorbing states
  nmatrix <- dtms_absorbing(matrix)

  # Fundamental matrix
  nstates <- dim(nmatrix)[1]
  Nmat <- solve(diag(1,nstates)-nmatrix)

  # Additional matrices
  Z <- cbind(diag(1,nstates),matrix(data=0,ncol=nabsorbing,nrow=nstates))
  one <- rep(1,nabsorbing+nstates)

  # e
  result <- Nmat %*% Z %*%(matrix %*% reward) %*% one
  colnames(result) <- "Reward"

  # Return result
  return(result)

}
