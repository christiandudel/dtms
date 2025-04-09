#' Generate the reward matrix for a Markov chain with rewards
#'
#' @description
#' This function generates a reward matrix which can be used with
#' \code{dtms_reward}.
#'
#' @param dtms dtms object, as created with \code{dtms}.
#' @param starting Character (optional), name or names of starting states. If NULL (default) any transition to the state or states specififed with \code{receiving} will get the reward.
#' @param receiving Character, name or names of states to which transitioning generates the reward. Can be both transient or absorbing states.
#' @param reward Numeric, reward value to be placed in matrix.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#'
#' @return A matrix with rewards.
#' @export
#'
#' @examples
#' ## Define model: Absorbing and transient states, time scale
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:20)
#' dtms_rewardmatrix(dtms=simple,receiving="B",reward=0.3)

dtms_rewardmatrix <- function(dtms,
                              starting=NULL,
                              receiving,
                              reward,
                              start_time=NULL,
                              end_time=NULL) {

  # Check
  dtms_proper(dtms)

  # Combine states and time
  transient_states <- dtms_combine(dtms$transient,dtms$timescale,sep=dtms$sep)
  absorbing <- paste(dtms$absorbing)
  all_states <- c(transient_states,absorbing)

  # Starting states, long names
  if(is.null(starting)) starting <- dtms$transient
  starting_time <- dtms$timescale
  ntime <- length(starting_time)
  if(!is.null(start_time)) starting_time <- starting_time[which(dtms$timescale==start_time):ntime]
  if(!is.null(end_time)) starting_time <- starting_time[1:which(dtms$timescale==end_time)]
  starting_states <- dtms_combine(starting,
                                  starting_time,
                                  sep=dtms$sep)

  # Receiving states, time
  receiving_time <- dtms$timescale
  ntime <- length(receiving_time)
  if(!is.null(start_time)) receiving_time <- receiving_time[which(dtms$timescale==start_time):ntime]
  if(!is.null(end_time)) receiving_time <- receiving_time[1:which(dtms$timescale==end_time)]

  # Receiving states, split between transient and absorbing
  rec_tra <- receiving[receiving%in%dtms$transient]
  rec_abs <- receiving[receiving%in%dtms$absorbing]

  if(length(rec_tra)>0) {
    rec_tra <- dtms_combine(receiving,
                            receiving_time,
                            sep=dtms$sep)
  }

  receiving_states <- c(rec_tra,rec_abs)

  # Number of states
  nall <- length(all_states)

  # Build empty matrix
  result <- matrix(data=0,ncol=nall,nrow=nall)
  rownames(result) <- all_states
  colnames(result) <- all_states

  # Fill
  result[starting_states,receiving_states] <- reward

  # Return result
  return(result)

}
