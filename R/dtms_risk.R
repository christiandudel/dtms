#' Title
#'
#' @param matrix
#' @param transient
#' @param timescale
#' @param dtms
#' @param risk
#' @param start_time
#' @param start_state
#' @param start_distr
#' @param end_time
#' @param sep
#'
#' @return
#' @export
#'
#' @examples
dtms_risk <- function(matrix,# Matrix with transition probabilities generated with dtms_matrix
                      transient=NULL, # Names of transient states
                      timescale=NULL, # Time scale
                      dtms=NULL,# DTMS model
                      risk, # name of state(s) for which risk is of interest
                      start_time=NULL,# Starting time, will be lowest time in 'dtms' if NULL
                      start_state=NULL,# Names of starting states, will be all transient states in 'dtms' if NULL
                      start_distr=NULL,# Distribution of starting states for average
                      end_time=NULL, # Time up to which lifetime risks are calculated, making them partial lifetime risks
                      sep="_") {# Separator for names

  # Use dtms if provided
  if(!is.null(dtms) & class(dtms)[2]=="dtms") {
    if(is.null(transient)) transient <- dtms$transient
    if(is.null(timescale)) {
      timescale <- dtms$timescale
      timescale <- timescale[-length(timescale)]
    }
  }

  # Starting states
  if(is.null(start_state)) {
    starting <- levels(interaction(transient,min(timescale),sep=sep))
  } else {
    starting <- levels(interaction(start_state,start_time,sep=sep))
  }

  # States of the transition matrix
  states <- rownames(matrix)

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

  # Which are absorbing
  whichabsorbing <- which(diag(newmatrix)==1)

  # Umat
  Umat <- newmatrix[-whichabsorbing,-whichabsorbing]

  # Get N
  nstates <- dim(Umat)[1]
  Nmat <- solve(diag(1,nstates)-Umat)

  # Get R
  R <- newmatrix[-whichabsorbing,whichabsorbing]

  # Get results
  results <- Nmat%*%R

  # Get results
  result <- rep(1,length(starting))
  whererisk <- !transient%in%risk
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

}  # End of function
