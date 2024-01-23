#' Title
#'
#' @param matrix
#' @param transient
#' @param timescale
#' @param timestep
#' @param dtms
#' @param start_time
#' @param start_state
#' @param start_distr
#' @param end_time
#' @param sep
#' @param correction
#' @param total
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
dtms_expectancy <- function(matrix,# Matrix with transition probabilities generated with dtms_matrix
                            transient=NULL, # Names of transient states
                            timescale=NULL, # Time scale
                            timestep=NULL, # Step length of time scale
                            dtms=NULL,# DTMS model
                            start_time=NULL,# Starting time, will be lowest time in 'dtms' if NULL
                            start_state=NULL,# Names of starting states, will be all transient states in 'dtms' if NULL
                            start_distr=NULL,# Distribution of starting states for average
                            end_time=NULL, # Time up to which expectancies are calculated, for partial life expectancies
                            sep="_", # Separator for names
                            correction=0.5, # Correction; if set to 'NULL' no correction
                            total=T, # Should total expectancy be added?
                            verbose = F) { # Give output

  # Use dtms if provided
  if(!is.null(dtms) & class(dtms)[2]=="dtms") {
    if(is.null(transient)) transient <- dtms$transient
    if(is.null(timescale)) {
      timescale <- dtms$timescale
      timescale <- timescale[-length(timescale)]
    }
    if(is.null(timestep)) timestep <- dtms$timestep
  }

  # Starting states
  if(is.null(start_state)) {
    starting <- levels(interaction(transient,min(timescale),sep=sep))
  } else {
    starting <- levels(interaction(start_state,start_time,sep=sep))
  }

  # Number of starting and receiving states
  nstart <- length(starting)
  ntransient <- length(transient)

  # Remove absorbing states
  matrix <- remove_absorbing(matrix)

  # Fundamental matrix
  nstates <- dim(matrix)[1]
  Nmat <- solve(diag(1,nstates)-matrix)

  # Matrix for results
  result <- matrix(data=NA,ncol=ntransient,nrow=nstart)
  rownames(result) <- paste0("start:",starting)
  colnames(result) <- transient

  # Get and place state expectancies
  allstates <- rownames(matrix)

  for(i in 1:ntransient) {

    # Get states
    selector <- unlist(lapply(strsplit(allstates,split=sep),function(y) any(transient[i]%in%y) ))

    # Use end time
    if(!is.null(end_time)) {
      times <- as.numeric(unlist(lapply(strsplit(allstates,split=sep),function(y) y[2] )))
      times <- times<=end_time
      times[!is.logical(times)] <- F
      selector <- selector & times
    }

    # Calculate results and place
    if(nstart>1) tmp <- rowSums(Nmat[starting,selector]) else tmp <- sum(Nmat[starting,selector])

    # Place
    result[,transient[i]] <- tmp
  }

  # Correction
  if(is.numeric(correction)) {
    short_start <- unlist(lapply(strsplit(starting,split=sep),function(x) x[1]))
    correction_matrix <- outer(short_start,transient,FUN='==')
    correction_matrix <- correction_matrix*correction
    result <- result - correction_matrix
    if(verbose) cat("(Applying correction)","\n\n")
  }

  # Calculate average if starting distribution is provided
  if(!is.null(start_distr)) {
    AVERAGE <- colSums(result*start_distr)
    result <- rbind(result,AVERAGE)
  }

  # Add row totals
  if(total) {
    TOTAL <- rowSums(result)
    result <- cbind(result,TOTAL)
  }

  # Adjust for time step
  if(timestep!=1) {
    result <- result*timestep
    if(verbose) cat("Adjusting for step length","\n\n")
  }

  # Return result
  return(result)

}  # End of function
