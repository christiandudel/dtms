### Create transition matrix without absorbing states
dtms_absorbing <- function(matrix) { # matrix=full transition matrix

  # Get states which are absorbing
  to_remove <- which(diag(matrix)==1)

  # Remove from matrix
  removed <- matrix[-to_remove,-to_remove]

  # Output reduced matrix
  return(removed)

}

### Check if object is proper dtms object
dtms_proper <- function(dtms) { # dtms=object to be checked

  # Error message
  message <- "Not a proper dtms object."

  # Check class
  if(!class(dtms)[2]=="dtms") stop(message)

  # Check names
  listnames <- c("transient" ,"absorbing" ,"timescale" ,"timestep" )
  if(!all(names(dtms)==listnames)) stop(message)

  # Check types (and length)
  if(!is.character(dtms$transient)) stop(message)
  if(!is.character(dtms$absorbing)) stop(message)
  if(!is.numeric(dtms$timescale) & length(dtms$timescale)>2) stop(message)
  if(!is.numeric(dtms$timestep) & length(dtms$timestep)==1) stop(message)

}

### Check if values of time scale are consecutive
dtms_consecutive <- function(data,idvar,timevar,timestep) {

  # Make sure no missing times and ids
  if(any(is.na(data[,timevar]))) stop("Missing values in time variable not allowed")
  if(any(is.na(data[,idvar]))) stop("Missing values in ID variable not allowed")

  # Get diff to next time step
  consecutive <- by(data[,timevar],data[,idvar],FUN=diff)

  # Add last obs
  consecutive <- lapply(consecutive,function(x) c(x,-1))

  # Unlist
  consecutive <- unlist(consecutive)

  # TRUE if equal to timestep, FALSE otherwise
  consecutive <- consecutive==timestep

  # Return
  return(consecutive)

}

### Rename variables
dtms_rename <- function(data,oldnames,newnames) {

  # Which names to change?
  changenames <- match(oldnames,names(data))

  # Change names
  names(data)[changenames] <- newnames

  # Return
  return(data)

}

### Combines values; e.g., to combine state names with time scale values
dtms_combine <- function(values1,values2,sep) {

  # Generate output vector
  output <- character(0)

  # Get values
  for(value in values1) {
    output <- c(output,paste(value,values2,sep=sep))
  }

  # Return
  return(output)

}
