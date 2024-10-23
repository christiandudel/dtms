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
  listnames <- c("transient" ,"absorbing" ,"timescale" ,"timestep" ,"sep")
  if(!all(names(dtms)==listnames)) stop(message)

  # Check types (and length)
  if(!is.character(dtms$transient)&!is.numeric(dtms$transient)) stop(message)
  if(!is.character(dtms$absorbing)&!is.numeric(dtms$absorbing)) stop(message)
  if(!is.character(dtms$sep)) stop(message)
  if(!is.numeric(dtms$timescale) | length(dtms$timescale)<2) stop(message)
  if(!is.numeric(dtms$timestep)) stop(message)

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
  result <- data.frame(true=consecutive%in%timestep,
                       numeric=consecutive)

  # Return
  return(result)

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

### Check if a short state name is in long name
dtms_in <- function(vector,name,sep) {
  res <- lapply(strsplit(vector,split=sep),function(y) any(name%in%y[[1]]) )
  res <- unlist(res)
  return(res)
}

### Get time from long name
dtms_gettime <- function(vector,sep) {
  res <- lapply(strsplit(vector,split=sep),function(y) y[2] )
  res <- unlist(res) |> as.numeric()
  return(res)
}

### Get state from long name
dtms_getstate <- function(vector,sep) {
  res <-lapply(strsplit(vector,split=sep),function(x) x[1])
  res <- unlist(res)
  return(res)
}

### Calculate power of matrix
dtms_mtexp <- function(matrix,n) {

  res <- diag(nrow=nrow(matrix))
  rep <- 0

  while(rep<n) {
    res <- res %*% matrix
    rep <- rep+1
  }

  return(res)

}

### Generate formula for estimation
dtms_formula <- function(controls, # Arguments the same as for dtms_fit
                         fromvar,
                         tovar,
                         full) {

  # If fromvar is NULL (not default, needs explicit call)
  if(is.null(fromvar)) fromvar <- "1"

  # Constrained model
  if(!full) {
    formula <- paste0(tovar,"~",fromvar)
    if(!is.null(controls)) {
      controls <- paste(controls,collapse="+")
      formula <- paste(formula,controls,sep="+")
    }
    formula <- stats::as.formula(formula)
    # Unconstrained/fully interacted model
  } else {
    formula <- paste0(tovar,"~",fromvar)
    if(!is.null(controls)) {
      varlist <- paste(controls,fromvar,sep="*")
      varlist <- paste(varlist,collapse="+")
      formula <- paste(formula,varlist,sep="+")
    }
    formula <- stats::as.formula(formula)
  }

  # Return
  return(formula)

}

### Get lagged state variable, accounting for gaps in the data
dtms_lag <- function(data,
                     dtms,
                     lag,
                     fromvar="from",
                     idvar="id",
                     timevar="time") {

  # Make data smaller
  data <- data[,c(idvar,timevar,fromvar)]

  # Get ID values
  idvalues <- data[,idvar] |> unique()

  # Full data
  fulldata <- expand.grid(dtms$timescale,idvalues,
                          stringsAsFactors=FALSE)
  names(fulldata) <- c(timevar,idvar)

  # Merge with data
  fulldata <- merge(fulldata,data,
                    by=c(idvar,timevar),
                    all=T)

  # shift state
  stateshift <- by(fulldata[,fromvar],
                   fulldata[,idvar],function(x) c(rep(NA,min(lag,length(x))),
                                                  x[-( max(0,(length(x)-(lag-1))): length(x))]))

  # shift time
  timeshift <- by(fulldata[,timevar],
                  fulldata[,idvar],function(x) c(rep(NA,min(lag,length(x))),
                                                 diff(x,lag=lag) ))

  # Unlist
  stateshift <- unlist(stateshift)
  timeshift <- unlist(timeshift)

  # drop wrong spacing
  stateshift[!timeshift%in%c(NA,dtms$timestep*lag)] <- NA

  # Merge back
  fulldata$stateshift <- stateshift
  data <- merge(data,fulldata,by=c(idvar,timevar,fromvar),sort=FALSE)

  # Return
  return(data$stateshift)

}

### Carry over absorbing values
dtms_carry <- function(x,
                       dtms) {

  # Check if action necessary
  if(any(x%in%dtms$absorbing)) {

    # Length
    n <- length(x)

    # From where to carry over
    whichfirst <- which(x%in%dtms$absorbing)[1]

    # What to carry over
    whichvalue <- x[whichfirst]

    # Carry over
    x[whichfirst:n] <- whichvalue

    # Return
    return(x)

  } else return(x)

}

### Move states forward
dtms_forward_help <- function(x, # Vector of states
                              state, # State name as character string
                              overwrite="missing", # transient, missing, absorbing, all
                              dtms=NULL) {

  # Length
  nx <- length(x)

  # Find appearances
  whichfirst <- which(x==state)

  # Stop if no appearance
  if(length(whichfirst)==0) return(x)

  # Get first
  whichfirst <- min(whichfirst)

  # Replace: both missing and absorbing
  if(overwrite=="all") {
    dochange <- whichfirst:nx
    x[dochange] <- state
  }

  # Replace: only missing
  if(overwrite=="missing") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    dochange <- whichfirst:nx
    dochange <- setdiff(dochange,whichabsorbing)
    x[dochange] <- state
  }

  # Replace: only absorbing
  if(overwrite=="absorbing") {
    whichmissing <- which(is.na(x))
    dochange <- whichfirst:nx
    dochange <- setdiff(dochange,whichmissing)
    x[dochange] <- state
  }

  # Do not replace absorbing or missing
  if(overwrite=="transient") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    whichmissing <- which(is.na(x))
    whichboth <- union(whichabsorbing,whichmissing)
    dochange <- whichfirst:nx
    dochange <- setdiff(dochange,whichboth)
    x[dochange] <- state
  }

  # Return
  return(x)

}

### Carry states backward
dtms_backward_help <- function(x, # Vector of states
                               state, # State name as character string
                               overwrite="missing", # transient, missing, absorbing, all
                               dtms=NULL) {

  # Find appearances
  whichfirst <- which(x==state)

  # Stop if no appearance
  if(length(whichfirst)==0) return(x)

  # Get first
  whichfirst <- min(whichfirst)

  # Replace: both missing and absorbing
  if(overwrite=="all") {
    dochange <- 1:whichfirst
    x[dochange] <- state
  }

  # Replace: only missing
  if(overwrite=="missing") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    dochange <- 1:whichfirst
    dochange <- setdiff(dochange,whichabsorbing)
    x[dochange] <- state
  }

  # Replace: only absorbing
  if(overwrite=="absorbing") {
    whichmissing <- which(is.na(x))
    dochange <- 1:whichfirst
    dochange <- setdiff(dochange,whichmissing)
    x[dochange] <- state
  }

  # Do not replace absorbing or missing
  if(overwrite=="transient") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    whichmissing <- which(is.na(x))
    whichboth <- union(whichabsorbing,whichmissing)
    dochange <- 1:whichfirst
    dochange <- setdiff(dochange,whichboth)
    x[dochange] <- state
  }

  # Return
  return(x)

}
