#' Title
#'
#' @param transient
#' @param time
#' @param dtms
#' @param data
#' @param start_time
#' @param start_state
#' @param fromvar
#' @param timevar
#' @param variables
#' @param sep
#'
#' @return
#' @export
#'
#' @examples
dtms_start <- function(transient=NULL, # Names of transient states
                       time=NULL, # Time scale
                       dtms=NULL,# DTMS model
                       data,# data frame with cleaned transition data
                       start_time=NULL,# Starting time, will be lowest time in 'dtms' if NULL
                       start_state=NULL,# Names of starting states, will be all transient states in 'dtms' if NULL
                       fromvar="from", # Variable with starting state
                       timevar="time", # Variable name of time
                       variables=NULL, # List with variable values
                       sep="_") { # Separator for state names

  # Use dtms if provided
  if(!is.null(dtms) & class(dtms)[2]=="dtms") {
    if(is.null(transient)) transient <- dtms$transient
    if(is.null(time)) {
      time <- dtms$time
      time <- time[-length(time)]
    }
  }

  # Starting time
  if(is.null(start_time)) start_time <- min(time)

  # Starting states
  if(is.null(start_state)) {
    starting <- levels(interaction(transient,start_time,sep=sep))
  } else {
    starting <- levels(interaction(start_state,start_time,sep=sep))
  }

  # Restrict data: Time
  data <- data[data[,timevar]==start_time,]

  # Apply restrictions, if any
  nvar <- length(variables)
  varnames <- names(variables)
  if(nvar>0) {
    for(var in 1:nvar) {
      data <- data[data[,varnames[var]]==variables[[var]],]
    }
  }

  # Tabulate
  result <- prop.table(table(data[,fromvar]))

  # Match with starting names
  name_order <- match(names(result),transient)
  result <- as.numeric(result)[name_order]
  names(result) <- starting

  # Return
  return(result)

} #End of function
