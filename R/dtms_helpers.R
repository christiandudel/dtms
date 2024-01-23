### Create transition matrix without absorbing states

remove_absorbing <- function(matrix) { # matrix=full transition matrix

  to_remove <- which(diag(matrix)==1)
  removed <- matrix[-to_remove,-to_remove]
  return(removed)

}

### Check if object is proper dtms object

proper_dtms <- function(dtms) { # Object to be checked

  message <- "Not a proper dtms object"

  # Check class
  if(!class(dtms)[2]=="dtms") stop(message)

  # Check names
  listnames <- c("transient" ,"absorbing" ,"timescale" ,"timestep" )
  if(!all(names(dtms)==listnames))

  # Types
  if(!is.character(dtms$transient)) stop(message)
  if(!is.character(dtms$absorbing)) stop(message)
  if(!is.numeric(dtms$timescale) & length(dtms$timescale)>2) stop(message)
  if(!is.numeric(dtms$timestep) & length(dtms$timestep)==1) stop(message)

}
