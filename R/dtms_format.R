dtms_format <- function(data, # data frame
                        dtms=NULL, # DTMS object
                        idvar="id",   # Variable name of person ID
                        timevar="time", # Variable name of time
                        statevar="state", # Variable name of current state
                        fromvar="from", # New name for variable 'statevar' (from part)
                        tovar="to", # New name for variable 'statevar' (to part)
                        timestep=NULL, # Time step size; e.g., 2 if biannual, taken from DTMS if not specified
                        keepnames=F, # Keep original names for id and time variable?
                        verbose=T) {

  # Require tidy
  require(dplyr)

  # Get timestep from dtms
  if(is.null(timestep)) {
    if(is.null(dtms)) stop("No timestep specified")
    timestep <- dtms$timestep
  }

  # Rearrange data
  res <- data  |>
    # Group by person ID
    group_by(!!sym(idvar)) |>
    # Check if observations are consecutive
    mutate(consec=lead(!!sym(timevar)),
           difz = consec-!!sym(timevar),
           # If consecutive, take leading state
           to=ifelse(difz==timestep,lead(!!sym(statevar)),NA)) |>
    # Ungroup
    ungroup |>
    # Drop added variables
    select(!(consec:difz))

  # Rename state variable
  res <- res |> rename_with(~ c(fromvar,tovar),c(statevar,"to"))


  # Change names of id and time to default
  if(!keepnames) {
    if(timevar!='time' & !'time'%in%names(res)) res <- res |> rename('time' = timevar) else if(verbose) cat("Kept original name for time \n")
    if(idvar!='id' & !'id'%in%names(res)) res <- res |> rename('id'   = idvar) else if(verbose) cat("Kept original name for id \n")
  }

  # Return result
  return(res)

}
