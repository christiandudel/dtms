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

  # Get timestep from dtms
  if(is.null(timestep)) {
    if(is.null(dtms)) stop("No timestep specified")
    timestep <- dtms$timestep
  }

  # Rearrange data
  res <- data  |>
    # Group by person ID
    dplyr::group_by(!!dplyr::sym(idvar)) |>
    # Check if observations are consecutive
    dplyr::mutate(consec=dplyr::lead(!!dplyr::sym(timevar)),
           difz = consec-!!dplyr::sym(timevar),
           # If consecutive, take leading state
           to=ifelse(difz==timestep,dplyr::lead(!!dplyr::sym(statevar)),NA)) |>
    # Ungroup
    dplyr::ungroup() |>
    # Drop added variables
    dplyr::select(!(consec:difz))

  # Rename state variable
  res <- res |> dplyr::rename_with(~ c(fromvar,tovar),c(statevar,"to"))


  # Change names of id and time to default
  if(!keepnames) {
    if(timevar!='time' & !'time'%in%names(res)) res <- res |> dplyr::rename('time' = timevar) else if(verbose) cat("Kept original name for time \n")
    if(idvar!='id' & !'id'%in%names(res)) res <- res |> dplyr::rename('id'   = idvar) else if(verbose) cat("Kept original name for id \n")
  }

  # Return result
  return(res)

}
