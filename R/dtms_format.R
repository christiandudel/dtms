#' Reshape data to transition format
#'
#' @description
#' Takes a data frame in long format and reshapes it into transition format.
#'
#' @details
#' The data frame supplied with the `data` argument has to be in long format:
#'
#' \tabular{lll}{
#' idvar \tab timevar \tab statevar \cr
#' 1 \tab 0 \tab A \cr
#' 1 \tab 1 \tab A \cr
#' 1 \tab 2 \tab B \cr
#' 1 \tab 3 \tab A \cr
#' 2 \tab 0 \tab B \cr
#' 2 \tab 1 \tab A \cr
#' ... \tab ... \tab ...
#' }
#'
#' If it is not in long format it has to be reshaped. `dtms_format` turns this into
#' a data frame in transition format:
#'
#' \tabular{llll}{
#' id \tab time \tab fromvar \tab tovar \cr
#' 1 \tab 0 \tab A \tab A \cr
#' 1 \tab 1 \tab A \tab B \cr
#' 1 \tab 2 \tab B \tab A \cr
#' 2 \tab 0 \tab B \tab A \cr
#' ... \tab ... \tab ... \tab ... \cr
#' }
#'
#' By default the variable names of the ID variable and the time variable are changed to `id` and `time`, as the other
#' functions of the package use these as default names. If renaming the variables is not possible because
#' these variable names already appear in the data then the original names are used.
#'
#' `dtms_format` by default drops gaps in the data, as no transitions are observed. For instance,
#' in the following example there is no observation at time 4, and thus no transition
#' is observed from t=3 to t=4; and no transition from t=4 to t=5:
#'
#' \tabular{lll}{
#' idvar \tab timevar \tab statevar \cr
#' 1 \tab 0 \tab A \cr
#' 1 \tab 1 \tab A \cr
#' 1 \tab 2 \tab B \cr
#' 1 \tab 3 \tab A \cr
#' 1 \tab 5 \tab A \cr
#' }
#'
#' In this example, `dtms_format` will return the following:
#'
#'  \tabular{llll}{
#' id \tab time \tab fromvar \tab tovar \cr
#' 1 \tab 0 \tab A \tab A \cr
#' 1 \tab 1 \tab A \tab B \cr
#' 1 \tab 2 \tab B \tab A \cr
#' }
#'
#' If `fill=T`, then `dtms_format` will return the following:
#'
#'  \tabular{llll}{
#' id \tab time \tab fromvar \tab tovar \cr
#' 1 \tab 0 \tab A \tab A \cr
#' 1 \tab 1 \tab A \tab B \cr
#' 1 \tab 2 \tab B \tab A \cr
#' 1 \tab 3 \tab A \tab NA \cr
#' 1 \tab 4 \tab NA \tab A \cr
#' }
#'
#' @param data Data frame in long format.
#' @param dtms Object of class 'dtms' (optional) as created with `dtms`. Needed when arguments 'timestep' and/or `timescale' are not specified. Default is NULL.
#' @param idvar Character (optional), name of variable in 'data' with unit identifier, default is 'id'.
#' @param timevar Character (optional), name of variable in 'data' with time scale, default is 'time'.
#' @param statevar Character (optional), name of variable in 'data' with state, default is 'state'.
#' @param fromvar Character (optional), name of variable in reshaped data with starting state, default is 'from'.
#' @param tovar Character (optional), name of variable in reshaped data with receiving state, default is 'to'.
#' @param timescale Numeric (optional), values of time scale. Needed when argument `dtms` is not specified.
#' @param timestep Numeric (optional), step length for time scale. Needed when argument 'dtms' is not specified. Default is NULL.
#' @param keepnames Logical, keep original names for id and time variable? Default is FALSE; i.e., not keeping original names.
#' @param fill Logical, fill implicit missings with NA? Default is FALSE.
#' @param verbose Logical, create output to console if changing variable names is not possible? Default is TRUE.
#'
#' @return A data set reshaped to transition format
#' @export
#'
#' @examples
#' # Define model
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:20)
#' # Transiton format
#' estdata <- dtms_format(data=simpledata,
#' dtms=simple,
#' idvar="id",
#' timevar="time",
#' statevar="state")

dtms_format <- function(data,
                        dtms=NULL,
                        idvar="id",
                        timevar="time",
                        statevar="state",
                        fromvar="from",
                        tovar="to",
                        timestep=NULL,
                        timescale=NULL,
                        keepnames=F,
                        fill=F,
                        verbose=T) {

  # Use dtms if provided
  if(!is.null(dtms)) {

    # Check
    proper_dtms(dtms)

    # Use values
    timescale <- dtms$timescale
    timestep <- dtms$timestep
  }

  # Fill data
  if(fill) {
    # Copy time to expand (complete does not work well with sym)
    data <- data |> dplyr::mutate(temporary_Time = !!dplyr::sym(timevar))
    # Expand
    data <- data |> tidyr::complete(!!dplyr::sym(idvar),temporary_Time=timescale)
    # Drop original time var and rename
    data <- data |>
      dplyr::select(!(!!dplyr::sym(timevar))) |>
      dplyr::rename_with(~ c(timevar),c("temporary_Time"))
  }

  # Rearrange data into transition format
  res <- data  |>
    # Sort data by ID and time
    dplyr::arrange(!!dplyr::sym(idvar),!!dplyr::sym(timevar)) |>
    # Group by ID
    dplyr::group_by(!!dplyr::sym(idvar)) |>
    # Check if observations are consecutive
    dplyr::mutate(consec=dplyr::lead(!!dplyr::sym(timevar)),
                  difz  =consec-!!dplyr::sym(timevar),
                  # If consecutive, take leading state
                  to    =ifelse(difz==timestep,
                            dplyr::lead(!!dplyr::sym(statevar)),
                            NA)) |>
    # Ungroup
    dplyr::ungroup() |>
    # Drop added variables
    dplyr::select(!(consec:difz))

  # Rename state variable
  res <- res |> dplyr::rename_with(~ c(fromvar,tovar),c(statevar,"to"))

  # Change names of id and time to default
  if(!keepnames) {
    if(timevar!='time' & !'time'%in%names(res)) res <- res |> dplyr::rename('time' = timevar) else
      if(verbose) cat("Kept original name for time \n")
    if(idvar!='id' & !'id'%in%names(res)) res <- res |> dplyr::rename('id'   = idvar) else
      if(verbose) cat("Kept original name for id \n")
  }

  # Return result
  return(res)

}
