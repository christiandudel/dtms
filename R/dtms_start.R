#' Tabulate starting distribution
#'
#' @description
#' Tabulates the starting distribution.
#'
#' @details
#' Per default, the starting distribution is the distribution of transient
#' states #' at the first value of the time scale in the data. This can be
#' changed to any value of the time scale, and any set of states. The
#' distribution can also be conditional on further covariate values which can be
#' specified with the argument `variables`.
#'
#' @param data Data in transition format, as created with `dtms_format` and cleaned with `dtms_clean`.
#' @param dtms DTMS object.
#' @param variables List (optional), a named list with covariate values which are used to restrict the data further.
#' @param transient Character (optional), names of all transient states.
#' @param timescale Numeric (optional), values of the time scale.
#' @param start_state Character (optional), names of starting states. If NULL (default), all transient states will be considered.
#' @param start_time Numeric (optional), start time. If NULL (default), the first value of the time scale will be used.
#' @param fromvar Character, name of variable with starting state. Default is `from`, which is the package default for all functions.
#' @param timevar Character, name of variable with time scale. Default is `time`, which is the package default for all functions.
#' @param sep Character, separator for constructed state names. Default is `_`, which is the package default for all functions.
#'
#' @return Returns a table of the starting distribution.
#' @export
#'
#' @examples
## Define model: Absorbing and transient states, time scale
#' hrs <- dtms(transient=c("Working","Non-working","Retired"),
#'             absorbing="Dead",
#'             timescale=50:99)
#' ## Reshape
#' estdata <- dtms_format(data=hrsdata,
#'                        dtms=hrs,
#'                        idvar="ID",
#'                        timevar="Age",
#'                        statevar="State")
#' ## Drop dead-to-dead transitions etc
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=hrs)
#' ## Starting distributions
#' # Men
#' Sm <- dtms_start(dtms=hrs,
#'                  data=estdata,
#'                  variables=list(Gender=0))
#' # Women
#' Sw <- dtms_start(dtms=hrs,
#'                  data=estdata,
#'                  variables=list(Gender=1))

dtms_start <- function(data,
                       dtms=NULL,
                       variables=NULL,
                       transient=NULL,
                       timescale=NULL,
                       start_state=NULL,
                       start_time=NULL,
                       fromvar="from",
                       timevar="time",
                       sep="_") {

  # Use dtms if provided
  if(!is.null(dtms)) {

    # Check
    dtms_proper(dtms)

    # Use values
    transient <- dtms$transient
    timescale <- dtms$timescale
  }

  # Starting state and time
  if(is.null(start_state)) start_state <- transient
  if(is.null(start_time)) start_time <- min(timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=sep)

  # Restrict data: States and time
  data <- data[data[,fromvar]%in%start_state,]
  data <- data[data[,timevar]==start_time,]

  # Apply restrictions, if any
  varnames <- names(variables)
  for(var in varnames) {
    data <- data[data[,var]==variables[[var]],]
  }

  # Tabulate
  result <- data[,fromvar] |> table() |> prop.table()

  # Match with starting names
  name_order <- match(names(result),start_state)
  result <- as.numeric(result)[name_order]
  names(result) <- starting

  # Return
  return(result)

}
