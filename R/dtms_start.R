#' Tabulate starting distribution
#'
#' @description
#' Tabulates the starting distribution.
#'
#' @details
#' Per default, the starting distribution is the distribution of transient
#' states at the first value of the time scale in the data. This can be
#' changed to any value of the time scale, and any set of states. The
#' distribution can also be conditional on further covariate values which can be
#' specified with the argument `variables`.
#'
#' `variables` takes a named list where each entry of the list is named like
#' the corresponding variable and with the values to be selected.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param variables List (optional), a named list with covariate values which are used to restrict the data.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If several values are specified, the average distribution over all these values is calculated. In this case the first value specified with this argument is used to construct the long state name. If NULL (default) first value of time scale will be used.
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is `from`.
#' @param timevar Character (optional), name of variable in `data` with time scale. Default is `time`.
#' @param weights Character (optional). Name of variable with survey weights.
#'
#' @return Returns a table of the starting distribution.
#' @export
#'
#' @examples
## Define model: Absorbing and transient states, time scale
#' work <- dtms(transient=c("Working","Non-working","Retired"),
#'              absorbing="Dead",
#'              timescale=50:99)
#' ## Reshape
#' estdata <- dtms_format(data=workdata,
#'                        dtms=work,
#'                        idvar="ID",
#'                        timevar="Age",
#'                        statevar="State")
#' ## Drop dead-to-dead transitions etc
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=work)
#' ## Starting distributions
#' # Men
#' Sm <- dtms_start(dtms=work,
#'                  data=estdata,
#'                  variables=list(Gender=0))
#' # Women
#' Sw <- dtms_start(dtms=work,
#'                  data=estdata,
#'                  variables=list(Gender=1))

dtms_start <- function(data,
                       dtms,
                       variables=NULL,
                       start_state=NULL,
                       start_time=NULL,
                       fromvar="from",
                       timevar="time",
                       weights=NULL) {

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time[1],sep=dtms$sep)

  # Restrict data: States and time
  data <- data[data[,fromvar]%in%start_state,]
  data <- data[data[,timevar]%in%start_time,]

  # Apply restrictions, if any
  varnames <- names(variables)
  for(var in varnames) {
    data <- data[data[,var]%in%variables[[var]],]
  }

  # Tabulate
  if(is.null(weights)) {
    tab <- data[,fromvar] |> table() |> prop.table()
  } else {
    # https://stackoverflow.com/questions/18585977/frequency-tables-with-weighted-data-in-r
    tmp <- stats::aggregate(x = data[,weights], by = list(data[,fromvar]), FUN = sum)
    tab <- tmp[,2]
    names(tab) <- tmp[,1]
    tab <- prop.table(tab)
  }

  # Match with starting names
  result <- numeric(length(start_state))
  names(result) <- start_state
  matchnames <- names(tab)[names(tab)%in%start_state]
  result[matchnames] <- tab[matchnames]

  # Fix if issues and procude warning
  wrong <- which(result%in%c(NA,NaN,Inf,-Inf))
  if(length(wrong)>0) {
    warning("Something might have gone wrong with the starting distribution")
    result[wrong] <- 0
    result <- result/sum(result)
  }

  # Proper names
  names(result) <- starting

  # Return
  return(result)

}
