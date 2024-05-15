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
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param variables List (optional), a named list with covariate values which are used to restrict the data.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param fromvar Character (optional), name of variable with starting state. Default is `from`.
#' @param timevar Character (optional), name of variable with time scale. Default is `time`.
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
                       dtms,
                       variables=NULL,
                       start_state=NULL,
                       start_time=NULL,
                       fromvar="from",
                       timevar="time") {

  # Check
  dtms_proper(dtms)

  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)

  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=dtms$sep)

  # Restrict data: States and time
  data <- data[data[,fromvar]%in%start_state,]
  data <- data[data[,timevar]%in%start_time,]

  # Apply restrictions, if any
  varnames <- names(variables)
  for(var in varnames) {
    data <- data[data[,var]==variables[[var]],]
  }

  # Tabulate
  tmp <- data[,fromvar] |> table() |> prop.table()

  # Match with starting names
  result <- numeric(length(start_state))
  names(result) <- start_state
  matchnames <- names(tmp)[names(tmp)%in%start_state]
  result[matchnames] <- tmp[matchnames]

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
