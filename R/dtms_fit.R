#' Estimate constrained discrete-time multistate model
#'
#' @description
#' This function estimates a constrained discrete-time multistate model
#' using multinomial logistic regression.
#'
#' @param data Data frame, as prepared with `dtms_format` and `dtms_clean`.
#' @param controls Character (optional), names of control variables
#' @param fromvar Character (optional), name of variable with starting state. Defaults to `from`, which is used as default in other functions of the package.
#' @param tovar Character (optional), name of variable with receiving state. Defaults to `to`, which is used as default in other functions of the package.
#' @param timevar Character (optional), name of variable with time scale. Defaults to `time`, which is used as default in other functions of the package.
#' @param formula Formula (optional). If no formula is specified, it will be build from the information specified with controls, fromvar, tovar, and timevar.
#' @param weights Character (optional). Name of variable with survey weights.
#' @param reference Numeric (optional). Reference level of multinomial logistic regression (VGAM)
#' @param method Character, chooses package for multinomial logistic regression, currently only `VGAM` supported.
#' @param ... Further arguments passed to estimation functions.
#'
#' @return Returns an object with class depending on the package used.
#' @export
#'
#' @examples
#' ## Define model: Absorbing and transient states, time scale
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:20)
#' ## Reshape to transition format
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' ## Clean
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' ## Fit model
#' fit <- dtms_fit(data=estdata)

dtms_fit <- function(data,      # Name of data frame, required
                     controls=NULL, # Name(s) of control variables, if any
                     formula=NULL,  # Alternatively, full formula
                     weights=NULL,  # Weights, if any
                     fromvar="from", # Name of variable with starting state
                     tovar="to",   # Name of variable with receiving state
                     timevar="time", # Name of variable with time
                     reference=1,   # Reference category
                     method="VGAM", # Function to use
                     ...) {  # Further arguments to be passed to estimation function

  # Checks
  if(!method=="VGAM") stop("Currently only VGAM supported")

  # Build formula if not specified
  if(is.null(formula)) {
    formula <- paste0(tovar,"~",fromvar)
    if(!is.null(timevar)) formula <- paste(formula,timevar,sep="+")
    if(!is.null(controls)) {
      controls <- paste(controls,collapse="+")
      formula <- paste(formula,controls,sep="+")
    }
    formula <- stats::as.formula(formula)
  }

  # Estimate
  if(method=="VGAM") {

    # Calculate
    dtms_fit <- VGAM::vgam(formula=formula,
                           family=VGAM::multinomial(refLevel=reference),
                           data=data,
                           weights=weights,
                           ...)
  }

  # Return results
  return(dtms_fit)

}
