#' Estimate constrained discrete-time multistate model
#'
#' @description
#' This function estimates a constrained discrete-time multistate model
#' using multinomial logistic regression.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param controls Character (optional), names of control variables
#' @param fromvar Character (optional), name of variable with starting state. Default is "from".
#' @param tovar Character (optional), name of variable with receiving state. Default is "to".
#' @param timevar Character (optional), name of variable with time scale. Default is "time".
#' @param formula Formula (optional). If no formula is specified, it will be build from the information specified with controls, fromvar, tovar, and timevar.
#' @param weights Character (optional). Name of variable with survey weights.
#' @param reference Numeric or character (optional). Reference level of multinomial logistic regression.
#' @param package Character, chooses package for multinomial logistic regression, currently `VGAM`, `nnet`, and `mclogit` are supported. Default is `VGAM`.
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

dtms_fit <- function(data,
                     controls=NULL,
                     formula=NULL,
                     weights=NULL,
                     fromvar="from",
                     tovar="to",
                     timevar="time",
                     reference=1,
                     package="VGAM",
                     ...) {

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

  # Get weights if specified
  if(!is.null(weights)) weights <- data[,weights]

  # Factors (needed by most packages)
  data[,fromvar] <- as.factor(data[,fromvar])
  data[,tovar] <- as.factor(data[,tovar])
  data[,tovar] <- stats::relevel(data[,tovar],ref=reference)

  # VGAM
  if(package=="VGAM") {

    # Estimate
    dtms_fit <- VGAM::vgam(formula=formula,
                           family=VGAM::multinomial(refLevel=reference),
                           data=data,
                           weights=weights,
                           ...)
  }

  #nnet
  if(package=="nnet") {

    # Estimate
    dtms_fit <- nnet::multinom(formula=formula,
                               data=data,
                               weights=weights,
                               ...)

  }

  #mclogit
  if(package=="mclogit") {

    # Estimate
    dtms_fit <- mclogit::mblogit(formula=formula,
                                 data=data,
                                 weights=weights,
                                 ...)

  }


  # Return results
  return(dtms_fit)

}
