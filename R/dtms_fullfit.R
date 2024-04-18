#' Estimate unconstrained discrete-time multistate model
#'
#' @description
#' This function estimates an unconstrained discrete-time multistate model
#' using multinomial logistic regression. This is achieved by interacting
#' the starting state with all predictors in the model. It is a wrapper for
#' \code{dtms_fit()} with slightly less arguments.
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
#' fit <- dtms_fullfit(data=estdata)

dtms_fullfit <- function(data,
                         controls=NULL,
                         formula=NULL,
                         weights=NULL,
                         fromvar="from",
                         tovar="to",
                         timevar="time",
                         reference=1,
                         package="VGAM",
                         ...) {

  dtms_fit(data=data,
           controls=controls,
           formula=formula,
           weights=weights,
           fromvar=fromvar,
           tovar=tovar,
           timevar=timevar,
           reference=reference,
           package=package,
           full=TRUE)

}
