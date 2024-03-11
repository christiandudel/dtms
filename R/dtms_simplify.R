#' Simplify state names
#'
#' @param probs Object with transition probabilities
#' @param fromvar Name of variable with starting state
#' @param tovar Name of variable with receiving state
#' @param sep Separator for long state name
#'
#' @return Data frame
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
#' ## Predict probabilities
#' probs    <- dtms_transitions(dtms=simple,
#'                              model = fit)
#' ## Simplify
#' probs |>  dtms_simplify()
#' ## NOT RUN: requires ggplot2
#' # library(ggplot2)
#' # probs |>  dtms_simplify() |>
#' #   ggplot(aes(x=time,y=P,color=to)) +
#' #   geom_line() +
#' #   facet_wrap(~from)

dtms_simplify <- function(probs,
                          fromvar="from",
                          tovar="to",
                          sep="_") {

  # Simplify names
  probs[,fromvar] <- dtms_getstate(probs[,fromvar],sep=sep)
  probs[,tovar] <- dtms_getstate(probs[,tovar],sep=sep)

  # Return
  return(probs)

}
