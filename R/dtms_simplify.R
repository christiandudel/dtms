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
