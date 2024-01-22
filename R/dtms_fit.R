#' Title
#'
#' @param fromvar
#' @param tovar
#' @param timevar
#' @param data
#' @param controls
#' @param formula
#' @param weights
#' @param reference
#' @param method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
dtms_fit <- function(fromvar="from", # Name of variable with starting state
                     tovar="to",   # Name of variable with receiving state
                     timevar="time", # Name of variable with time
                     data,      # Name of data frame, required
                     controls=NULL, # Name(s) of control variables, if any
                     formula=NULL,  # Alternatively, full formula
                     weights=NULL,  # Weights, if any
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
    dtms <- VGAM::vgam(formula=formula,
                       family=VGAM::multinomial(refLevel=reference),
                       data=data,
                       weights=weights,
                       ...)
  }

  # Return results
  return(dtms)

}
