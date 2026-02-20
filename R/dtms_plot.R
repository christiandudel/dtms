#' Plotting transition probabilities
#'
#' @description
#' A simple function for plotting transition probabilities with base R. This is
#' fast, but it is much easier to produce nicer looking results with
#' dtms_simplify.
#'
#' @param probs Object with transition probabilities as created with \code{dtms_transitions}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable in `probs` with starting state. Default is `from`.
#' @param tovar Character (optional), name of variable in `probs` with receiving state. Default is `to`.
#' @param timevar Character (optional), name of variable in `probs` with time scale. Default is `time`.
#' @param Pvar Character (optional), name of variable in `probs` with transition probabilities. Default is `P`.
#' @param ... Further arguments passed to plot().
#'
#' @return No return value, called for side effects
#' @export
#'
#' @examples
#' ## Model setup
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
#' ## Plot
#' dtms_plot(probs=probs,
#'           dtms=simple)

dtms_plot <- function(probs,
                      dtms,
                      fromvar="from",
                      tovar="to",
                      timevar="time",
                      Pvar="P",
                      ...) {

  # Check
  if(!is.null(dtms)) dtms_proper(dtms)

  # Simplify data
  probs <- dtms_simplify(probs,
                         fromvar=fromvar,
                         tovar=tovar,
                         sep=dtms$sep)

  # All states
  if(!is.null(dtms)) {
    transient <- dtms$transient
    allstates <- c(dtms$transient,dtms$absorbing)
  } else {
    transient <- unique(probs[,fromvar])
    allstates <- unique(probs[,tovar])
  }

  # Number of rows and columns
  nrows <- length(transient)
  ncols <- length(allstates)

  # Number of panels for layout
  npanels <- ncols*nrows

  # Layout
  graphics::layout(matrix(1:npanels,nrow=nrows))

  # Max probs
  maxprob <- max(probs[,Pvar])
  maxprob <- round(maxprob,digits=1)

  # Plot
  for(fromstate in transient) {

    for(tostate in allstates) {

      tmp <- subset(probs,
                    probs[,fromvar]%in%fromstate &
                      probs[,tovar]%in%tostate)

      plot(x=tmp[,timevar],
           y=tmp[,Pvar],
           type="l",
           ylim=c(0,maxprob),
           main=paste(fromstate,"to",tostate),
           ...)

    }

  }

  # Back to normal layout
  graphics::layout(1)

}
