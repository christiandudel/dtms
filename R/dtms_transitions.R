#' Predict transition probabilities
#'
#' @description
#' `dtms_transitions` predicts transition probabilities based on a model
#' estimated with `dtms_fit`.
#'
#' @details
#' Predicted transition probabilities are returned as a data frame, and not
#' as a transition matrix. While the latter is required for applying Markov
#' chain methods, the data frame is more convenient for viewing and
#' analyzing the transition probabilities themselves.
#'
#' Depending on the model specification, the prediction of transition
#' probabilities will require values for predictor variables which can be
#' specified with the arguments `constant` and `varying` for time-constant and
#' time-varying variables, respectively. In both cases, named lists have to be
#' used, where each entry name must correspond to a variable name in the model.
#' For time-constant variables, each list entry is of length one and provides
#' a value for the corresponding time-constant variable. For time-varying
#' variables, each entry must have the length of the time scale minus one, and
#' provide a value for each (potential) transition in the model; i.e., starting
#' from time t=0, starting from time t=1, etc., until time t=T-1.
#'
#' The argument `separator` sets the separator used to create state names. State
#' names are either a combination of the name of a transient state and a value
#' of the time scale, or the name of an absorbing state.
#'
#' @param model Model for transition probabilities estimated with `dtms_fit`.
#' @param dtms DTMS object as created with `dtms`.
#' @param constant List (optional) with values for time-constant predictors (see details).
#' @param varying List (optional) with values for time-varying predictors (see details).
#' @param transient Character (optional), names of transient states in the model.
#' @param absorbing Character (optional), names of absorbing states in the model.
#' @param timescale Numeric (optional), values of the time scale.
#' @param timestep Numeric (optional), step length of the time scale.
#' @param timevar Character, name of variable with time scale in output data. Default is `time`.
#' @param fromvar Character, name of variable with sending state in output data. Default is `from`.
#' @param tovar Character, name of variable with receiving state in output data. Default is `to`.
#' @param Pvar Character, name of variable with transition probabilities in output data. Default is `P`.
#' @param sep Character, separator between state name and value of time scale. Default is `_`, see details.
#'
#' @return A data frame with transition probabilities.
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

dtms_transitions <- function(model,
                             dtms=NULL,
                             constant=NULL,
                             varying=NULL,
                             transient=NULL,
                             absorbing=NULL,
                             timescale=NULL,
                             timestep=NULL,
                             timevar="time",
                             fromvar="from",
                             tovar="to",
                             Pvar="P",
                             sep="_") {

  # Use dtms if provided
  if(!is.null(dtms)) {

    # Check
    dtms_proper(dtms)

    # Use values
    timescale <- dtms$timescale
    absorbing <- dtms$absorbing
    transient <- dtms$transient
    timestep <- dtms$timestep
  }

  # Adjust time scale (transitions in the model)
  timescale <- timescale[-length(timescale)]

  # Get full state space
  all_states <- paste(c(transient,absorbing))

  # Create empty frame
  model_frame <- expand.grid(from=transient,
                             time=timescale)

  # Get names right
  names(model_frame) <- c(fromvar,timevar)

  # Add time-constant covariates
  nvar <- length(constant)
  if(nvar>0) {
    for(var in 1:nvar) {
      varname <- names(constant)[var] # Variable name
      value <- constant[[var]] # Value
      model_frame[varname] <- value # Assign
    }
  }

  # Add time-varying covariates
  nvar <- length(varying)
  if(nvar>0) {
    for(var in 1:nvar) {
      varname <- names(varying)[var] # Variable name
      value <- varying[[var]] # Value
      if(length(value)!=length(timescale)) stop("Wrong number of time-varying values")
      assign_values <- match(model_frame[,timevar],timescale) # Match to time variable
      model_frame[varname] <- value[assign_values] # Assign
    }
  }

  # Predict (might need adjustment for other packages)
  if(class(model)=="vgam") {
  model_frame[,all_states] <- stats::predict(model,model_frame,"response")[,all_states]
  } else stop("Currently only vgam is supported")

  # Values of starting state
  model_frame$from <- paste(model_frame$from,model_frame$time,sep=sep)

  # Reshape
  model_frame <- reshape(model_frame,
                   varying=all_states,
                   idvar=fromvar,
                   timevar=tovar,
                   times=all_states,
                   direction="long",
                   v.names=Pvar)

  # Values of receiving variable
  rightrows <- model_frame[,tovar]%in%transient
  oldvalues <- model_frame[rightrows,tovar]
  timevalues <- model_frame[rightrows,timevar]+timestep
  model_frame[rightrows,tovar] <- paste(oldvalues,
                                        timevalues,
                                        sep=sep)

  # Return
  return(model_frame)

}
