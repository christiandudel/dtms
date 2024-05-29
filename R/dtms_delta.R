#' Calculate delta
#'
#' @description
#' Calculates delta, either to compare transition probabilities from two
#' different models, or to assess how including lags changes transition
#' probabilites.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param controls Character (optional), names of control variables
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable in `data` with receiving state. Default is "to".
#' @param formula Formula (optional). If no formula is specified, it will be build from the information specified with controls, fromvar, tovar, and timevar.
#' @param full Logical (optional), estimate fully interacted model? Default is FALSE.
#' @param package Character, chooses package for multinomial logistic regression, currently `VGAM`, `nnet`, and `mclogit` are supported. Default is `VGAM`.
#' @param reference Numeric or character (optional). Reference level of multinomial logistic regression.
#' @param weights Character (optional). Name of variable with survey weights.
#' @param ... Further arguments passed to estimation functions.
#'
#' @return Vector of values of delta
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
#' fit <- dtms_delta(data=estdata)

dtms_delta <- function(data,
                       dtms=NULL,
                       model1=NULL,
                       model2=NULL,
                       lags=1:5,
                       controls=NULL,
                       weights=NULL,
                       fromvar="from",
                       tovar="to",
                       timevar="time",
                       idvar="id",
                       reference=1,
                       package="VGAM",
                       full=FALSE,
                       keepNA=TRUE,
                       ...) {

  # Compare models
  if(!is.null(model1) & !is.null(model2)) {

    # Predict based on model class
    if(inherits(model1,c("vgam","mclogit"))) {
      pred1 <- stats::predict(model1,data,"response")
    }

    if(inherits(model2,c("vgam","mclogit"))) {
      pred2 <- stats::predict(model2,data,"response")
    }

    if(inherits(model1,"nnet")) {
      pred1 <- stats::predict(model1,data,"probs")
    }

    if(inherits(model2,"nnet")) {
      pred2 <- stats::predict(model2,data,"probs")
    }

    # Calculate delta
    delta <- mean(rowSums(0.5*abs(pred1-pred2)))

    # Return
    return(delta)

  }

  # List for formulas
  formulist <- list()

  # Model without history
  formulist[[1]] <- dtms_formula(controls=controls,
                                 fromvar=NULL,
                                 tovar=tovar,
                                 full=full)

  # Markov model
  formulist[[2]] <- dtms_formula(controls=controls,
                                 fromvar=fromvar,
                                 tovar=tovar,
                                 full=full)

  # Loop over lags
  for(addlag in lags) {

    # Name of lag variable
    varname <- paste0(fromvar,"l",addlag)

    # Add to data
    data[,varname] <- dtms_lag(data=data,
                               dtms=dtms,
                               lag=addlag,
                               fromvar=fromvar,
                               idvar=idvar,
                               timevar=timevar)

    # Keep NA?
    if(keepNA) data[is.na(data[,varname]),varname] <- "NA"

    # As factor
    data[,varname] <- as.factor(data[,varname])

    # Add lag to controls
    controls <- c(controls,varname)

    # Generate formula
    nlag <- which(lags==addlag)
    formulist[[2+nlag]] <- dtms_formula(controls=controls,
                                      fromvar=fromvar,
                                      tovar=tovar,
                                      full=full)

  }

  # Keep NA?
  if(!keepNA) data <- na.omit(data)

  # Number of models
  nmodels <- length(formulist)

  # Get weights if specified
  if(!is.null(weights)) weights <- data[,weights]

  # Factors (needed by some packages)
  data[,fromvar] <- as.factor(data[,fromvar])
  data[,tovar] <- as.factor(data[,tovar])
  data[,tovar] <- stats::relevel(data[,tovar],ref=reference)

  # List for regression results
  fitlist <- list()

  # Making list assignment below simple
  count <- 1

  # Loop over models
  for(model in formulist) {

    # VGAM
    if(package=="VGAM") {

      # Estimate
      fitlist[count] <- VGAM::vgam(formula=model,
                                   family=VGAM::multinomial(refLevel=reference),
                                   data=data,
                                   weights=weights,
                                   ...)
    }

    #nnet
    if(package=="nnet") {

      # Estimate
      fitlist[count] <- nnet::multinom(formula=model,
                                       data=data,
                                       weights=weights,
                                       ...)

    }

    #mclogit
    if(package=="mclogit") {

      # Estimate
      fitlist[count] <- mclogit::mblogit(formula=model,
                                         data=data,
                                         weights=weights,
                                         ...)

    }

    count <- count + 1

  }

  # List for predicted results
  predictlist <- list()
  count <- 1

  # Predict
  for(fit in fitlist) {

    if(inherits(fit,c("vgam","mclogit"))) {
      predictlist[[count]] <- stats::predict(fit,data,"response")
    }

    if(inherits(fit,"nnet")) {
      predictlist[[count]] <- stats::predict(fit,data,"probs")
    }

    count <- count + 1

  }

  # Get differences
  resultlength <- nmodels-1
  result <- numeric(resultlength)

  for(i in 1:resultlength) {

    result[i] <- mean(rowSums(abs(predictlist[[1]] - predictlist[[i+1]])))

  }

  # Return results
  return(fitted)

}
