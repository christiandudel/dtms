#' Calculate delta
#'
#' @description
#' Calculates delta, either to compare transition probabilities from two
#' different models, or to assess how including lags changes transition
#' probabilities.
#'
#' @details
#' Delta is the weighted average absolute difference between the predicted
#' transition probabilities from two multistate models. It can attain values
#' between 0 and 1, where 0 indicates perfect similarity and 1 indicates that
#' the two models always give predictions at the opposite extremes; i.e., for
#' all predicted probabilities, one model predicts a probability of 0 and the
#' other predicts a probability of 1.
#'
#' This function is designed to use delta to assess the impact of including
#' different lags of the state variable in the model.
#'
#' To compare two different models, the arguments `data`, `model1`, and `model2`
#' are needed. `data` specifies the data frame used for predicting transition
#' probabilities. It needs to have all variables required for predicting based
#' on both `model1` and `model2`. The latter two arguments are the names of
#' multistate models estimated with \code{dtms_fit}.
#'
#' To compare how the inclusion of different lags of the state variable affects
#' predictions, a model needs to be specified using `data` and `dtms`, as well
#' as potential covariates with `controls`. The argument `lags` sets which lags
#' are included These are always including lower lags; e.g., a model including
#' the state at t-3 also has the state at t-2, at t-1, and at t. All resulting
#' models are compared to a model which does not control for the current or any
#' past state. If `lags=NULL` the Markov model is comapred to this model not
#' accounting for the current state or any past states.
#'
#' The argument `keepNA` controls how missing values are handled. These will
#' often occur for lagged states. For instance, for the first transition
#' observed for an individual, the state at time t is known, but not at time
#' t-1. In this case, if a first-order lag is used, this observation could either
#' be dropped; or, a missing value of the state at time t-1 could be included
#' as a predictor. `keepNA=TRUE` will do the latter, while if `FALSE`, all
#' observations with missing states are dropped. This is done for all
#' models, irrespective of the lag, such that they are based on exactly the
#' same observations.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param model1 Name of object containing a model estimated with \code{dtms_fit}.
#' @param model2 Name of object containing a model estimated with \code{dtms_fit}.
#' @param lags Numeric (optional), vector containing the lags as positive integers.
#' @param keepNA Logical (optional), keep missing values of lags as predictor value? Default is TRUE.
#' @param controls Character (optional), names of control variables
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable in `data` with receiving state. Default is "to".
#' @param idvar Character (optional), name of variable in `data` with unit ID. Default is "id".
#' @param timevar Character (optional), name of variable in `data` with time scale. Default is "time".
#' @param full Logical (optional), estimate fully interacted model? Default is FALSE.
#' @param package Character, chooses package for multinomial logistic regression, currently `VGAM`, `nnet`, and `mclogit` are supported. Default is `VGAM`.
#' @param reference Numeric or character (optional). Reference level of multinomial logistic regression.
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
#' ## Fit models
#' fit1 <- dtms_fit(data=estdata,controls="time")
#' fit2 <- dtms_fullfit(data=estdata,controls="time")
#'
#' ## Compare
#' dtms_delta(data=estdata,model1=fit1,model2=fit2)

dtms_delta <- function(data,
                       dtms=NULL,
                       model1=NULL,
                       model2=NULL,
                       lags=1:5,
                       controls=NULL,
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
  if(!keepNA) data <- stats::na.omit(data)

  # Number of models
  nmodels <- length(formulist)

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
      fitlist[[count]] <- VGAM::vgam(formula=model,
                                     family=VGAM::multinomial(refLevel=reference),
                                     data=data,
                                     ...)
    }

    #nnet
    if(package=="nnet") {

      # Estimate
      fitlist[[count]] <- nnet::multinom(formula=model,
                                         data=data,
                                         ...)
    }

    #mclogit
    if(package=="mclogit") {

      # Estimate
      fitlist[[count]] <- mclogit::mblogit(formula=model,
                                           data=data,
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

  # Names
  names(result) <- c("Markov",paste("Lag",lags))

  # Return results
  return(result)

}
