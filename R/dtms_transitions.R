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
#' If `vcov=TRUE` the full variance-covariance matrix of the transition
#' probabilities will be returned instead of the transition probabilities.
#'
#' The argument `dropvar` controls whether the covariate values used for
#' prediction are dropped. If `FALSE` each row of the resulting data frame will
#' have the covariate values #' which were used to predict the corresponding
#' probability.
#'
#' The argument `separator` sets the separator used to create state names. State
#' names are either a combination of the name of a transient state and a value
#' of the time scale, or the name of an absorbing state.
#'
#' @param model Model estimated with \code{dtms_fit}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param constant List (optional) with values for time-constant predictors (see details).
#' @param varying List (optional) with values for time-varying predictors (see details).
#' @param se Logical (optional), return standard errors of predicted probabilites. Default is `TRUE`.
#' @param vcov Logical (optional), return variance-covariance matrix of predicted probabilites. Default is `FALSE`.
#' @param dropvar Logical (optional), should covariate values used for prediction be returned (see details). Default is `TRUE`.
#' @param fromvar Logical (optional), should covariate values be kept in output? Default is `FALSE`.
#' @param fromvar Character (optional), name of variable with starting state. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state. Default is `to`.
#' @param timevar Character (optional), name of variable with time scale. Default is `time`.
#' @param Pvar Character (optional), name of variable with transition probabilities. Default is `P`.
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
                             dtms,
                             constant=NULL,
                             varying=NULL,
                             dropvar=TRUE,
                             timevar="time",
                             fromvar="from",
                             tovar="to",
                             Pvar="P",
                             se=FALSE,
                             vcov=FALSE) {

  # Check
  dtms_proper(dtms)

  # Adjust time scale (transitions in the model)
  timescale_reduced <- dtms$timescale[-length(dtms$timescale)]

  # Get full state space
  all_states <- paste(c(dtms$transient,dtms$absorbing))

  # Create empty frame
  model_frame <- expand.grid(from=dtms$transient,
                             time=timescale_reduced,
                             stringsAsFactors=FALSE)

  # Get names right
  names(model_frame) <- c(fromvar,timevar)

  # Add time-constant covariates
  varnames <- names(constant)
  for(var in varnames) {
    model_frame[var] <- constant[[var]]
  }

  # Add time-varying covariates
  varnames <- names(varying)
  for(var in varnames) {
    # Get values
    value <- varying[[var]]
    # Check if enough values
    if(length(value)!=length(timescale_reduced)) stop("Wrong number of time-varying values")
    # Match to time variable
    assign_values <- match(model_frame[,timevar],timescale_reduced)
    # Assign values
    model_frame[var] <- value[assign_values]
  }

  # Predict
  if(inherits(model,c("vgam","mclogit"))) {
    model_frame[,all_states] <- stats::predict(model,model_frame,"response")[,all_states]
  }

  if(inherits(model,"nnet")) {
    model_frame[,all_states] <- stats::predict(model,model_frame,"probs")[,all_states]
  }

  # Values of starting state
  model_frame[,fromvar] <- paste(model_frame[,fromvar],model_frame[,timevar],sep=dtms$sep)

  # Reshape
  model_frame <- stats::reshape(model_frame,
                                varying=all_states,
                                idvar=fromvar,
                                timevar=tovar,
                                times=all_states,
                                direction="long",
                                v.names=Pvar)

  # SE/CI/vcov
  if(se|vcov) {

    # Simplify starting state (needed for model.matrix below)
    model_frame[,fromvar] <- dtms_simplify(model_frame)$from

    # Coefficients
    if(inherits(model,"mclogit")) {
      C <- stats::coef(model)
      Cstates <- dtms_getstate(names(C),sep="~")
      Cstates <- unique(Cstates)
      C <- matrix(data=C,
                  ncol=length(dtms$transient),
                  byrow=T)
    }

    if(inherits(model,"nnet")) {
      C <- stats::coef(model)
      Cstates <- rownames(C)
      C <- t(C)
    }

    # vcov of coefficients
    Vml <- stats::vcov(model)

    if(inherits(model,"nnet")) Vnames <- dtms_getstate(rownames(Vml),sep=c(":"))
    if(inherits(model,"mclogit")) {
      ordernames <- sort(rownames(Vml))
      Vml <- Vml[ordernames,ordernames]
      Vnames <- dtms_getstate(rownames(Vml),sep=c("~"))
    }

    # Number of probabilities, coefs, states
    nprobs <- dim(model_frame)[1]
    ncoef <- dim(Vml)[1]
    nstates <- length(all_states)

    # Model matrix
    form <- stats::formula(model)
    mm <- stats::model.matrix(object=form,data=model_frame)

    # Scores
    dscores <- matrix(data=1,
                      ncol=nstates,
                      nrow=nprobs)
    colnames(dscores) <- sort(all_states)
    dscores[,Cstates] <- exp(mm%*%C)

    # Full score (minus one because reference category)
    fullscore <- rowSums(dscores)

    # Parts of full derivative (n'*z-n*z')/z^2
    Z <- matrix(data=fullscore,
                nrow=nprobs,
                ncol=ncoef)

    N <- stats::model.matrix(object=~to,data=model_frame)
    N[,1] <- N[,1]-rowSums(N[,-1])
    N <- rowSums(N*dscores)
    N <- matrix(data=N,
                nrow=nprobs,
                ncol=ncoef)

    varvalues <- do.call("cbind",replicate(nstates-1,mm,simplify=FALSE))
    Zdash <- varvalues*dscores[,Vnames]

    Ndash <- outer(model_frame$to,Vnames,FUN=`==`)
    Ndash <- Ndash*N*varvalues

    # Matrix of derivatives
    G <- (Ndash*Z-N*Zdash)/Z^2

    # Vcov matrix
    Vp <- G%*%Vml%*%t(G)

    # Return vcov matrix
    if(vcov) return(Vp)

    # Full starting state
    model_frame[,fromvar] <- paste(model_frame[,fromvar],model_frame[,timevar],sep=dtms$sep)
    model_frame$se <- sqrt(diag(Vp))

  }

  # Values of receiving state (state name + time)
  rightrows <- model_frame[,tovar]%in%dtms$transient
  oldvalues <- model_frame[rightrows,tovar]
  timevalues <- model_frame[rightrows,timevar]+dtms$timestep
  model_frame[rightrows,tovar] <- paste(oldvalues,
                                        timevalues,
                                        sep=dtms$sep)

  # Drop row names
  rownames(model_frame) <- NULL

  # Drop covariate values for prediction
  if(dropvar) {
    model_frame <- model_frame[,names(model_frame)%in%c(fromvar,tovar,timevar,Pvar,"se")]
  }

  # Class
  class(model_frame) <- c("dtms_probs","data.frame")

  # Return
  return(model_frame)

}
