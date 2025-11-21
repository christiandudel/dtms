#' Nonparametric estimates of transition probabilities
#'
#' @description
#' This function calculates nonparametric estimates of transition probabilities.
#' Standard errors assume that all observations are independent.
#'
#' @details
#' The argument `data` takes a data set in transition format. Predicted
#' transition probabilities are returned as a data frame, and not
#' as a transition matrix. While the latter is required for applying Markov
#' chain methods, the data frame is more convenient for viewing and
#' analyzing the transition probabilities themselves. Standard errors are
#' approximated using binomial standard errors. In case of small cell counts
#' this might be inaccurate.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable with starting state in `data`. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state in `data`. Default is `to`.
#' @param timevar Character (optional), name of variable with time scale in `data`. Default is `time`.
#' @param Pvar Character (optional), name of variable with transition probabilities in the returned data frame. Default is `P`.
#' @param weights Character (optional). Name of variable with survey weights.
#' @param se Logical (optional), return standard errors of predicted probabilities. Default is `TRUE`.
#' @param ci Logical (optional), return confidence intervals? See details. Default is FALSE.
#' @param alpha Numeric (optional), if ci=TRUE, what confidence level is used? Default is 0.05.
#'
#' @returns A data frame with transition probabilities.
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
#' ## Nonparametric transition probabilities
#' probs <- dtms_nonparametric(data=estdata,
#'                             dtms=simple)

dtms_nonparametric <- function(data,
                               dtms,
                               fromvar="from",
                               tovar="to",
                               timevar="time",
                               Pvar="P",
                               weights=NULL,
                               se=TRUE,
                               ci=FALSE,
                               alpha=0.05) {

  # Check
  dtms_proper(dtms)

  # Adjust time scale (transitions in the model)
  timescale_reduced <- dtms$timescale[-length(dtms$timescale)]
  ntime <- length(timescale_reduced)

  # Get full state space
  all_states <- paste(c(dtms$transient,dtms$absorbing))

  # Create empty frame
  model_frame <- expand.grid(from=dtms$transient,
                             to=all_states,
                             time=timescale_reduced,
                             stringsAsFactors=FALSE)

  # Get names right
  names(model_frame) <- c(fromvar,tovar,timevar)

  # Weights per transition
  if(is.null(weights)) data$COUNT <- 1 else
    data <- dtms_rename(data,weights,"COUNT")

  # Warning if missing values
  if(any(is.na(data[,c(fromvar,tovar,timevar)]))) warning("Missing values dropped")

  # Aggregate (denominators)
  formal1 <- paste0("COUNT~",fromvar,"+",timevar)
  formal1 <- stats::as.formula(formal1)
  denominators <- stats::aggregate(formal1,data,FUN=sum,drop=F)

  # Aggregate (numerators)
  formal2 <- paste0("COUNT~",fromvar,"+",tovar,"+",timevar)
  formal2 <- stats::as.formula(formal2)
  numerators <- stats::aggregate(formal2,data,FUN=sum,drop=F)

  # Merge
  probs <- merge(numerators,denominators,by=c(fromvar,timevar))
  model_frame <- merge(model_frame,probs,by=c(fromvar,tovar,timevar))

  # Replace missing with 0
  model_frame$COUNT.x[is.na(model_frame$COUNT.x)] <- 0

  # Calculate
  model_frame[,Pvar] <- model_frame$COUNT.x/model_frame$COUNT.y

  # Standard error/confidence interval/vcov?
  if(se|ci) {

    P <- model_frame$P
    N <- model_frame$COUNT.y
    error <- sqrt( (P*(1-P))/N)

    if(se) model_frame$se <- error
    if(ci) {
      z <- (1-alpha/2)
      z <- stats::qnorm(z)
      model_frame$cilow <- model_frame[,Pvar]-z*error
      model_frame$ciup <- model_frame[,Pvar]+z*error

    }

  }

  # Warning if empty cells etc cause missing values
  if(any(is.na(model_frame[,Pvar]))) warning("Some probabilities are missing")

  # Values of starting state
  model_frame[,fromvar] <- paste(model_frame[,fromvar],model_frame[,timevar],sep=dtms$sep)

  # Values of receiving state (state name + time)
  rightrows <- model_frame[,tovar]%in%dtms$transient
  oldvalues <- model_frame[rightrows,tovar]
  timevalues <- model_frame[rightrows,timevar]+dtms$timestep
  model_frame[rightrows,tovar] <- paste(oldvalues,
                                        timevalues,
                                        sep=dtms$sep)

  # Only keep relevant variables
  model_frame <- model_frame[,names(model_frame)%in%c(fromvar,tovar,timevar,Pvar,"se","CIlow","CIup")]

  # Class
  class(model_frame) <- c("dtms_probs","data.frame")

  # Return
  return(model_frame)

}
