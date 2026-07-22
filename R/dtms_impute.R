#' Impute missing transition probabilities
#'
#' @description
#' This function imputes missing transition probabilities.
#'
#' @details
#' Missing transition probabilities might occur after using
#' \code{dtms_nonparametic} if for a given combination of starting state and
#' value of the time scale there are no observations. For instance, for some
#' state `A` at time 10, there might be no units in this state which could
#' transition to another state.
#'
#' There are three options to impute missing transition probabilities which
#' are selected with the argument `method`: `even`, `absorbing`, and `previous`.
#' `even` assume that all states, both transient and absorbing, have the same
#' probability of being reached. For instance, if there are two transient states
#' `A` and `B` and an absorbing state `X`, then it is assumed that the
#' transition probabilities to these states all equal 1/3. `absorbing`
#' assumes that the probability of being absorbed equals 1. If there are
#' several absorbing states, then the first absorbing state as specified in
#' the dtms object through the argument `dtms` is used. `previous` takes the
#' transition probabilities from the same starting state at the previous value
#' of the time scale. For instance, if transition probabilities are missing
#' for state `A` at time 10, then the transition probabilities from state `A` at
#' time 9 are used. This is applied recursively; e.g., if transition
#' probabilities are missing at times 8, 9, and 11, then first transition
#' probabilities are imputed for time 8 using those from time 7; next, missing
#' probabilities are imputed for time 9 using those from time 8 and thus from
#' time 8; the missing probabilities for time 11 are taken from time 10.
#'
#' The argument `reshapesep` is used for reshaping in an internal call and
#' passed to the argument `sep` of the base R function \code{reshape}.
#' It does not affect the computation.
#' Depending on the names of the variables and states, however, it might be
#' necessary to choose something different than the default.
#'
#' @param probs Data frame with transition probabilities, as created with \code{dtms_transitions} or \code{dtms_nonparametric}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable with starting state in `data`. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state in `data`. Default is `to`.
#' @param timevar Character (optional), name of variable with time scale in `data`. Default is `time`.
#' @param Pvar Character (optional), name of variable with transition probabilities in the returned data frame. Default is `P`.
#' @param method Character (optional), method for imputation, see details. Default is `even` .
#' @param reshapesep Character (optional), separator for internal reshaping, see details. Default is `.`.
#'
#' @returns A data frame with transition probabilities.
#' @export
#'
#' @examples
#' #' ## Define model: Absorbing and transient states, time scale
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
#'                             dtms=simple,
#'                             se=FALSE)
#' ## Set some transitions to NA
#' probs[probs$from=="A_9","P"] <- NA
#' ## Impute
#' dtms_impute(probs=probs,
#'             dtms=simple)

dtms_impute <- function(probs,
                        dtms,
                        fromvar="from",
                        tovar="to",
                        timevar="time",
                        Pvar="P",
                        method="even",
                        reshapesep=".") {

  # Check
  dtms_proper(dtms)

  # Reshape data: starting state to rows
  tmp <- dtms_simplify(probs[,c(fromvar,tovar,timevar,Pvar)])
  tmp <- stats::reshape(tmp,
                        direction="wide",
                        idvar=c(fromvar,timevar),
                        timevar=tovar,
                        v.names=Pvar,
                        sep=reshapesep)

  # Impute: even
  if(method=="even") {
    # Which rows/starting states need fixing
    whichrows <- apply(tmp[,-c(1:2)],1,function(x) any(is.na(x)))
    # Number of columns to fill
    manycols <- dim(tmp)[2]-2
    # Values
    fillvec <- rep(1/manycols,manycols)
    # Place in rows
    tmp[whichrows,-c(1:2)] <- fillvec
  }

  # Impute: absorbing
  if(method=="absorbing") {
    # Which rows/starting states need fixing
    whichrows <- apply(tmp[,-c(1:2)],1,function(x) any(is.na(x)))
    # Number of columns to fill
    manycols <- dim(tmp)[2]-2
    # Mostly zeros
    fillvec <- rep(0,manycols)
    # Where in vector is (first) absorbing state
    firstabsorb <- paste0(Pvar,".",dtms$absorbing[1])
    whichabsorb <- which(names(tmp[,-c(1:2)])==firstabsorb)
    # Absorbing state gets 1
    fillvec[whichabsorb] <- 1
    # Place in rows
    tmp[whichrows,-c(1:2)] <- fillvec
  }

  # Impute: previous
  if(method=="previous") {

    # Sort data
    tmp <- tmp[order(tmp[,fromvar], tmp[,timevar]),]

    # Which rows/starting states need fixing
    whichrows <- apply(tmp[,-c(1:2)],1,function(x) any(is.na(x)))
    rownumbers <- (1:dim(tmp)[1])[whichrows]

    for(number in rownumbers) {

      # Get current state and time
      currentstate <- tmp[number,fromvar]
      currenttime <- tmp[number,timevar]

      # Get previous time
      previoustime <- which(dtms$timescale==currenttime)
      previoustime <- previoustime - 1
      if(previoustime<=0) next
      previoustime <- dtms$timescale[previoustime]
      whichrow <- which(tmp[,fromvar]==currentstate & tmp[,timevar]==previoustime)

      # Get values
      fillvec <- tmp[whichrow,-c(1,2)]

      # Place in rows
      tmp[whichrows,-c(1:2)] <- fillvec

    }

  }

  # Reshape back to long
  tmp <- stats::reshape(tmp,
                        direction="long",
                        varying=names(tmp)[-c(1:2)],
                        sep=reshapesep,
                        timevar=tovar)

  # Get variables right
  tmp <- tmp[,-dim(tmp)[2]]
  tmp[,fromvar] <- paste0(tmp[,fromvar],dtms$sep,tmp[,timevar])
  tmp[,tovar][tmp[,tovar]%in%dtms$transient] <- paste0(tmp[,tovar][tmp[,tovar]%in%dtms$transient],
                                                       dtms$sep,
                                                       tmp[,timevar][tmp[,tovar]%in%dtms$transient]+1)

  # Remove P var in probs
  dropvar <- which(names(probs)==Pvar)
  probs <- probs[,-dropvar]

  # Replace existing probabilities
  probs <- merge(probs,tmp,
                 by=c(fromvar,tovar,timevar),
                 sort=FALSE)

  # Return
  return(probs)

}
