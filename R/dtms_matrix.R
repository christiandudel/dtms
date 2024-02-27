#' Creates a transition matrix from transition probabilities
#'
#' @description
#' This function creates a transiton matrix based on transition probabilities
#' predicted using the function `dtms_transitions`.
#'
#' @param probs Data frame with transition probabilities, as created with `dtms_transitions`.
#' @param dtms DTMS object, as created with `dtms`.
#' @param transient Character (optional), names of transient states
#' @param absorbing Character (optional), names of absorbing states
#' @param timescale Numeric (optional), values of time scale
#' @param fromvar Character, name of variable with sending state. Default is `from`, as is used per default with `dtms_transitions`.
#' @param tovar Character, name of variable with receiving state. Default is `to`, as is used per default with `dtms_transitions`.
#' @param Pvar Character, name of variable with transition probabilities. Default is `P`, as is used per default with `dtms_transitions`.
#' @param enforcedeath Logical, make sure that every unit moves to absorbing state after last value of time scale. Default is TRUE.
#' @param sep Character, separator between state name and value of time scale. Default is `_`.
#' @param rescale Logical, rescale transition probabilities to sum to one? Default is TRUE.
#'
#' @return Returns a transition matrix.
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:20)
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
#' ## Get transition matrix
#' Tp <- dtms_matrix(dtms=simple,
#'                   probs=probs)

dtms_matrix <- function(probs,
                        dtms=NULL,
                        transient=NULL,
                        absorbing=NULL,
                        timescale=NULL,
                        fromvar="from",
                        tovar="to",
                        Pvar="P",
                        enforcedeath=T,
                        sep="_",
                        rescale=T) {

  # Use dtms if provided
  if(!is.null(dtms)) {

    # Check
    dtms_proper(dtms)

    # Use values
    timescale <- dtms$timescale
    absorbing <- dtms$absorbing
    transient <- dtms$transient
  }

  # Combine states and time
  transient_states <- dtms_combine(transient,timescale,sep=sep)
  absorbing <- paste(absorbing)
  all_states <- c(transient_states,absorbing)

  # Get names in probs right
  probs <- dtms_rename(probs,c(fromvar,tovar,Pvar),c("from","to","P"))

  # Subset
  getthem <- probs$from%in%transient_states & probs$to%in%all_states
  probs <- subset(probs,subset=getthem)

  # Total number of transient and absorbing states
  s_states <- length(transient_states)
  a_states <- length(absorbing)
  n_states <- length(all_states)

  # Reshape
  Tm <- stats::reshape(probs[,c("from","to","P")],
                  timevar="to",
                  idvar="from",
                  direction="wide")

  # Edit a bit
  Tm[is.na(Tm)] <- 0
  keepnames <- Tm$from
  Tm <- Tm[,-1]

  # Generate matrix
  Tm <- as.matrix(Tm)
  rownames(Tm) <- keepnames

  # Column names
  oldnames <- strsplit(colnames(Tm),split="[.]")
  oldnames <- lapply(oldnames,function(x) x[2])
  colnames(Tm) <- unlist(oldnames)

  # Add "missing" starting states, if any
  addnames <- rownames(Tm)[!rownames(Tm)%in%colnames(Tm)]
  nadd <- length(addnames)
  if(nadd>0) {
    add <- matrix(data=0,ncol=nadd,nrow=dim(Tm)[1])
    colnames(add) <- addnames
    rownames(add) <- rownames(Tm)
    Tm <- cbind(Tm,add)
  }

  # Add potentially missing final states
  addnames <- colnames(Tm)[!colnames(Tm)%in%rownames(Tm)]
  nadd <- length(addnames)
  if(nadd>0) {
    add <- matrix(data=0,nrow=nadd,ncol=dim(Tm)[2])
    rownames(add) <- addnames
    colnames(add) <- colnames(Tm)
    Tm <- rbind(Tm,add)
  }

  # Add death (the column should already be there)
  Tm <- rbind(Tm,rep(0,n_states))
  rownames(Tm)[(s_states+1):n_states] <- absorbing

  # The dead stay dead (hopefully)
  if(a_states==1) Tm[absorbing,absorbing] <- 1
  if(a_states>1) diag(Tm[absorbing,absorbing]) <- 1

  # Sort a little
  Tm <- Tm[all_states,all_states]

  # Numbers please
  class(Tm) <- "numeric"

  # Make sure everyone dies at the end
  if(enforcedeath==T) {
    last_states <- paste(transient,max(timescale),sep=sep)
    if(length(absorbing)==1) {
      Tm[last_states,] <- 0
      Tm[last_states,absorbing] <- 1
    }
    if(length(absorbing)>1) {
      Tm[last_states,] <- 0
      Tm[last_states,absorbing[1]] <- 1
    }
  }

  # Rescale
  if(rescale) Tm <- t(apply(Tm,1,function(x) x/sum(x)))

  # Class
  class(Tm)[2] <- "dtms_matrix"

  # Return
  return(Tm)

}
