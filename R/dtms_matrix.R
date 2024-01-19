dtms_matrix <- function(transient=NULL, # character vector of transient states
                        absorbing=NULL, # character vector of absorbing states
                        time=NULL, # Time steps
                        dtms=NULL, # DTMS object
                        probs, # probs frame with probabilities
                        fromvar="from", # Variable name of starting state in 'probs'
                        tovar="to", # Variable name of receiving state in 'probs'
                        Pvar="P", # Variable name of transition probability in 'probs'
                        enforcedeath=T, # Make sure that at t=T everyone dies?
                        sep="_", # State/time separator
                        rescale=T) { # Rescale transition probabilities to sum to 1

  # Require tidy
  require(dplyr)
  require(magrittr)

  # Use dtms if provided
  if(!is.null(dtms) & class(dtms)[2]=="dtms") {
    transient <- dtms$transient
    absorbing <- dtms$absorbing
    time <- dtms$time
    time <- time[-length(time)]
  }

  # Combine states and time
  transient_states <- levels(interaction(transient,time,sep=sep))
  absorbing <- paste(absorbing)
  all_states <- c(transient_states,absorbing)

  # Get names in probs right
  probs <- probs %>% rename_with(~ c("from","to","P"),c(fromvar,tovar,Pvar))

  # Subset
  probs <- probs %>% filter(from%in%transient_states & to%in%all_states)

  # Total number of transient and absorbing states
  s_states <- length(transient_states)
  a_states <- length(absorbing)
  n_states <- length(all_states)

  # Reshape
  Tm <- probs %>%
    select(from,to,P) %>%
    pivot_wider(id_cols     = c(from),
                names_from  = to,
                values_from = P)

  # Edit a bit
  Tm[is.na(Tm)] <- 0
  keepnames <- Tm$from
  Tm <- Tm[,-1]

  # Generate matrix
  Tm <- as.matrix(Tm)
  rownames(Tm) <- keepnames

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
  rownames(Tm)[(s_states+1):n_states] <- absorbing # HIER WEITER => Character string?

  # The dead stay dead (hopefully)
  if(a_states==1) Tm[absorbing,absorbing] <- 1
  if(a_states>1) diag(Tm[absorbing,absorbing]) <- 1

  # Sort a little
  Tm <- Tm[all_states,all_states]

  # Numbers please
  class(Tm) <- "numeric"

  # Make sure everyone dies at the end
  if(enforcedeath==T) {
    last_states <- paste(transient,max(time),sep=sep)
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

  # Return
  return(Tm)

}
