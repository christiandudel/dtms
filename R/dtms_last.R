dtms_last <- function(matrix,# Matrix with transition probabilities generated with dtms_matrix
                      transient=NULL, # Names of transient states
                      time=NULL, # Time scale
                      timestep=NULL, # Step length
                      dtms=NULL,# DTMS model
                      risk, # name of state(s) for which risk is of interest
                      risk_to=NULL, # States to which the states specified in 'risk' are left
                      start_time=NULL,# Starting time, will be lowest time in 'dtms' if NULL
                      start_state=NULL,# Names of starting states, will be all transient states in 'dtms' if NULL
                      start_distr=NULL,# Distribution of starting states for average
                      end_time=NULL, # Time up to which lifetime risks are calculated, making them partial lifetime risks
                      method="mid", # Mid-interval transitions or end of interval transitions, "mid" or "end"
                      sep="_",# Separator
                      total=T, # Should total risk be added?
                      rescale=T){ # Rescale probs to sum to one? If F, totals will be equivalent to lifetime risk


  # Require package Biodem
  require(Biodem)

  # Use dtms if provided
  if(!is.null(dtms) & class(dtms)[2]=="dtms") {
    if(is.null(transient)) transient <- dtms$transient
    if(is.null(time)) {
      time <- dtms$time
      time <- time[-length(time)]
    }
    if(is.null(timestep)) timestep <- dtms$timestep
  }

  # Starting states
  if(is.null(start_state)) {
    starting <- levels(interaction(transient,min(time),sep=sep))
  } else {
    starting <- levels(interaction(start_state,start_time,sep=sep))
  }

  # States of the transition matrix
  states <- rownames(matrix)
  nstates <- length(states)

  # Select subset
  selectorD <- unlist(lapply(strsplit(states,sep),function(z) any(risk%in%z) ))
  if(is.null(risk_to)) selectorU <- !selectorD else selectorU <- unlist(lapply(strsplit(states,sep),function(z) any(risk_to%in%z) ))

  # Get tau
  if(is.null(end_time)) tau <- length(time)-1 else tau <- which(end_time==time)

  # Generate t
  if(is.null(start_time)) t <- 0 else t <- which(start_time==time)-1

  # Partition transition matrix
  P_E <- matrix(data=0,nrow=nstates,ncol=nstates)
  P_S <- matrix(data=0,nrow=nstates,ncol=nstates)
  P_E[selectorD,selectorU] <- matrix[selectorD,selectorU]
  P_S[!selectorD,] <- matrix[!selectorD,]
  P_S[selectorD,!selectorU] <- matrix[selectorD,!selectorU]

  # Special case E_0
  ones <- matrix(data=1,nrow=nstates,ncol=nstates)
  results <- vector("list",1)
  names(results) <- "E_0"
  results[["E_0"]] <- t(mtx.exp(P_S,tau-t+1)%*%ones) * mtx.exp(matrix,t)
  colnames(results[["E_0"]]) <- states

  # Loop for other E_x
  step <- 0.5
  past.steps <- c(0)

  while(!step>=(tau-t+1)) {
    past.steps <- c(past.steps,step)
    e <- step-0.5
    tmp <- vector("list",1)
    names(tmp) <- paste("E",step,sep="_")
    tmp[[paste("E",step,sep="_")]] <- t(mtx.exp(matrix,e)%*%P_E%*%mtx.exp(P_S,tau-t-e)%*%ones) * mtx.exp(matrix,t)
    colnames(tmp[[paste("E",step,sep="_")]]) <- states
    results <- c(results,tmp)
    step <- step+1
  }

  # Get distribution
  tmp <- unlist(lapply(results, function(z) colSums(z)[starting]))

  # Result object
  result <- matrix(data=tmp,ncol=length(past.steps),nrow=length(starting))
  rownames(result) <- paste0("start:",starting)

  # Get times right (interval length, mid interval vs end of interval)
  steps <- past.steps*timestep
  if(method=="end") steps[-1] <- steps[-1]+0.5*timestep
  colnames(result) <- paste(steps)

  # Average
  if(!is.null(start_distr)) {

    # Overall average
    AVERAGE <- start_distr%*%result
    rownames(AVERAGE) <- "AVERAGE"
    result <- rbind(result,AVERAGE)

    # Conditional on not starting in state in risk set
    whererisk <- !transient%in%risk
    tmp_distr <- start_distr[whererisk]/sum(start_distr[whererisk])
    tmp_names <- paste0("start:",names(tmp_distr))
    tmp <- tmp_distr%*%result[tmp_names,]
    result <- rbind(result,tmp)
    rownames(result)[dim(result)[1]] <- "AVERAGE(COND.)"

  }

  # Rescale
  if(rescale) {
    result <- result[,-1]
    result <- t(apply(result,1,function(x) x/sum(x)))
  }

  # Total
  if(total) {
    TOTAL <- rowSums(result)
    result <- cbind(result,TOTAL)
    if(rescale) colnames(result)[dim(result)[2]] <- "TOTAL(RESCALED)"
  }

  # Output
  return(result)

} # End of function
