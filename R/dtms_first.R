dtms_first  <- function(matrix,# Matrix with transition probabilities generated with dtms_matrix
                        transient=NULL, # Names of transient states
                        time=NULL, # Time scale
                        timestep=NULL, # Step length of the model
                        dtms=NULL,# DTMS model
                        risk, # name of state(s) for which risk is of interest
                        start_time=NULL,# Starting time, will be lowest time in 'dtms' if NULL
                        start_state=NULL,# Names of starting states, will be all transient states in 'dtms' if NULL
                        start_distr=NULL,# Distribution of starting states for average
                        end_time=NULL, # Time up to which lifetime risks are calculated, making them partial lifetime risks
                        method="mid", # Mid-interval transitions or end of interval transitions, "mid" or "end"
                        sep="_",# Separator
                        total=T, # Should total risk be added?
                        rescale=T, # Rescale probs to sum to one? If F, totals will be equivalent to lifetime risk
                        maxiter=200) { # Maximum iterations, safeguard to avoid non-ending loop

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
  selectorU <- unlist(lapply(strsplit(states,sep),function(y) any(risk%in%y) ))
  P_E <- matrix(data=0,nrow=nstates,ncol=nstates)
  P_E[!selectorU,!selectorU] <- matrix[!selectorU,!selectorU]

  # Generate t
  if(is.null(start_time)) t <- 0 else t <- which(start_time==time)-1

  # Generate max
  if(is.null(end_time)) maxtime <- length(time)-1
  if(!is.null(end_time)) if(!is.infinite(end_time))  {
    maxtime <- which(end_time==time)-1
    maxiter <- max(maxtime,maxiter)
  } else maxtime <- Inf

  # Generate W_t_0 and W_t_0.5, initial conditions
  results <- vector("list",2)
  names(results) <- c("W_0","W_0.5")
  results[["W_0"]] <- mtx.exp(matrix,t)
  if(t==0) {
    results[["W_0.5"]] <- matrix(data=0,ncol=nstates,nrow=nstates)
    diag(results[["W_0.5"]][!selectorU,!selectorU]) <- 1
  } else results[["W_0.5"]] <- mtx.exp(P_E,t)

  # Variables
  upcoming <- 1.5
  past.steps <- c(0,0.5)
  steps <- 1
  end <- 0

  # Loop to generate results
  while(end<5) {
    past.steps <- c(past.steps,upcoming)
    tmp <- vector("list",1)
    tmp[[1]] <- t(mtx.exp(P_E,upcoming-0.5)%*% matrix(data=1,nrow=nstates,ncol=nstates))*results[["W_0.5"]]
    names(tmp) <- paste("W",upcoming,sep="_")
    results <- c(results,tmp)
    steps <- steps+1
    if(steps>maxtime) end <- 5
    if(is.infinite(maxtime) & identical(results[[paste("W",upcoming-1,sep="_")]], results[[paste("W",upcoming,sep="_")]])) end <- end+1
    if(steps>maxiter) end <- 5
    upcoming <- upcoming+1
  }

  # Generate V
  results_V <- vector("list",length(past.steps)-1)
  for(i in 1:length(results_V)) {
    names(results_V)[i] <- paste("V",past.steps[i],sep="_")
    results_V[[paste("V",past.steps[i],sep="_")]] <- results[[paste("W",past.steps[i],sep="_")]]-results[[paste("W",past.steps[i+1],sep="_")]]
    colnames(results_V[[paste("V",past.steps[i],sep="_")]]) <- states
  }

  # Distribution
  tmp <- unlist(lapply(results_V, function(y) colSums(y)[starting]))

  # Result object
  result <- matrix(data=tmp,ncol=length(past.steps)-1,nrow=length(starting))
  rownames(result) <- paste0("start:",starting)

  # Get times right (interval length, mid interval vs end of interval)
  steps <- past.steps*timestep
  if(method=="end") steps[-1] <- steps[-1]+0.5*timestep
  colnames(result) <- paste(steps)[-length(steps)]

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
