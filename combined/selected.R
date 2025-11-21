dtms <- function(transient,
                 absorbing,
                 timescale,
                 timestep=NULL,
                 sep="_") {
  # Guess time step?
  if(is.null(timestep)) {
    # Step lengths as specified by time scale
    step <- timescale |> diff() |> unique()
    # Number of different step lengths
    nstep <- length(step)
    # Guess timestep if possible, else error
    if(nstep==1) timestep <- step else
      stop("Not able to guess step length of time scale")
  }
  # Combine everything in a list
  result <- list(transient=transient,
                 absorbing=absorbing,
                 timescale=timescale,
                 timestep=timestep,
                 sep=sep)
  # Assign class
  class(result)[2] <- "dtms"
  # Return
  return(result)
}
dtms_format <- function(data,
                        dtms,
                        idvar="id",
                        timevar="time",
                        statevar="state",
                        fromvar="from",
                        tovar="to",
                        absorbing=TRUE,
                        keepnames=FALSE,
                        fill=FALSE,
                        verbose=TRUE,
                        steplength=FALSE,
                        stepvar="steplength") {
  # Transform to data frame, e.g., if tibble
  if(class(data)[1]!="data.frame") data <- as.data.frame(data)
  # Check
  dtms_proper(dtms)
  # Fill data
  if(fill) {
    # Get ID values
    idvalues <- data[,idvar] |> unique()
    # Full data
    fulldata <- expand.grid(dtms$timescale,idvalues,
                            stringsAsFactors=FALSE)
    names(fulldata) <- c(timevar,idvar)
    # Merge with data
    data <- merge(fulldata,data,
                  by=c(idvar,timevar),
                  all=T)
    # Drop temporary data
    rm(fulldata)
  }
  # Sort data
  dataorder <- order(data[,idvar],
                     data[,timevar])
  data <- data[dataorder,]
  # Check if next value is valid
  consecutive <- dtms_consecutive(data=data,
                                  idvar=idvar,
                                  timevar=timevar,
                                  timestep=dtms$timestep)
  # Absorbing states carry forward
  if(absorbing) {
    tmp <- tapply(data[,statevar],data[,idvar],function(x) dtms_carry(x=x,dtms=dtms))
    data[,statevar] <- unlist(tmp)
  }
  # Get next state
  data[,tovar] <- NA
  tovalues <- c(data[-1,statevar],NA)
  data[consecutive$true,tovar] <- tovalues[consecutive$true]
  if(steplength) data[consecutive$true,stepvar] <- consecutive$numeric[consecutive$true]
  # Rename from variable
  data <- dtms_rename(data,statevar,fromvar)
  # Change names of id and time if possible
  if(!keepnames) {
    if(timevar!='time' & !'time'%in%names(data))
      data <- dtms_rename(data,timevar,"time") else
      if(verbose) cat("Kept original name for time \n")
    if(idvar!='id' & !'id'%in%names(data))
      data <- dtms_rename(data,idvar,"id") else
      if(verbose) cat("Kept original name for id \n")
  }
  # Class
  class(data) <- c("dtms_data","data.frame")
  # Return result
  return(data)
}
dtms_clean <- function(data,
                       dtms,
                       fromvar="from",
                       tovar="to",
                       timevar="time",
                       dropTime=TRUE,
                       dropState=TRUE,
                       dropNA=TRUE,
                       dropAbs=TRUE,
                       verbose=TRUE) {
  # Check
  dtms_proper(dtms)
  # Drop observations not in state space
  if(dropState) {
    allstates <- c(dtms$transient,dtms$absorbing,NA)
    whichrows <- unlist(data[,fromvar])%in%allstates & unlist(data[,tovar])%in%allstates
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows not in state space\n")
    }
  }
  # Drop observations not in time range
  if(dropTime) {
    whichrows <- unlist(data[,timevar])%in%dtms$timescale
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows not in time range\n")
    }
  }
  # Drop missing values
  if(dropNA) {
    whichrows <- !is.na(unlist(data[,fromvar])) & !is.na(unlist(data[,tovar]))
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows starting or ending in NA\n")
    }
  }
  # Drop transitions starting in absorbing states
  if(dropAbs) {
    whichrows <- !unlist(data[,fromvar])%in%dtms$absorbing
    data <- data[whichrows,]
    if(verbose) {
      count <- sum(!whichrows)
      cat("Dropping ",count," rows starting in absorbing state\n")
    }
  }
  # Return
  return(data)
}
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
dtms_aggregate <- function(data,
                           weights=NULL,
                           idvar="id",
                           countvar="count") {
  # Transform to data frame, e.g., if tibble
  if(class(data)[1]!="data.frame") data <- as.data.frame(data)
  # Get variable names without weights, collapse for formula
  controls <- names(data)
  controls <- controls[which(controls!=idvar)]
  if(!is.null(weights)) controls <- controls[which(controls!=weights)]
  controls <- paste(controls,collapse="+")
  # If no weights
  if(is.null(weights)) {
    data[,countvar] <- 1
    aggformula <- paste0(countvar,"~",controls)
  # If weights
  } else {
    aggformula <- paste0(weights,"~",controls)
  }
  # Formula for aggregate
  aggformula <- stats::as.formula(aggformula)
  # Warning if missing values
  drops <- data |>
    stats::na.omit() |> dim()
  drops <- dim(data)[1]-drops[1]
  if(drops>0) warning(paste("Dropping",drops,"rows with missing values"))
  # Aggregate
  tmp <- stats::aggregate(by=aggformula,
                          x=data,
                          FUN=sum)
  # Return
  return(tmp)
}
dtms_fit <- function(data,
                     controls=NULL,
                     formula=NULL,
                     weights=NULL,
                     fromvar="from",
                     tovar="to",
                     reference=1,
                     package="VGAM",
                     full=FALSE,
                     ...) {
  # Require package used for estimation (requireNamespace does not help here)
  require(package,character.only=TRUE,quietly=TRUE)
  # Build formula if not specified
  if(is.null(formula)) formula <- dtms_formula(controls=controls,
                                               fromvar=fromvar,
                                               tovar=tovar,
                                               full=full)
  # Make sure environment for formula is correct (ugh)
  environment(formula) <- environment()
  # Get weights if specified
  if(!is.null(weights)) weights <- data[,weights]
  # Factors (needed by most packages)
  data[,fromvar] <- as.factor(data[,fromvar])
  data[,tovar] <- as.factor(data[,tovar])
  data[,tovar] <- stats::relevel(data[,tovar],ref=reference)
  # VGAM
  if(package=="VGAM") {
    # Estimate
    fitted <- VGAM::vgam(formula=formula,
                         family=VGAM::multinomial(refLevel=reference),
                         data=data,
                         weights=weights,
                         ...)
  }
  #nnet
  if(package=="nnet") {
    # Estimate
    fitted <- nnet::multinom(formula=formula,
                             data=data,
                             weights=weights,
                             ...)
  }
  #mclogit
  if(package=="mclogit") {
    # Estimate
    fitted <- mclogit::mblogit(formula=formula,
                               data=data,
                               weights=weights,
                               ...)
  }
  # Return results
  return(fitted)
}
dtms_transitions <- function(model,
                             dtms,
                             controls=NULL,
                             dropvar=TRUE,
                             timevar="time",
                             fromvar="from",
                             tovar="to",
                             Pvar="P",
                             se=TRUE,
                             vcov=FALSE,
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
                             time=timescale_reduced,
                             stringsAsFactors=FALSE)
  # Get names right
  names(model_frame) <- c(fromvar,timevar)
  # Deal with controls
  varnames <- names(controls)
  for(var in varnames) {
    # Get values
    value <- controls[[var]]
    nvalue <- length(value)
    # Check
    if(!nvalue%in%c(1,ntime,ntime+1)) stop("Wrong number of time-varying values")
    # Act depending on how many values
    if(nvalue==1) model_frame[var] <- value
    if(nvalue==ntime) {
      assign_values <- match(model_frame[,timevar],timescale_reduced)
      model_frame[var] <- value[assign_values]
    }
    if(nvalue==ntime+1) {
      value <- value[-nvalue]
      assign_values <- match(model_frame[,timevar],timescale_reduced)
      model_frame[var] <- value[assign_values]
    }
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
  if(se|vcov|ci) {
    # Simplify starting state (needed for model.matrix below)
    model_frame[,fromvar] <- dtms_simplify(model_frame)$from
    # Coefficients
    if(inherits(model,"mclogit")) {
      C <- stats::coef(model)
      Cstates <- dtms_getstate(names(C),sep="~")
      Cstates <- unique(Cstates)
      C <- matrix(data=C,
                  ncol=length(all_states)-1,
                  byrow=T)
    }
    if(inherits(model,"vgam")) {
      C <- stats::coef(model)
      C <- matrix(C,
                  ncol=length(all_states)-1,
                  byrow=T)
      Cstates <- model@extra$colnames.y[-model@extra$use.refLevel]
    }
    if(inherits(model,"nnet")) {
      C <- stats::coef(model)
      Cstates <- rownames(C)
      C <- t(C)
    }
    # vcov of coefficients
    Vml <- stats::vcov(model)
    if(inherits(model,"mclogit")) {
      ordernames <- sort(rownames(Vml))
      Vml <- Vml[ordernames,ordernames]
      Vnames <- dtms_getstate(rownames(Vml),sep=c("~"))
    }
    if(inherits(model,"vgam")) {
      # Get nice names (assigned below)
      nicenames <- dtms_getstate(rownames(Vml),sep=c(":"))
      nicenames <- unique(nicenames)
      nicenames <- sort(dtms_combine(Cstates,nicenames,sep=":"))
      # Reorder
      newnames <- unlist(lapply(strsplit(colnames(Vml),split=":"),function(x) paste0(x[2],":",x[1])))
      colnames(Vml) <- rownames(Vml) <- newnames
      ordernames <- sort(rownames(Vml))
      Vml <- Vml[ordernames,ordernames]
      # Assign nice names
      colnames(Vml) <- rownames(Vml) <- nicenames
      # State names
      Vnames <- dtms_getstate(rownames(Vml),sep=c(":"))
    }
    if(inherits(model,"nnet")) Vnames <- dtms_getstate(rownames(Vml),sep=c(":"))
    # Number of probabilities, coefs, states
    nprobs <- dim(model_frame)[1]
    ncoef <- dim(Vml)[1]
    nstates <- length(all_states)
    # Model matrix
    form <- stats::formula(model)
    mm <- stats::model.matrix(object=form,data=model_frame)
    # Scores (denominator for predicted prob)
    dscores <- matrix(data=1,
                      ncol=nstates,
                      nrow=nprobs)
    colnames(dscores) <- sort(all_states)
    dscores[,Cstates] <- exp(mm%*%C)
    # Full score (numerator for predicted prob)
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
    G <- (Ndash*Z-N*Zdash)/(Z^2)
    # Re-order
    if(inherits(model,"mclogit")) colnames(G) <- dtms_combine(Cstates,colnames(mm),sep="~")
    if(inherits(model,"vgam")) colnames(G) <- dtms_combine(Cstates,colnames(mm),sep=":")
    if(inherits(model,"nnet")) colnames(G) <- dtms_combine(Cstates,colnames(mm),sep=":")
    G <- G[,rownames(Vml)]
    # Vcov matrix
    Vp <- G%*%Vml%*%t(G)
    # Return vcov matrix
    if(vcov) return(Vp)
    # Full starting state
    model_frame[,fromvar] <- paste(model_frame[,fromvar],model_frame[,timevar],sep=dtms$sep)
    # SE
    if(se) model_frame$se <- sqrt(diag(Vp))
    # CI
    if(ci) {
      z <- (1-alpha/2)
      z <- stats::qnorm(z)
      error <- sqrt(diag(Vp))
      model_frame$cilow <- model_frame[,Pvar]-z*error
      model_frame$ciup <- model_frame[,Pvar]+z*error
    }
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
    model_frame <- model_frame[,names(model_frame)%in%c(fromvar,tovar,timevar,Pvar,"se","cilow","ciup")]
  }
  # Class
  class(model_frame) <- c("dtms_probs","data.frame")
  # Return
  return(model_frame)
}
dtms_matrix <- function(probs,
                        dtms=NULL,
                        fromvar="from",
                        tovar="to",
                        Pvar="P",
                        enforcedeath=T,
                        rescale=T,
                        reshapesep=":") {
  # Check
  dtms_proper(dtms)
  # Combine states and time
  transient_states <- dtms_combine(dtms$transient,dtms$timescale,sep=dtms$sep)
  absorbing <- paste(dtms$absorbing)
  all_states <- c(transient_states,absorbing)
  # Get variable names in probs right
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
                  direction="wide",
                  sep=reshapesep)
  # Edit a bit
  Tm[is.na(Tm)] <- 0
  keepnames <- Tm$from
  Tm <- Tm[,-1]
  # Generate matrix
  Tm <- as.matrix(Tm)
  rownames(Tm) <- keepnames
  # Column names
  oldnames <- strsplit(colnames(Tm),split=paste0("[",reshapesep,"]"))
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
    last_states <- paste(dtms$transient,max(dtms$timescale),sep=dtms$sep)
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
  class(Tm) <- c("dtms_matrix","matrix")
  # Return
  return(Tm)
}
dtms_expectancy <- function(matrix,
                            dtms,
                            risk=NULL,
                            start_distr=NULL,
                            start_time=NULL,
                            start_state=NULL,
                            end_time=NULL,
                            correction=0.5,
                            total=TRUE,
                            fundamental=FALSE,
                            verbose=FALSE) {
  # Check
  dtms_proper(dtms)
  # Starting state and time
  if(is.null(start_state)) start_state <- dtms$transient
  if(is.null(start_time)) start_time <- min(dtms$timescale)
  # Starting states, long names
  starting <- dtms_combine(start_state,start_time,sep=dtms$sep)
  # Number of starting and receiving states
  nstart <- length(starting)
  ntransient <- length(dtms$transient)
  # Remove absorbing states
  matrix <- dtms_absorbing(matrix)
  # All transient states
  allstates <- rownames(matrix)
  # Fundamental matrix
  nstates <- dim(matrix)[1]
  Nmat <- solve(diag(1,nstates)-matrix)
  # Correction
  if(is.numeric(correction)) {
    # Adjust
    diag(Nmat) <- diag(Nmat) - correction
    # Output
    if(verbose) cat("(Applying correction)","\n\n")
  }
  # Only return fundamental matrix?
  if(fundamental) {
    return(Nmat)
  }
  # Variant 1: Expectation of all transient states
  if(is.null(risk)) {
    # Matrix for results
    result <- matrix(data=NA,ncol=ntransient,nrow=nstart)
    rownames(result) <- paste0("start:",starting)
    colnames(result) <- dtms$transient
    for(i in 1:ntransient) {
      # Get states
      selector <- dtms_in(allstates,dtms$transient[i],dtms$sep)
      # Use end_time if specified
      if(!is.null(end_time)) {
        times <- dtms_gettime(allstates,dtms$sep)
        times <- times<=end_time
        times[!is.logical(times)] <- F
        selector <- selector & times
      }
      # Calculate results and place
      if(nstart>1) tmp <- rowSums(Nmat[starting,selector]) else tmp <- sum(Nmat[starting,selector])
      # Place
      result[,dtms$transient[i]] <- tmp
    }
  }
  # Variant 2: Expectation in one state by time scale
  if(!is.null(risk)) {
    # Check
    if(length(risk)!=1) stop("Only one state allowed for 'risk'")
    # Get time right
    first <- which(dtms$timescale==start_time)
    if(is.null(end_time)) last <- length(dtms$timescale) else
      last <- which(dtms$timescale==end_time)
    times <- dtms$timescale[first:last]
    ntimes <- length(times)
    # Get right columns from fundamental matrix
    selector1 <- dtms_in(allstates,risk,dtms$ep)
    selector2 <- dtms_gettime(allstates,dtms$sep)%in%times
    selector <- selector1 & selector2
    # Get result
    tmp <- rowSums(Nmat[,selector])
    # Matrix with results
    result <- matrix(data=tmp,ncol=ntimes,nrow=nstart,byrow=T)
    rownames(result) <- paste0("start:",starting)
    colnames(result) <- paste(times)
  }
  # Calculate average if starting distribution is provided
  if(!is.null(start_distr) & is.null(risk)) {
    # Check if matching
    if(length(start_distr)!=dim(result)[1]) stop("Starting distribution too long or short")
    # Match to starting/row ordering of result
    start_distr <- start_distr[match(names(start_distr),starting)]
    # Calculate
    AVERAGE <- colSums(result*start_distr)
    # Put into matrix for results
    result <- rbind(result,AVERAGE)
  }
  # Add row totals
  if(total & is.null(risk)) {
    TOTAL <- rowSums(result)
    result <- cbind(result,TOTAL)
  }
  # Adjust for time step
  if(dtms$timestep!=1) {
    result <- result*dtms$timestep
    if(verbose) cat("Adjusting for step length","\n\n")
  }
  # Return result
  return(result)
}
dtms_absorbing <- function(matrix) { # matrix=full transition matrix
  # Get states which are absorbing
  to_remove <- which(diag(matrix)==1)
  # Remove from matrix
  removed <- matrix[-to_remove,-to_remove]
  # Output reduced matrix
  return(removed)
}
dtms_proper <- function(dtms) { # dtms=object to be checked
  # Error message
  message <- "Not a proper dtms object."
  # Check class
  if(!class(dtms)[2]=="dtms") stop(message)
  # Check names
  listnames <- c("transient" ,"absorbing" ,"timescale" ,"timestep" ,"sep")
  if(!all(names(dtms)==listnames)) stop(message)
  # Check types (and length)
  if(!is.character(dtms$transient)&!is.numeric(dtms$transient)) stop(message)
  if(!is.character(dtms$absorbing)&!is.numeric(dtms$absorbing)) stop(message)
  if(!is.character(dtms$sep)) stop(message)
  if(!is.numeric(dtms$timescale) | length(dtms$timescale)<2) stop(message)
  if(!is.numeric(dtms$timestep)) stop(message)
}
dtms_consecutive <- function(data,idvar,timevar,timestep) {
  # Make sure no missing times and ids
  if(any(is.na(data[,timevar]))) stop("Missing values in time variable not allowed")
  if(any(is.na(data[,idvar]))) stop("Missing values in ID variable not allowed")
  # Get diff to next time step
  consecutive <- by(data[,timevar],data[,idvar],FUN=diff)
  # Add last obs
  consecutive <- lapply(consecutive,function(x) c(x,-1))
  # Unlist
  consecutive <- unlist(consecutive)
  # TRUE if equal to timestep, FALSE otherwise
  result <- data.frame(true=consecutive%in%timestep,
                       numeric=consecutive)
  # Return
  return(result)
}
dtms_rename <- function(data,oldnames,newnames) {
  # Which names to change?
  changenames <- match(oldnames,names(data))
  # Change names
  names(data)[changenames] <- newnames
  # Return
  return(data)
}
dtms_combine <- function(values1,values2,sep) {
  # Generate output vector
  output <- character(0)
  # Get values
  for(value in values1) {
    output <- c(output,paste(value,values2,sep=sep))
  }
  # Return
  return(output)
}
dtms_in <- function(vector,name,sep) {
  res <- lapply(strsplit(vector,split=sep),function(y) any(name%in%y[[1]]) )
  res <- unlist(res)
  return(res)
}
dtms_gettime <- function(vector,sep) {
  res <- lapply(strsplit(vector,split=sep),function(y) y[2] )
  res <- unlist(res) |> as.numeric()
  return(res)
}
dtms_getstate <- function(vector,sep) {
  res <-lapply(strsplit(vector,split=sep),function(x) x[1])
  res <- unlist(res)
  return(res)
}
dtms_mtexp <- function(matrix,n) {
  res <- diag(nrow=nrow(matrix))
  rep <- 0
  while(rep<n) {
    res <- res %*% matrix
    rep <- rep+1
  }
  return(res)
}
dtms_formula <- function(controls, # Arguments the same as for dtms_fit
                         fromvar,
                         tovar,
                         full) {
  # If fromvar is NULL (not default, needs explicit call)
  if(is.null(fromvar)) fromvar <- "1"
  # Constrained model
  if(!full) {
    formula <- paste0(tovar,"~",fromvar)
    if(!is.null(controls)) {
      controls <- paste(controls,collapse="+")
      formula <- paste(formula,controls,sep="+")
    }
    formula <- stats::as.formula(formula)
    # Unconstrained/fully interacted model
  } else {
    formula <- paste0(tovar,"~",fromvar)
    if(!is.null(controls)) {
      varlist <- paste(controls,fromvar,sep="*")
      varlist <- paste(varlist,collapse="+")
      formula <- paste(formula,varlist,sep="+")
    }
    formula <- stats::as.formula(formula)
  }
  # Return
  return(formula)
}
dtms_lag <- function(data,
                     dtms,
                     lag,
                     fromvar="from",
                     idvar="id",
                     timevar="time") {
  # Make data smaller
  data <- data[,c(idvar,timevar,fromvar)]
  # Get ID values
  idvalues <- data[,idvar] |> unique()
  # Full data
  fulldata <- expand.grid(dtms$timescale,idvalues,
                          stringsAsFactors=FALSE)
  names(fulldata) <- c(timevar,idvar)
  # Merge with data
  fulldata <- merge(fulldata,data,
                    by=c(idvar,timevar),
                    all=T)
  # shift state
  stateshift <- by(fulldata[,fromvar],
                   fulldata[,idvar],function(x) c(rep(NA,min(lag,length(x))),
                                                  x[-( max(0,(length(x)-(lag-1))): length(x))]))
  # shift time
  timeshift <- by(fulldata[,timevar],
                  fulldata[,idvar],function(x) c(rep(NA,min(lag,length(x))),
                                                 diff(x,lag=lag) ))
  # Unlist
  stateshift <- unlist(stateshift)
  timeshift <- unlist(timeshift)
  # drop wrong spacing
  stateshift[!timeshift%in%c(NA,dtms$timestep*lag)] <- NA
  # Merge back
  fulldata$stateshift <- stateshift
  data <- merge(data,fulldata,by=c(idvar,timevar,fromvar),sort=FALSE)
  # Return
  return(data$stateshift)
}
dtms_carry <- function(x,
                       dtms) {
  # Check if action necessary
  if(any(x%in%dtms$absorbing)) {
    # Length
    n <- length(x)
    # From where to carry over
    whichfirst <- which(x%in%dtms$absorbing)[1]
    # What to carry over
    whichvalue <- x[whichfirst]
    # Carry over
    x[whichfirst:n] <- whichvalue
    # Return
    return(x)
  } else return(x)
}
dtms_forward_help <- function(x, # Vector of states
                              state, # State name as character string
                              overwrite="missing", # transient, missing, absorbing, all
                              dtms=NULL) {
  # Length
  nx <- length(x)
  # Find appearances
  whichfirst <- which(x==state)
  # Stop if no appearance
  if(length(whichfirst)==0) return(x)
  # Get first
  whichfirst <- min(whichfirst)
  # Replace: both missing and absorbing
  if(overwrite=="all") {
    dochange <- whichfirst:nx
    x[dochange] <- state
  }
  # Replace: only missing
  if(overwrite=="missing") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    dochange <- whichfirst:nx
    dochange <- setdiff(dochange,whichabsorbing)
    x[dochange] <- state
  }
  # Replace: only absorbing
  if(overwrite=="absorbing") {
    whichmissing <- which(is.na(x))
    dochange <- whichfirst:nx
    dochange <- setdiff(dochange,whichmissing)
    x[dochange] <- state
  }
  # Do not replace absorbing or missing
  if(overwrite=="transient") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    whichmissing <- which(is.na(x))
    whichboth <- union(whichabsorbing,whichmissing)
    dochange <- whichfirst:nx
    dochange <- setdiff(dochange,whichboth)
    x[dochange] <- state
  }
  # Return
  return(x)
}
dtms_backward_help <- function(x, # Vector of states
                               state, # State name as character string
                               overwrite="missing", # transient, missing, absorbing, all
                               dtms=NULL) {
  # Find appearances
  whichfirst <- which(x==state)
  # Stop if no appearance
  if(length(whichfirst)==0) return(x)
  # Get last
  whichfirst <- max(whichfirst)
  # Replace: both missing and absorbing
  if(overwrite=="all") {
    dochange <- 1:whichfirst
    x[dochange] <- state
  }
  # Replace: only missing
  if(overwrite=="missing") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    dochange <- 1:whichfirst
    dochange <- setdiff(dochange,whichabsorbing)
    x[dochange] <- state
  }
  # Replace: only absorbing
  if(overwrite=="absorbing") {
    whichmissing <- which(is.na(x))
    dochange <- 1:whichfirst
    dochange <- setdiff(dochange,whichmissing)
    x[dochange] <- state
  }
  # Do not replace absorbing or missing
  if(overwrite=="transient") {
    whichabsorbing <- which(x%in%dtms$absorbing)
    whichmissing <- which(is.na(x))
    whichboth <- union(whichabsorbing,whichmissing)
    dochange <- 1:whichfirst
    dochange <- setdiff(dochange,whichboth)
    x[dochange] <- state
  }
  # Return
  return(x)
}
dtms_duration_help <- function(states, # Vector of states of one unit
                               time, # Vector of time scale values
                               dtms, # dtms object
                               ignoreleft) { # TRUE or FALSE, as per dtms_duration
  # Get lengths
  lengths <- rle(states)
  duration <- unlist(lapply(lengths$length, function(x) 1:x))
  # Change if left censoring
  if(!ignoreleft & time[1]!=dtms$timescale[1]) duration[1:lengths[1]] <- NA
  # Remove gaps
  timediffs <- diff(time)
  timediffs <- c(FALSE,!timediffs%in%dtms$timestep)
  # Go through all gaps
  if(any(timediffs)) {
    # Find entries
    dropwhich <- which(timediffs)
    cumlengths <- cumsum(lengths$length)
    cumwhere <- c(1,cumlengths[-length(cumlengths)]+1)
    # Loop and replace
    for(drop in dropwhich) {
      # Which entries need to be replaced
      whichlengths <- which(cumlengths>=drop)[c(1,2)]
      from <- drop
      to <- cumwhere[whichlengths[2]]-1
      if(is.na(to)) to <- cumwhere[whichlengths[1]]
      # Replace
      duration[from:to] <- NA
    }
  }
  # Return
  return(duration)
}
dtms_occurrence_help <- function(states, # Vector of states of one unit
                                time, # Vector of time scale values
                                dtms, # dtms object
                                ignoreleft) { # TRUE or FALSE, as per dtms_occurence
  # Change if left censoring
  if(!ignoreleft & time[1]!=dtms$timescale[1]) {
    result <- rep(NA,length(states))
    return(result)
  }
  # Get states of spells
  spells <- rle(states)$values
  # Count occurence of spells
  if(any(is.na(spells))) occurences <- table(spells,useNA="always") else occurences <- table(spells)
  # Start counting from 1
  expanded_occurences <- lapply(occurences,function(x) 1:x)
  # Get names, vectorize
  names_expanded <- names(expanded_occurences)
  occurences <- unlist(expanded_occurences)
  # Match ordering
  ordering <- unlist(lapply(names_expanded,function(x) {
    if(is.na(x)) which(is.na(spells)) else which(x==spells)
  }))
  # Get counts in right order
  result <- numeric(length(occurences))
  result[ordering] <- occurences
  # Repeat counts for each observation in spell
  result <- inverse.rle(list(values=result,lengths=rle(states)$lengths))
  # Handle gaps: find if any
  timediffs <- diff(time)
  timediffs <- c(!timediffs%in%dtms$timestep)
  # If gap fill with NA starting from first gap
  if(any(timediffs)) {
      dropwhich <- min(which(timediffs))
      nresult <- length(result)
      result[(dropwhich+1):nresult] <- NA
    }
  # Return result
  return(result)
}
