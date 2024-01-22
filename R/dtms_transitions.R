#' Title
#'
#' @param transient
#' @param absorbing
#' @param time
#' @param dtms
#' @param model
#' @param constant
#' @param varying
#' @param sep
#' @param timestep
#' @param timevar
#' @param fromvar
#' @param tovar
#' @param Pvar
#'
#' @return
#' @export
#'
#' @examples
dtms_transitions <- function(transient=NULL,
                             absorbing=NULL,
                             time=NULL,
                             dtms=NULL,
                             model,
                             constant=NULL,
                             varying=NULL,
                             sep="_",
                             timestep=NULL,
                             timevar="time",
                             fromvar="from",
                             tovar="to",
                             Pvar="P") {

  # Use dtms if provided
  if(!is.null(dtms) & class(dtms)[2]=="dtms") {
    if(is.null(transient)) transient <- dtms$transient
    if(is.null(absorbing)) absorbing <- dtms$absorbing
    if(is.null(time)) {
      time <- dtms$time
      time <- time[-length(time)]
    }
    if(is.null(timestep)) {
      timestep <- dtms$timestep
    }
  }

  # Get full state space
  all_states <- paste(c(transient,absorbing))

  # Create empty frame
  model_frame <- expand.grid(from=transient,
                             time=time)

  # Get names right
  names(model_frame) <- c(fromvar,timevar)

  # Add time-constant covariates
  nvar <- length(constant)
  if(nvar>0) {
    for(var in 1:nvar) {
      varname <- names(constant)[var] # Variable name
      value <- constant[[var]] # Value
      model_frame[varname] <- value # Assign
    }
  }

  # Add time-varying covariates
  nvar <- length(varying)
  if(nvar>0) {
    for(var in 1:nvar) {
      varname <- names(varying)[var] # Variable name
      value <- varying[[var]] # Value
      if(length(value)!=length(time)) stop("Wrong number of time-varying values")
      assign_values <- match(model_frame[,timevar],time) # Match to time variable
      model_frame[varname] <- value[assign_values] # Assign
    }
  }

  # Predict
  model_frame[,all_states] <- stats::predict(model,model_frame,"response")[,all_states]

  # Reshape
  model_frame <- model_frame |> tidyr::pivot_longer(cols=all_states,
                                              names_to=tovar,
                                              values_to=Pvar)

  # Get state names right with time
  model_frame <- model_frame |> dplyr::mutate(newfrom=ifelse(!!dplyr::sym(fromvar)%in%transient,
                                                       paste(!!dplyr::sym(fromvar),!!dplyr::sym(timevar),sep=sep),
                                                       !!dplyr::sym(fromvar)),
                                        newto=ifelse(!!dplyr::sym(tovar)%in%transient,
                                                     paste(!!dplyr::sym(tovar),!!dplyr::sym(timevar)+timestep,sep=sep),
                                                     !!dplyr::sym(tovar)))

  # Keep and rename variables
  model_frame <- model_frame[,c("newfrom","newto",timevar,Pvar)]
  names(model_frame)[1:2] <- c(fromvar,tovar)

  # Return
  return(model_frame)

}
