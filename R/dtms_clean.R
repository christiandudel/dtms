dtms_clean <- function(data, # Transition format data
                       dtms, # dtms object
                       fromvar="from", # Variable with starting state
                       tovar="to", # Variable with receiving state
                       timevar="time", # Variable name of time
                       dropTime=T, # Drop obs with time not in the model
                       dropNA=T, # Drop gaps, last obs, ...
                       dropAbs=T) { # Drop obs starting from absorbing state

  if(dropTime) {
    whichrows <- unlist(data[,timevar])%in%dtms$time
    count <- sum(!whichrows)
    cat("Dropping ",count," rows because not in time range\n")
    data <- data[whichrows,]
  }

  if(dropNA) {
    whichrows <- !is.na(unlist(data[,fromvar])) & !is.na(unlist(data[,tovar]))
    count <- sum(!whichrows)
    cat("Dropping ",count," rows because gap, last obs, ...\n")
    data <- data[whichrows,]
  }

  if(dropAbs) {
    whichrows <- !unlist(data[,fromvar])%in%dtms$absorbing
    count <- sum(!whichrows)
    cat("Dropping ",count," rows because starting in absorbing state\n")
    data <- data[whichrows,]
  }

}
