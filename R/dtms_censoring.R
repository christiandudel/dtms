#' Get information on left censoring, right censoring, and gaps in data
#'
#' @param data Data frame
#' @param dtms dtms object
#' @param fromvar Variable with starting state
#' @param tovar Variable with receiving state
#' @param timevar Variable with time scale
#' @param idvar Variable with unit ID
#' @param add Add indicators to data set?
#' @param addtype If add=T, what type of information should be added? Either "id" or "obs".
#' @param print Print short table?
#' @param printlong Print long table?
#' @param varnames If variables are added, what should be their names
#'
#' @return Table or data frame
#' @export

dtms_censoring <- function(data,
                           dtms,
                           fromvar="from",
                           tovar="to",
                           timevar="time",
                           idvar="id",
                           add=F,
                           addtype="id",
                           print=T,
                           printlong=F,
                           varnames=c("LEFT","GAP","RIGHT")) {

  # Check dtms
  dtms_proper(dtms)

  # Sort data
  dataorder <- order(data[,idvar],
                     data[,timevar])
  data <- data[dataorder,]

  # Check for left censoring
  left <- by(data[,timevar],
             data[,idvar],
             FUN=function(x) min(x)>min(dtms$timescale))

  # Check for gaps
  gap <- by(data[,timevar],
            data[,idvar],
            FUN=function(x) any(diff(x)!=dtms$timestep))

  # Check for right censoring
  right1 <- by(data[,tovar],
              data[,idvar],
              FUN=function(x) !x[length(x)]%in%dtms$absorbing)

  right2 <- by(data[,timevar],
               data[,idvar],
               FUN=function(x) !max(x)==max(dtms$timescale))

  # Vectorize and combine
  left <- as.logical(left)
  gap <- as.logical(gap)
  right <- as.logical(right1) & as.logical(right2)

  # Print
  if(print) {
    cat("Units with left censoring: ",sum(left),"\n")
    cat("Units with gaps: ",sum(gap),"\n")
    cat("Units with right censoring: ",sum(right),"\n")
  }
  if(printlong) {
    cat("Cross tabulation:","\n")
    print(table(left,right,gap))
  }

  # Add indicators to data by ID
  if(add & addtype=="id") {

    # Frame for merging
    tmp <- data.frame(id=unique(data[,idvar]),
                      left=left,
                      gap=gap,
                      right=right)

    names(tmp) <- c(idvar,varnames)

    # Merge with data
    data <- merge(data,tmp)

    # Return
    return(data)

  }

  # Add indicators to data by observation
  if(add & addtype=="obs") {

    # Getting vector for left censoring
    left <- by(data[,timevar],
               data[,idvar],
               FUN=function(x) {
                 first <- x[1]>min(dtms$timescale)
                 rest <- rep(FALSE,length(x)-1)
                 return(c(first,rest))
                 })

    # Getting vector for gaps
    gap <- by(data[,timevar],
              data[,idvar],
              FUN=function(x) c(diff(x)!=dtms$timestep,FALSE))

    # Getting vectors for right censoring
    right1 <- by(data[,tovar],
                 data[,idvar],
                 FUN=function(x) {
                   last <- !x[length(x)]%in%dtms$absorbing
                   rest <- rep(FALSE,length(x)-1)
                   return <- c(rest,last)
                   })

    right2 <- by(data[,timevar],
                 data[,idvar],
                 FUN=function(x) {
                   last <- !max(x)==max(dtms$timescale)
                   rest <- rep(FALSE,length(x)-1)
                   return <- c(rest,last)
                   })

    # Vectorize and combine
    left <- unlist(left)
    gap <- unlist(gap)
    right <- unlist(right1) & unlist(right2)

    # Put into data
    data[,varnames[1]] <- left
    data[,varnames[2]] <- gap
    data[,varnames[3]] <- right

    # Return
    return(data)

  }

}
