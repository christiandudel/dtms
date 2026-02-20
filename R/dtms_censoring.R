#' Left censoring, right censoring, and gaps in data
#'
#' @description
#' This function provides an overview of censoring and gaps in the data. It can
#' do so in several ways: by providing counts of units with left censoring,
#' right censoring, and gaps; by providing a cross-tabulation of the number of
#' units with left censoring and/or right censoring and/or gaps; and by
#' returning a data frame with added indicators on censoring and gaps.
#'
#' @details
#' Added variables can be at the unit level or at the observation level. This
#' is controlled by the argument "addtype". If it is set to "id" then the unit
#' level is used. In this case the added variables are the same for
#' each observation  of a unit. For instance, if a unit experiences any gap,
#' then the added variable has the value TRUE for all observations of that unit.
#' If "addtype" is set to "obs" the observation level is used and the indicators
#' are only set to TRUE if they apply to a specific observation. For instance,
#' if a unit experience right censoring, only the last observation will have
#' TRUE as the value for the right-censoring indicator; i.e., showing that after
#' this last observation there is right censoring. This can be helpful for
#' analyses to understand censoring better.
#'
#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable in `data` with starting state. Default is "from".
#' @param tovar Character (optional), name of variable in `data` with receiving state. Default is "to".
#' @param timevar Character (optional), name of variable in `data` with time scale. Default is "time".
#' @param idvar Character (optional), name of variable in `data` with unit ID. Default is "id".
#' @param print Logical (optional), print counts? Default is TRUE.
#' @param printlong Logical (optional), print cross-tabulation? Default is FALSE.
#' @param add Logical (optional), add indicators to data set? Default is FALSE. If TRUE the data frame specified with \code{data} is returned with added columns.
#' @param addtype Character (optional), what type of information should be added if add=TRUE. Either "id" or "obs", see details. Default is "id".
#' @param varnames Character vector (optional), names of added variables if add=TRUE. Default is "c("LEFT","GAP","RIGHT")".
#'
#' @return Table or data frame.
#' @export
#'
#' @examples
#' ## Define model: Absorbing and transient states, time scale
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:19)
#' # Reshape to transition format
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' ## Clean
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' ## Censoring
#' dtms_censoring(data=estdata,
#'                dtms=simple)

dtms_censoring <- function(data,
                           dtms,
                           fromvar="from",
                           tovar="to",
                           timevar="time",
                           idvar="id",
                           print=TRUE,
                           printlong=FALSE,
                           add=FALSE,
                           addtype="id",
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

    # Left censoring
    left <- by(data[,timevar],
               data[,idvar],
               FUN=function(x) {
                 first <- x[1]>min(dtms$timescale)
                 rest <- rep(FALSE,length(x)-1)
                 return(c(first,rest))
                 })

    # GGaps
    gap <- by(data[,timevar],
              data[,idvar],
              FUN=function(x) c(diff(x)!=dtms$timestep,FALSE))

    # Right censoring: Last state is not absorbing
    right1 <- by(data[,tovar],
                 data[,idvar],
                 FUN=function(x) {
                   last <- !x[length(x)]%in%dtms$absorbing
                   rest <- rep(FALSE,length(x)-1)
                   return <- c(rest,last)
                   })

    # Right censoring: Last obs is not at end of timescale
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
