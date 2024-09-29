#' Reshape data to transition format
#'
#' @description
#' Takes a data frame in long format and reshapes it into transition format.
#'
#' @details
#' The data frame supplied with the `data` argument has to be in long format,
#' where X is a time-constant covariate and Z(t) is a time-dependent covariate:
#'
#' \tabular{lllll}{
#' idvar \tab timevar \tab statevar \tab X \tab Z(t) \cr
#' 1 \tab 0 \tab A \tab x_1 \tab z_1(0) \cr
#' 1 \tab 1 \tab A \tab x_1 \tab z_1(1) \cr
#' 1 \tab 2 \tab B \tab x_1 \tab z_1(2) \cr
#' 1 \tab 3 \tab A \tab x_1 \tab z_1(3) \cr
#' 2 \tab 0 \tab B \tab x_2 \tab z_2(0) \cr
#' 2 \tab 1 \tab A \tab x_2 \tab z_2(1) \cr
#' ... \tab ... \tab ... \tab ... \tab ...
#' }
#'
#' If it is not in long format it has to be reshaped. The state variable
#' should provide the states as character strings or numbers; factors are not
#' supported.
#'
#' `dtms_format` turns the data set above into a data frame in transition
#' format:
#'
#' \tabular{llllll}{
#' id \tab time \tab fromvar \tab tovar \tab X \tab Z(t) \cr
#' 1 \tab 0 \tab A \tab A \tab x_1 \tab z_1(0) \cr
#' 1 \tab 1 \tab A \tab B \tab x_1 \tab z_1(1)  \cr
#' 1 \tab 2 \tab B \tab A \tab x_1 \tab z_1(2)  \cr
#' 2 \tab 0 \tab B \tab A \tab x_2 \tab z_2(0)  \cr
#' ... \tab ... \tab ... \tab ... \tab ... \tab ...
#' }
#'
#' Covariates do not need to be specified and are handled implicitly. The
#' transition from time t to t+1 takes covariate values from time t. By default
#' the variable names of the ID variable and the time variable are changed to
#' `id` and `time`, as the other functions of the package use these as default
#' names. If renaming the variables is not possible because these variable
#' names already appear in the data then the original names are used. If
#' `keepnames=TRUE` the original names for `id` and `time` are kept.
#'
#' `dtms_format` by default drops gaps in the data, as no transitions are
#' observed. For instance, in the following example there is no observation at
#' time 4, and thus no transition is observed from t=3 to t=4; and no
#' transition from t=4 to t=5:
#'
#' \tabular{lll}{
#' idvar \tab timevar \tab statevar \cr
#' 1 \tab 0 \tab A \cr
#' 1 \tab 1 \tab A \cr
#' 1 \tab 2 \tab B \cr
#' 1 \tab 3 \tab A \cr
#' 1 \tab 5 \tab A \cr
#' }
#'
#' In this example, `dtms_format` will return the following:
#'
#'  \tabular{llll}{
#' id \tab time \tab fromvar \tab tovar \cr
#' 1 \tab 0 \tab A \tab A \cr
#' 1 \tab 1 \tab A \tab B \cr
#' 1 \tab 2 \tab B \tab A \cr
#' }
#'
#' If `fill=T`, then `dtms_format` will return the following:
#'
#'  \tabular{llll}{
#' id \tab time \tab fromvar \tab tovar \cr
#' 1 \tab 0 \tab A \tab A \cr
#' 1 \tab 1 \tab A \tab B \cr
#' 1 \tab 2 \tab B \tab A \cr
#' 1 \tab 3 \tab A \tab NA \cr
#' 1 \tab 4 \tab NA \tab A \cr
#' }
#'
#' The argument \code{absorbing} controls if the first observed absorbing state
#' is carried over to later observations. For instance, if a unit is first
#' observed to be in transient state 'A' at time 1, then in absorbing state 'X'
#' at time 2, and then in transient state 'A' at time 3, \code{absorbing=TRUE}
#' will lead to replacement of the state at time 3 with 'X'.
#'
#' @param data Data frame in long format.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param idvar Character (optional), name of variable in `data` with unit ID. Default is "id".
#' @param timevar Character (optional), name of variable in `data` with time scale. Default is "time".
#' @param statevar Character (optional), name of variable in `data` with state, default is 'state'.
#' @param fromvar Character (optional), name of variable`with starting state in reshaped data. Default is "from".
#' @param tovar Character (optional), name of variable with receiving state in reshaped data. Default is "to".
#' @param absorbing Logical (optional), use first observed absorbing state consistently? See details. Default is TRUE.
#' @param keepnames Logical (optional), keep original names for id and time variable? Default is FALSE; i.e., not keeping original names.
#' @param fill Logical (optional), fill implicit missing values with explicit NA? Default is FALSE.
#' @param verbose Logical (optional), create output to console if changing variable names is not possible? Default is TRUE.
#' @param steplength Logical (optional), if true, the time to the next state is returned as a variable. Default is FALSE.
#' @param stepvar Character (optional), if `steplength==TRUE`, this specifies the name of the variable with the step length. Default is `steplength`.
#'
#' @return A data set reshaped to transition format
#' @export
#'
#' @examples
#' # Define model
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:20)
#' # Transiton format
#' estdata <- dtms_format(data=simpledata,
#' dtms=simple,
#' idvar="id",
#' timevar="time",
#' statevar="state")

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
