#' Aggregate data
#'
#' @description
#' This function takes any data set and returns a new data frame which only
#' includes the unique rows from the original data set and indicates how
#' often these rows appear in the original data.
#'
#' @details
#' Currently, missing values are not supported and will be dropped; consider
#' using factors if you want to keep them. The variable provided with the
#' argument `idvar` will be dropped from the aggregated data. If `weights` is
#' specified, the counts will be placed in a variable with the same name. If
#' `countvar` is used, any existing variable in the original data with the
#' same name will be replaced.
#'
#' @param data Data frame.
#' @param weights Character (optional). Name of variable with weights.
#' @param idvar Character (optional). Name of variable in `data` with unit ID. Default is "id".
#' @param countvar Character (optional). Name of new variable in data which provides the counts. Default is "count".
#'
#' @returns An aggregated data frame
#' @export
#'
#' @examples
#' ## Define model: Absorbing and transient states, time scale
#' simple <- dtms(transient=c("A","B"),
#'                absorbing="X",
#'                timescale=0:20)
#' ## Reshape to transition format
#' estdata <- dtms_format(data=simpledata,
#'                        dtms=simple,
#'                        idvar="id",
#'                        timevar="time",
#'                        statevar="state")
#' ## Clean
#' estdata <- dtms_clean(data=estdata,
#'                       dtms=simple)
#' ## Aggregate
#' aggdata <- dtms_aggregate(estdata)
#' ## Fit model
#' fit <- dtms_fit(data=aggdata,
#'                 weights="count")
#' ## Basic resampling (can be extended for basic bootstrap)
#' aggdata$newcount <- rmultinom(1,size=sum(aggdata$count),prob=aggdata$count)
#' fit <- dtms_fit(data=aggdata,
#'                 weights="count")

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
