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
  aggformula <- as.formula(aggformula)

  # Warning if missing values
  drops <- data |>
    na.omit() |>
    attributes(x=_) |>
    getElement(object=_,"na.action") |>
    length()
  if(drops>0) warning(paste("Dropping",drops,"rows with missing values"))

  # Aggregate
  tmp <- aggregate(by=aggformula,
                   x=data,
                   FUN=sum)

  # Return
  return(tmp)

}
