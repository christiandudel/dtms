#' Carry states forward
#'
#' @description
#' This function carries a state forward after its first occurrence.
#'
#' @details
#' This function carries a state forward after its first occurrence.
#' For instance, carrying the state "A" forward in the sequence `B, B, A, B, B`
#' will give the sequence `B, B, A, A, A`.
#'
#' The data frame has to be in long format and not in transition format; i.e.,
#' usually, this function will be applied before \code{dtms_format}.
#'
#' The argument `overwrite` is used to control what type of information is
#' replaced. If `overwrite==transient`, then only transient states are replaced
#' while missing values and absorbing states remain unchanged. For example,
#' carrying forward state "A" in the sequence `B, B, A, B, NA, X, X` with X
#' being an absorbing state will give `B, B, A, A, NA, X, X`. If
#' `overwrite==missing` then in addition to transient states also missing values
#' are replaced and for the example sequence `B, B, A, A, A, X, X` would be
#' returned. If `overwrite==absorbing` then in addition to transient states
#' absorbing states will be replaced; for the example sequence the result would
#' be `B, B, A, A, NA, A, A`. Finally, if `overwrite==all` then all values in
#' the sequence will be replaced: `B, B, A, A, A, A, A`.
#'
#' @seealso [func(dtms_backward)]
#'
#' @param data A data frame in long format.
#' @param state Character, name of the state to be carried forward.
#' @param statevar Character (optional), name of the variable in the data frame with the states. Default is "state".
#' @param idvar Character (optional), name of variable with unit ID. Default is `id`.
#' @param timevar Character (optional), name of variable with time scale. Default is `time`.
#' @param dtms dtms object (optional), as created with \code{dtms}. Not required if `overwrite==transient`.
#' @param overwrite Character (optional), one of `transient`, `missing`, `absorbing`, and `all`, see details. Default is `transient`.
#' @param vector Logical (optional), return vector (if TRUE) or data frame (if FALSE). Default is FALSE
#'
#' @return The data frame specified with `data` and the edited state variable (if `vector=FALSE`) or a vector (if `vector=TRUE`).
#' @export
#'
#' @examples
#' simple <- dtms(transient=c("A","B"),
#' absorbing="X",
#' timescale=0:19)
#'
#' dtms_forward(data=simpledata,
#'              state="A",
#'              dtms=simple,
#'              overwrite="transient")

dtms_forward <- function(data,
                         state,
                         statevar="state",
                         idvar="id",
                         timevar="time",
                         dtms=NULL,
                         overwrite="missing",
                         vector=FALSE) {

  # Check dtms only if specified
  if(!is.null(dtms)) dtms_proper(dtms)

  # Sort data
  dataorder <- order(data[,idvar],
                     data[,timevar])
  data <- data[dataorder,]

  # Apply helper
  tmp <- tapply(data[,statevar],
                data[,idvar],
                function(x) dtms_forward_help(x=x,
                                              state=state,
                                              overwrite=overwrite,
                                              dtms=dtms))

  # Return vector
  if(vector) return(unlist(tmp))

  # Assign new values
  data[,statevar] <- unlist(tmp)

  # Return
  return(data)

}
