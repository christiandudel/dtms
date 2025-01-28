#' Carry states forward
#'
#' @description
#' This function carries a state forward after its first occurrence.
#'
#' @details
#' This function carries a state forward after its first occurrence.
#' For instance, carrying the state "A" forward in the sequence `B, B, A, B, B`
#' will give the sequence `B, B, A, A, A`. The sequence `C, B, C, A, B, A, A, B`
#' will give `C, B, C, A, A, A, A, A`.
#'
#' This function works with data frames in transition format and in long format.
#' The default is transition format, using the arguments `fromvar` and `tovar`.
#' If, however, the argument `statevar` is specified, it is used instead.
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
#' @param fromvar Character (optional), name of variable with starting state. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state. Default is `to`.
#' @param statevar Character (optional), name of the variable in the data frame in long format with the states. Default is NULL.
#' @param idvar Character (optional), name of variable with unit ID. Default is `id`.
#' @param timevar Character (optional), name of variable with time scale. Default is `time`.
#' @param dtms dtms object (optional), as created with \code{dtms}. Not required if `overwrite==transient`.
#' @param overwrite Character (optional), one of `transient`, `missing`, `absorbing`, and `all`, see details. Default is `transient`.
#' @param vector Logical (optional), return vector (if TRUE) or data frame (if FALSE). Default is FALSE. Argument is only used if argument `statevar` is specified.
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
#'              statevar="state",
#'              state="A",
#'              dtms=simple,
#'              overwrite="transient")

dtms_forward <- function(data,
                         state,
                         fromvar="from",
                         tovar="to",
                         statevar=NULL,
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

  # Apply helper for transition format
  if(is.null(statevar)) {

    # Move starting state
    data[,tovar][data[,fromvar]==state] <- state

    # Carry forward starting state
    tmp1 <- tapply(data[,fromvar],
                   data[,idvar],
                   function(x) dtms_forward_help(x=x,
                                                 state=state,
                                                 overwrite=overwrite,
                                                 dtms=dtms))
    # Carry forward receiving state
    tmp2 <- tapply(data[,tovar],
                   data[,idvar],
                   function(x) dtms_forward_help(x=x,
                                                 state=state,
                                                 overwrite=overwrite,
                                                 dtms=dtms))

    # Assign new values
    data[,fromvar] <- unlist(tmp1)
    data[,tovar] <- unlist(tmp2)


  } else {

    # Apply helper for long format
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

  }

  # Return
  return(data)

}
