#' @param data Data frame in transition format, as created with \code{dtms_format}.
#' @param dtms dtms object, as created with \code{dtms}.
#' @param fromvar Character (optional), name of variable with starting state. Default is `from`.
#' @param tovar Character (optional), name of variable with receiving state. Default is `to`.
#' @param timevar Character (optional), name of variable with time scale. Default is `time`.
#' @param idvar Character (optional), name of variable with unit ID. Default is `id`.
#' @param print Logical (optional), print counts? Default is TRUE.
#' @param printlong Logical (optional), print cross-tabulation? Default is FALSE.
#' @param add Logical (optional), add indicators to data set? Default is FALSE. If TRUE the data frame specified with \code{data} is returned with added columns.
#' @param addtype Character (optional), what type of information should be added if add=TRUE. Either `id` or `obs`, see details. Default is `id`.
#' @param varnames Character vector (optional), names of added variables if add=T. Default is `c("LEFT","GAP","RIGHT")`.
#' @param distr An object of class "dtms_distr" created with \code{dtms_visits}, \code{dtms_first}, or \code{dtms_last}.
#' @param weightvar Character (optional), name of variable with weights. Default is NULL.
#' @param matrix Matrix with transition probabilities, as generated with \code{dtms_matrix}.
#' @param risk Character (otpional), name of one transient state. If specified expectancies are only shown for this state but by values of the time scale.
#' @param risk Character, name of state(s) for which risk is of interest.
#' @param start_distr Numeric (optional), distribution of starting states. If specified, average expectancy over all starting states will be calculated. Only applied if risk=NULL.
#' @param start_state Character (optional), name of starting states. If NULL (default) all transient states will be used.
#' @param start_time Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.
#' @param end_time Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.
#' @param correction Numeric (optional), correction for expectancy when starting state and state under consideration match, see details. Defaults to 0.5.
#' @param total Logical (optional), calculate total expectancy. Default is TRUE. Only applied if \code{risk=NULL}.
#' @param verbose Logical (optional), print some information on what is computed. Default is FALSE.
#' @param fundamental Logical (optional), return fundamental matrix? Default is FALSE.
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#' @param method Character (optional), do transitions happen mid-interval (`mid`, default) or at the end of the interval (`end`), see details.
#' @param rescale Logical (optional), should distribution be rescaled to sum to 1? See details. Default is TRUE.
#' @param Pvar Character (optional), name of variable with transition probabilities. Default is `P`.
#' @param enforcedeath Logical (optional), make sure that every unit moves to absorbing state after last value of time scale? Default is TRUE.
#' @param probs Object with transition probabilities as created with \code{dtms_transitions}.
#' @param variables List (optional), a named list with covariate values which are used to restrict the data.
#' @param model Model estimated with \code{dtms_fit}.
