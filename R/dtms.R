#' Create dtms object
#'
#' @description
#' This function creates an object of class 'dtms' to be passed to other
#' functions of the package.
#'
#' @details
#' \code{dtms} provides an abstract definition of a multistate model, including
#' the names of the transient states, the names of the absorbing states, the
#' values the time scale can take, and the step length of the time scale.
#'
#' The names of the absorbing and transient states should be provided as
#' character strings. However, numeric values also work. Factors are not
#' supported
#'
#' @param transient A character vector of names of the transient states in the state space.
#' @param absorbing A character vector of names of the absorbing states in the state space.
#' @param timescale A numeric vector with the time scale, including the starting time and the final time.
#' @param timestep Numeric (optional), step length of the time scale, will be guessed if NULL (default).
#' @param sep Character (optional), separator between short state name and value of time scale. Default is `_`.
#'
#' @return Returns an object of class 'dtms'
#' @export
#'
#' @examples
#' dtms(transient=c("A","B"),
#'      absorbing="X",
#'      timescale=1:10)

dtms <- function(transient,
                 absorbing,
                 timescale,
                 timestep=NULL,
                 sep="_") {

  # Guess time step?
  if(is.null(timestep)) {

    # Step lengths as specified by time scale
    step <- timescale |> diff() |> unique()

    # Number of different step lengths
    nstep <- length(step)

    # Guess timestep if possible, else error
    if(nstep==1) timestep <- step else
      stop("Not able to guess step length of time scale")

  }

  # Combine everything in a list
  result <- list(transient=transient,
                 absorbing=absorbing,
                 timescale=timescale,
                 timestep=timestep,
                 sep=sep)

  # Assign class
  class(result)[2] <- "dtms"

  # Return
  return(result)

}
