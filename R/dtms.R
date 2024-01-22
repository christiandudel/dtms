#' Create dtms object
#'
#' @param transient A character vector of names of the transient states in the state space
#' @param absorbing A character vector of names of the absorbing states in the state space
#' @param time A numeric vector with the time scale, including the starting time and the final time
#' @param timestep Step length of the time scale, will be guessed if NULL
#'
#' @return Returns an object of class 'dtms' to be passed to other functions of the package
#' @export
#'
#' @examples
#' dtms(transient=c("A","B"),
#'      absorbing="X",
#'      time=1:10)

dtms <- function(transient,
                 absorbing,
                 time,
                 timestep=NULL) {

  # Guess time step?
  if(is.null(timestep)) {

    # Step lengths as specified
    step <- time |> diff() |> unique()

    # Number of different step lengths
    nstep <- length(step)

    # Guess timestep if possible, else error
    if(nstep==1) timestep <- step else
      stop("Not able to guess step length of time scale")

  }

  # Combine everything in a list
  result <- list(transient=transient,
                 absorbing=absorbing,
                 time=time,
                 timestep=timestep)

  # Assign class
  class(result)[2] <- "dtms"

  # Return
  return(result)

}
