#' @export

### Create transition matrix without absorbing states

remove_absorbing <- function(matrix) { # matrix=full transition matrix

  to_remove <- which(diag(matrix)==1)
  removed <- matrix[-to_remove,-to_remove]
  return(removed)

}
