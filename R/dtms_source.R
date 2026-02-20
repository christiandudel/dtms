#' Source code only from specific lines
#'
#' @description
#' This is a helper function which allows to source only specific lines from
#' a file.
#'
#' @source User Matthew Plourde on Stackoverflow https://stackoverflow.com/questions/12214963/source-only-part-of-a-file
#'
#' @param file Name of the file to read (character).
#' @param start First line to read (numeric).
#' @param end Last line to read (numeric).
#' @param local Logical or an environment. Default is 'FALSE'. See function \code{source} for details.
#'
#' @seealso [func(source)]
#'
#' @return No return value, called for side effects
#' @export
#'
dtms_source <- function(file,start,end,local=FALSE) {

  file.lines <- scan(file=file,
                     what=character(),
                     skip=start-1,
                     nlines=end-start+1,
                     sep='\n')

  file.lines.collapsed <- paste(file.lines, collapse='\n')

  source(textConnection(file.lines.collapsed), local=local)

}
