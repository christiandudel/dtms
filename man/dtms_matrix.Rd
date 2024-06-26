% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dtms_matrix.R
\name{dtms_matrix}
\alias{dtms_matrix}
\title{Creates a transition matrix from transition probabilities}
\usage{
dtms_matrix(
  probs,
  dtms = NULL,
  fromvar = "from",
  tovar = "to",
  Pvar = "P",
  enforcedeath = T,
  rescale = T
)
}
\arguments{
\item{probs}{Data frame with transition probabilities, as created with \code{dtms_transitions}.}

\item{dtms}{dtms object, as created with \code{dtms}.}

\item{fromvar}{Character (optional), name of variable in `probs` with starting state. Default is "from".}

\item{tovar}{Character (optional), name of variable in `probs` with receiving state. Default is "to".}

\item{Pvar}{Character (optional), name of variable in `probs` with transition probabilities. Default is `P`.}

\item{enforcedeath}{Logical (optional), make sure that every unit moves to absorbing state after last value of time scale? Default is TRUE.}

\item{rescale}{Logical (optional), rescale transition probabilities to sum to 1? Default is TRUE.}
}
\value{
Returns a transition matrix.
}
\description{
This function creates a transiton matrix based on transition probabilities
predicted using the function `dtms_transitions`.
}
\examples{
simple <- dtms(transient=c("A","B"),
absorbing="X",
timescale=0:20)
## Reshape to transition format
estdata <- dtms_format(data=simpledata,
                       dtms=simple,
                       idvar="id",
                       timevar="time",
                       statevar="state")
## Clean
estdata <- dtms_clean(data=estdata,
                      dtms=simple)
## Fit model
fit <- dtms_fit(data=estdata)
## Predict probabilities
probs    <- dtms_transitions(dtms=simple,
                             model = fit)
## Get transition matrix
Tp <- dtms_matrix(dtms=simple,
                  probs=probs)
}
