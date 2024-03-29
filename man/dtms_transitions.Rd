% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dtms_transitions.R
\name{dtms_transitions}
\alias{dtms_transitions}
\title{Predict transition probabilities}
\usage{
dtms_transitions(
  model,
  dtms,
  constant = NULL,
  varying = NULL,
  timevar = "time",
  fromvar = "from",
  tovar = "to",
  Pvar = "P"
)
}
\arguments{
\item{model}{Model estimated with \code{dtms_fit}.}

\item{dtms}{dtms object, as created with \code{dtms}.}

\item{constant}{List (optional) with values for time-constant predictors (see details).}

\item{varying}{List (optional) with values for time-varying predictors (see details).}

\item{timevar}{Character (optional), name of variable with time scale. Default is `time`.}

\item{fromvar}{Character (optional), name of variable with starting state. Default is `from`.}

\item{tovar}{Character (optional), name of variable with receiving state. Default is `to`.}

\item{Pvar}{Character (optional), name of variable with transition probabilities. Default is `P`.}
}
\value{
A data frame with transition probabilities.
}
\description{
`dtms_transitions` predicts transition probabilities based on a model
estimated with `dtms_fit`.
}
\details{
Predicted transition probabilities are returned as a data frame, and not
as a transition matrix. While the latter is required for applying Markov
chain methods, the data frame is more convenient for viewing and
analyzing the transition probabilities themselves.

Depending on the model specification, the prediction of transition
probabilities will require values for predictor variables which can be
specified with the arguments `constant` and `varying` for time-constant and
time-varying variables, respectively. In both cases, named lists have to be
used, where each entry name must correspond to a variable name in the model.
For time-constant variables, each list entry is of length one and provides
a value for the corresponding time-constant variable. For time-varying
variables, each entry must have the length of the time scale minus one, and
provide a value for each (potential) transition in the model; i.e., starting
from time t=0, starting from time t=1, etc., until time t=T-1.

The argument `separator` sets the separator used to create state names. State
names are either a combination of the name of a transient state and a value
of the time scale, or the name of an absorbing state.
}
\examples{
## Define model: Absorbing and transient states, time scale
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
}
