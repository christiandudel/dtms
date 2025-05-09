% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dtms_visits.R
\name{dtms_visits}
\alias{dtms_visits}
\title{Calculate the distribution of the time spent in a subset of states}
\usage{
dtms_visits(
  matrix,
  dtms,
  risk,
  start_time = NULL,
  start_state = NULL,
  start_distr = NULL,
  end_time = NULL,
  method = "mid",
  total = F
)
}
\arguments{
\item{matrix}{Matrix with transition probabilities, as generated with \code{dtms_matrix}.}

\item{dtms}{dtms object, as created with \code{dtms}.}

\item{risk}{Character, name of state(s) for which risk is of interest.}

\item{start_time}{Numeric (optional), value of time scale for start. If NULL (default) first value of time scale will be used.}

\item{start_state}{Character (optional), name of starting states. If NULL (default) all transient states will be used.}

\item{start_distr}{Numeric (optional), distribution of starting states. If specified, average distribution over all starting states will be calculated.}

\item{end_time}{Numeric (optional), last value of time scale to consider. If NULL (default) all values of time scale starting from start_time will be used.}

\item{method}{Character (optional), do transitions happen mid-interval (`mid`, default) or at the end of the interval (`end`), see details.}

\item{total}{Logical, should total of distribution be shown (always sums to 1)? Default is FALSE.}
}
\value{
A table with the distribution of time spent in a subset of states.
}
\description{
Calculates the distribution of the time spent in a state or a subset of states.
}
\details{
The state(s) which count to the time are specified with the argument `risk`.
If several states are specified, the resulting distribution refers to the
lifetime spent in any of the specified states.

In a discrete-time model, the time spent in a state depends on assumptions
about when transitions happen. Currently, this functions supports two
variants which can be specified with the argument `method`: mid-interval
transitions can be selected with the option `mid` and imply that transitions
happen at the middle of the time interval; and the option `end` assumes
that instead transitions happen at the end of the interval. In this latter
case the distribution of the time spent in a state is equivalent to the
number of visits to that state. The calculation takes the step length of
the time scale into account as specified by the `dtms` object. If the
step length is not one fixed value, the first entry of `dtms$timestep` will
be used.

If a distribution of the starting states is provided with `start_distr` the
output table has two additional rows. One shows the distribution
unconditional on the starting state. The other shows the distribution
conditional on not starting in any state of the risk set.
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
## Get transition matrix
Tp <- dtms_matrix(dtms=simple,
                  probs=probs)
## Get starting distribution
S <- dtms_start(dtms=simple,
                data=estdata)
## Distribution of visits
dtms_visits(dtms=simple,
            matrix=Tp,
            risk="A",
            start_distr=S,
            total=TRUE)
}
\seealso{
\code{\link{dtms_distr_summary}} to help with summarizing the resulting distribution.
}
