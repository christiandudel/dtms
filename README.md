
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dtms

<!-- badges: start -->
<!-- badges: end -->

The package dtms implements discrete-time multistate models in R with a
user-friendly workflow. It also comes with many tools to analyze the
results of multistate models.

## Installation

You can install the development version of dtms like this:

``` r
library(devtools)
install_github("christiandudel/dtms")
```

## Example 1

This is a very basic example.

``` r
## Load package
library(dtms)

## Define model: Absorbing and transient states, time scale
# simple <- dtms(transient=c("A","B"),
#                absorbing="X",
#                time=0:20)
# 
# ## Quick look at data
# simpledata
# 
# ## Reshape to transition format
# estdata <- dtms_format(data=simpledata,
#                        dtms=simple,
#                        idvar="id",
#                        timevar="time",
#                        statevar="state")
# 
# ## Clean
# estdata <- dtms_clean(data=estdata,
#                       dtms=simple)
# 
# # Fit model 
# fit <- dtms_fit(data=estdata)
# 
# ## Predict probabilities
# probs    <- dtms_transitions(dtms=simple,
#                              model = fit)
# 
# ## Get transition matrix 
# Tp <- dtms_matrix(dtms=simple,
#                   probs=probs)
# 
# ## Get starting distribution 
# S <- dtms_start(dtms=simple,
#                 data=estdata)
# 
# ## State expectancies 
# dtms_expectancy(dtms=simple,
#                 matrix=Tp,
#                 start_distr=S)
# 
# ## Lifetime risk 
# dtms_risk(dtms=simple,
#           matrix=Tp,
#           risk="A")
# 
# ## Distribution of visits
# dtms_visits(dtms=simple,
#             matrix=Tp,
#             risk="A",
#             start_distr=S)
# 
# ## First/last visit
# dtms_first(dtms=simple,
#            matrix=Tp,
#            risk="A",
#            start_distr=S)
# 
# dtms_last(dtms=simple,
#           matrix=Tp,
#           risk="A",
#           start_distr=S,
#           rescale=T,
#           total=F)
```
