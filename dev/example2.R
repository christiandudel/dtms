## Load package
library(dtms)
library(ggplot2)
library(devtools)
source("R/dtms_helpers.R")

## Define model: Absorbing and transient states, time scale
hrs <- dtms(transient=c("Working","Non-working","Retired"),
            absorbing="Dead",
            timescale=50:99)

## Quick look at data
head(hrsdata)

## Reshape
estdata <- dtms_format(data=hrsdata,
                       dtms=hrs,
                       idvar="ID",
                       timevar="Age",
                       statevar="State")

## Drop dead-to-dead transitions etc
estdata <- dtms_clean(data=estdata,
                      dtms=hrs)

## Overview
summary(estdata)

## Basic censoring
dtms_censoring(data=estdata,
               dtms=hrs)

## More advanced censoring example
estdata <- dtms_censoring(data=estdata,
                          dtms=hrs,
                          add=T,
                          addtype="obs")

estdata |>
  subset(subset=to!="Dead",select=c(RIGHT,to)) |>
  table() |>
  prop.table(margin=1)

## Add age squared
estdata$time2 <- estdata$time^2

## Fit model
fit <- dtms_fit(data=estdata,
                controls=c("Gender","time2"),
                package="mclogit")

## Transition probabilities by gender

# Men
probs_m <- dtms_transitions(dtms=hrs,
                            model = fit,
                            constant = list(Gender=0),
                            varying = list(time2 = (50:98)^2),
                            CI=TRUE)

# Women
probs_w <- dtms_transitions(dtms=hrs,
                            model = fit,
                            constant = list(Gender=1),
                            varying = list(time2 = (50:98)^2))

# Overview
summary(probs_m)
summary(probs_w)

# Plotting, men as example
probs_m |>  dtms_simplify() |>
  ggplot(aes(x=time,y=P,color=to)) +
  geom_ribbon(aes(ymin = CIlow, ymax = CIup,fill=to),alpha=0.5) +
  geom_line() +
  facet_wrap(~from)


## Transition matrices
Tm <- dtms_matrix(dtms=hrs,
                  probs=probs_m)

Tw <- dtms_matrix(dtms=hrs,
                  probs=probs_w)

## Starting distributions
Sm <- dtms_start(dtms=hrs,
                 data=estdata,
                 variables=list(Gender=0))

Sw <- dtms_start(dtms=hrs,
                 data=estdata,
                 variables=list(Gender=1))

## State expectancies
dtms_expectancy(dtms=hrs,
                matrix=Tm,
                start_distr=Sm)

dtms_expectancy(dtms=hrs,
                matrix=Tw,
                start_distr=Sw)

## A variant: ignoring retirement as a starting state (shown only for men)
limited <- c("Working","Non-working")

Smwr <- dtms_start(dtms=hrs,
                   data=estdata,
                   start_state=limited,
                   variables=list(Gender=0))

dtms_expectancy(dtms=hrs,
                matrix=Tm,
                start_state=limited,
                start_distr=Smwr)

## Lifetime risk of reaching retirement
dtms_risk(dtms=hrs,
          matrix=Tm,
          risk="Retired",
          start_distr=Sm)

dtms_risk(dtms=hrs,
          matrix=Tw,
          risk="Retired",
          start_distr=Sw)

## Distribution of visits
visitsm <- dtms_visits(dtms=hrs,
                       matrix=Tm,
                       risk="Retired",
                       start_distr=Sm)

visitsw <- dtms_visits(dtms=hrs,
                       matrix=Tw,
                       risk="Retired",
                       start_distr=Sw,
                       method="end")

summary(visitsm)
summary(visitsw)

## First visit
firstm <- dtms_first(dtms=hrs,
                     matrix=Tm,
                     risk="Retired",
                     start_distr=Sm)

firstw <- dtms_first(dtms=hrs,
                     matrix=Tw,
                     risk="Retired",
                     start_distr=Sw)

summary(firstm)
summary(firstw)

## Last exit

# Leaving employment to any state
last1m <- dtms_last(dtms=hrs,
                    matrix=Tm,
                    risk="Working",
                    start_distr=Sm)

last1w <- dtms_last(dtms=hrs,
                    matrix=Tw,
                    risk="Working",
                    start_distr=Sw)

summary(last1m)
summary(last1w)

# Leaving employment for retirement
last2m <- dtms_last(dtms=hrs,
                    matrix=Tm,
                    risk="Working",
                    risk_to="Retired",
                    start_distr=Sm)

last2w <- dtms_last(dtms=hrs,
                    matrix=Tw,
                    risk="Working",
                    risk_to="Retired",
                    start_distr=Sw)

summary(last2m)
summary(last2w)

# Bootstrap example
bootfun <- function(data,dtms) {

  fit <- dtms_fit(data=data,
                  controls=c("Gender","time2"),
                  package="mclogit")

  probs_m <- dtms_transitions(dtms=dtms,
                              model = fit,
                              constant = list(Gender=0),
                              varying = list(time2 = (50:98)^2),
                              CI=TRUE)

  probs_w <- dtms_transitions(dtms=dtms,
                              model = fit,
                              constant = list(Gender=1),
                              varying = list(time2 = (50:98)^2))

  Tm <- dtms_matrix(dtms=dtms,
                    probs=probs_m)

  Tw <- dtms_matrix(dtms=dtms,
                    probs=probs_w)

  Sm <- dtms_start(dtms=dtms,
                   data=data,
                   variables=list(Gender=0))

  Sw <- dtms_start(dtms=dtms,
                   data=data,
                   variables=list(Gender=1))

  res1 <- dtms_expectancy(dtms=dtms,
                  matrix=Tm,
                  start_distr=Sm)

  res2 <- dtms_expectancy(dtms=dtms,
                  matrix=Tw,
                  start_distr=Sw)

  rbind(res1,res2)

}

bootresults <- dtms_boot(data=estdata,
                         dtms=hrs,
                         fun=bootfun,
                         idvar="id",
                         rep=5,
                         method="block",
                         seed=321,
                         parallel=TRUE)

wle_men <- unlist(lapply(bootresults,function(x) x[4,1]))
wle_women <- unlist(lapply(bootresults,function(x) x[8,1]))

quantile(wle_men,probs=c(0.025,0.975))
quantile(wle_women,probs=c(0.025,0.975))

## Comparing constrained vs unconstrained model: AIC
fit1 <- dtms_fit(data=estdata,
                 controls=c("Gender","time2"),
                 package="mclogit")

fit2 <- dtms_fullfit(data=estdata,
                     controls=c("Gender","time2"),
                     package="mclogit")

AIC(fit1)
AIC(fit2)

## Comparing constrained vs unconstrained model: results
probs1 <- dtms_transitions(dtms=hrs,
                           model = fit1,
                           constant = list(Gender=0),
                           varying = list(time2 = (50:98)^2),
                           CI=TRUE)
T1 <- dtms_matrix(dtms=hrs,
                  probs=probs1)

probs2 <- dtms_transitions(dtms=hrs,
                           model = fit2,
                           constant = list(Gender=0),
                           varying = list(time2 = (50:98)^2),
                           CI=TRUE)
T2 <- dtms_matrix(dtms=hrs,
                  probs=probs2)

dtms_expectancy(dtms=hrs,
                matrix=T1)

dtms_expectancy(dtms=hrs,
                matrix=T2)
