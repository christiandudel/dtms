## Load package
library(devtools)
#load_all()
library(dtms)
library(ggplot2)
source("R/dtms_helpers.R")

## Define model: Absorbing and transient states, time scale
work <- dtms(transient=c("Working","Non-working","Retired"),
            absorbing="Dead",
            timescale=50:99)

## Quick look at data
head(workdata)

## Reshape
estdata <- dtms_format(data=workdata,
                       dtms=work,
                       idvar="ID",
                       timevar="Age",
                       statevar="State")

## Drop dead-to-dead transitions etc
estdata <- dtms_clean(data=estdata,
                      dtms=work)

## Overview
summary(estdata)

## Basic censoring
dtms_censoring(data=estdata,
               dtms=work)

## More advanced censoring example
estdata <- dtms_censoring(data=estdata,
                          dtms=work,
                          add=T,
                          addtype="obs")

estdata |>
  subset(subset=to!="Dead",select=c(RIGHT,to)) |>
  table() |>
  prop.table(margin=1)

## Add age squared (used below)
estdata$time2 <- estdata$time^2

## Fit model with spline
fit <- dtms_fit(data=estdata,
                controls=c("Gender","s(time)"),
                package="VGAM")

## Transition probabilities by gender

# Men
probs_m <- dtms_transitions(dtms=work,
                            model = fit,
                            controls = list(Gender=0,
                                            time  =50:98,
                                            time2 =(50:98)^2),
                            CI=TRUE)

# Women
probs_w <- dtms_transitions(dtms=work,
                            model = fit,
                            controls = list(Gender=1,
                                            time  =50:98,
                                            time2 =(50:98)^2))

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
Tm <- dtms_matrix(dtms=work,
                  probs=probs_m)

Tw <- dtms_matrix(dtms=work,
                  probs=probs_w)

## Starting distributions
Sm <- dtms_start(dtms=work,
                 data=estdata,
                 variables=list(Gender=0))

Sw <- dtms_start(dtms=work,
                 data=estdata,
                 variables=list(Gender=1))

## State expectancies
dtms_expectancy(dtms=work,
                matrix=Tm,
                start_distr=Sm)

dtms_expectancy(dtms=work,
                matrix=Tw,
                start_distr=Sw)

## A variant: ignoring retirement as a starting state (shown only for men)
limited <- c("Working","Non-working")

Smwr <- dtms_start(dtms=work,
                   data=estdata,
                   start_state=limited,
                   variables=list(Gender=0))

dtms_expectancy(dtms=work,
                matrix=Tm,
                start_state=limited,
                start_distr=Smwr)

## Lifetime risk of reaching retirement
dtms_risk(dtms=work,
          matrix=Tm,
          risk="Retired",
          start_distr=Sm)

dtms_risk(dtms=work,
          matrix=Tw,
          risk="Retired",
          start_distr=Sw)

## Distribution of visits
visitsm <- dtms_visits(dtms=work,
                       matrix=Tm,
                       risk="Retired",
                       start_distr=Sm)

visitsw <- dtms_visits(dtms=work,
                       matrix=Tw,
                       risk="Retired",
                       start_distr=Sw,
                       method="end")

summary(visitsm)
summary(visitsw)

## First visit
firstm <- dtms_first(dtms=work,
                     matrix=Tm,
                     risk="Retired",
                     start_distr=Sm)

firstw <- dtms_first(dtms=work,
                     matrix=Tw,
                     risk="Retired",
                     start_distr=Sw)

summary(firstm)
summary(firstw)

## Last exit

# Leaving employment to any state
last1m <- dtms_last(dtms=work,
                    matrix=Tm,
                    risk="Working",
                    start_distr=Sm)

last1w <- dtms_last(dtms=work,
                    matrix=Tw,
                    risk="Working",
                    start_distr=Sw)

summary(last1m)
summary(last1w)

# Leaving employment for retirement
last2m <- dtms_last(dtms=work,
                    matrix=Tm,
                    risk="Working",
                    risk_to="Retired",
                    start_distr=Sm)

last2w <- dtms_last(dtms=work,
                    matrix=Tw,
                    risk="Working",
                    risk_to="Retired",
                    start_distr=Sw)

summary(last2m)
summary(last2w)

# Bootstrap example
bootfun <- function(data,dtms) {

  fit <- dtms_fit(data=data,
                  controls=c("Gender","time","time2"),
                  package="mclogit")

  probs_m <- dtms_transitions(dtms=dtms,
                              model = fit,
                              controls = list(Gender=0,
                                              time  =50:98,
                                              time2 =(50:98)^2))

  probs_w <- dtms_transitions(dtms=dtms,
                              model = fit,
                              controls = list(Gender=1,
                                              time  =50:98,
                                              time2 =(50:98)^2))

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
                         dtms=work,
                         fun=bootfun,
                         idvar="id",
                         rep=5,
                         method="block",
                         parallel=TRUE)

summary(bootresults)

wle_men <- unlist(lapply(bootresults,function(x) x[4,1]))
wle_women <- unlist(lapply(bootresults,function(x) x[8,1]))

quantile(wle_men,probs=c(0.025,0.975))
quantile(wle_women,probs=c(0.025,0.975))

## Comparing constrained vs unconstrained model: likelihood ratio test
## Restricting to men, because we know the true model
estdatam <- subset(estdata,Gender==0)

fit1 <- dtms_fit(data=estdatam,
                 #controls=c("Gender","time2"),
                 controls=c("time","time2"),
                 package="mclogit")

fit2 <- dtms_fullfit(data=estdatam,
                     #controls=c("Gender","time2"),
                     controls=c("time","time2"),
                     package="mclogit")

llfit1 <- logLik(fit1)
llfit2 <- logLik(fit2)

lldiff <- -2 * (llfit1[1]-llfit2[1])
dftest <- attr(llfit2,"df")-attr(llfit1,"df")

pchisq(lldiff, df = dftest, lower.tail = FALSE)


## Comparing constrained vs unconstrained model: results
probs1 <- dtms_transitions(dtms=work,
                           model = fit1,
                           controls = list(Gender=0,
                                           time  =50:98,
                                           time2 =(50:98)^2))
T1 <- dtms_matrix(dtms=work,
                  probs=probs1)

probs2 <- dtms_transitions(dtms=work,
                           model = fit2,
                           controls = list(Gender=0,
                                            time  =50:98,
                                            time2 =(50:98)^2))
T2 <- dtms_matrix(dtms=work,
                  probs=probs2)

dtms_expectancy(dtms=work,
                matrix=T1)

dtms_expectancy(dtms=work,
                matrix=T2)

## Correct lifetime risk
riskdata <- dtms_forward(data=estdata,
                         state="Retired",
                         dtms=work)

riskfit <- dtms_fit(data=riskdata,
                controls=c("Gender","time","time2"),
                package="mclogit")

riskprobs <- dtms_transitions(dtms=work,
                            model = riskfit,
                            controls = list(Gender=0,
                                            time  =50:98,
                                            time2 =(50:98)^2),
                            CI=TRUE)


riskTp <- dtms_matrix(dtms=work,
                      probs=riskprobs)

dtms_risk(dtms=work,
          matrix=riskTp,
          risk="Retired")
