## Load package
#library(dtms)
library(devtools)
load_all()
#source("R/dtms_helpers.R")

## Look at data
head(simpledata)

## Number of units
simpledata$id |> unique() |> length()

## Number of observations
dim(simpledata)

## Define model: Absorbing and transient states, time scale
simple <- dtms(transient=c("A","B"),
               absorbing="X",
               timescale=0:19)

# Reshape to transition format
estdata <- dtms_format(data=simpledata,
                       dtms=simple,
                       idvar="id",
                       timevar="time",
                       statevar="state",
                       steplength=TRUE)

## Look at reshaped data
head(estdata)

## Missing values?
estdata$to |> table(useNA="always")


## Clean
estdata <- dtms_clean(data=estdata,
                      dtms=simple)

## Censoring
dtms_censoring(data=estdata,
               dtms=simple)

# More advanced censoring example
censdata <- dtms_censoring(data=estdata,
                          dtms=simple,
                          add=T,
                          addtype="obs")

censdata |>
  subset(subset=to!="X",select=c(RIGHT,to)) |>
  table() |>
  prop.table(margin=1)

## Summary of data
summary(estdata)

## Fit basic model (only starting state as predictor)
fit <- dtms_fit(data=estdata,
                controls="time",
                package="VGAM")

## Fitting model with aggregated data
aggdata <- dtms_aggregate(estdata)
fit_agg <- dtms_fit(data=aggdata,
                 weights="count",
                controls="time",
                package="VGAM")


## Predict probabilities
probs <- dtms_transitions(dtms=simple,
                          controls=list(time=simple$timescale),
                          model = fit)

## Get transition matrix
Tp <- dtms_matrix(dtms=simple,
                  probs=probs)

## Summary of probabilities
summary(probs)

## Simple plot
library(ggplot2)
probs |>  dtms_simplify() |>
  ggplot(aes(x=time,y=P,color=to)) +
  geom_line() +
  facet_wrap(~from)

## Simple base plot
plot(probs,dtms=simple)

## Get starting distribution
S <- dtms_start(dtms=simple,
                data=estdata)

## State expectancies
dtms_expectancy(dtms=simple,
                matrix=Tp,
                start_distr=S)

dtms_expectancy(dtms=simple,
                matrix=Tp,
                start_time=5,
                end_time=7)

dtms_expectancy(dtms=simple,
                risk="A",
                matrix=Tp,
                start_distr=S)

# Time to being absorbed
dtms_absorbed(matrix=Tp,
              dtms=simple,
              start_distr=S)

dtms_absorbed(matrix=Tp,
              dtms=simple,
              start_time=5,
              end_time=5)

summary(dtms_absorbed(matrix=Tp,dtms=simple,start_distr=S)) # Life expectancy

# Survivorship function
dtms_survivor(matrix=Tp,
              dtms=simple,
              start_distr=S)

dtms_survivor(matrix=Tp,
              dtms=simple,
              start_distr=S,
              start_time=5,
              end_time=7)

## Niave lifetime risk
dtms_risk(dtms=simple,
          matrix=Tp,
          risk="A")

dtms_risk(dtms=simple,
          matrix=Tp,
          start_time=5,
          end_time=7,
          start_state="B",
          risk="A")

## Distribution of visits
example <- dtms_visits(dtms=simple,
                       matrix=Tp,
                       risk="A",
                       start_distr=S)
summary(example)

dtms_visits(dtms=simple,
            matrix=Tp,
            risk="A",
            start_time=5,
            end_time=7)

## Distribution of waiting time to last exit
example2 <- dtms_last(dtms=simple,
                      matrix=Tp,
                      risk="A",
                      start_distr=S,
                      rescale=T,
                      total=F)
summary(example2)

dtms_last(dtms=simple,
          matrix=Tp,
          risk="A",
          start_time=5,
          end_time=7)

dtms_first(dtms=simple,
          matrix=Tp,
          risk="A",
          start_time=5,
          end_time=7)


## Check relevance of lags
lags <- dtms_delta(data=estdata,
                   dtms=simple,
                   lags=1:4,
                   controls="time")

## Compare full model with constrained model
fullfit <- dtms_fullfit(data=estdata,
                controls="time",
                package="mclogit")

dtms_delta(data=estdata,
           model1=fit,
           model2=fullfit)

## Correct lifetime risk
riskdata <- dtms_forward(data=estdata,
                         state="A",
                         dtms=simple)

riskfit <- dtms_fit(data=riskdata,
                controls="time",
                package="mclogit")

riskprobs <- dtms_transitions(dtms=simple,
                          controls=list(time=simple$timescale),
                          model = riskfit)

riskTp <- dtms_matrix(dtms=simple,
                  probs=riskprobs)

dtms_risk(dtms=simple,
          matrix=Tp,
          risk="A")

dtms_risk(dtms=simple,
          matrix=riskTp,
          risk="A")

## Duration indicator

durationdata <- dtms_duration(data=simpledata,
                              dtms=simple)

head(durationdata)

## Occurrence indicator
occurrencedata <- dtms_occurrence(data=simpledata,
                                  dtms=simple)

head(occurrencedata)
