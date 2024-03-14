library(dtms)
library(devtools)
simple <- dtms(transient=c("A","B"),
               absorbing="X",
               timescale=0:19)
estdata <- dtms_format(data=simpledata,
                       dtms=simple,
                       idvar="id",
                       timevar="time",
                       statevar="state")
estdata <- dtms_clean(data=estdata,
                      dtms=simple)
## Fit model
fit1 <- dtms_fit(data=estdata,package="nnet")
fit2 <- dtms_fit(data=estdata,package="VGAM")
fit3 <- dtms_fit(data=estdata,package="mclogit",random=~1|id)

## Variance-covariance matrices
vcov(fit1)
vcov(fit2)
vcov(fit3)
