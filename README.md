
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dtms - An R package for discrete-time multistate models

<!-- badges: start -->

[![R-CMD-check](https://github.com/r-lib/usethis/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/r-lib/usethis/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Authors

Christian Dudel, <dudel@demogr.mpg.de>

Peng Li, <li@demogr.mpg.de>

## Overview

The package dtms is a user-friendly implementation of discrete-time
multistate models in R. It comes with many tools to analyze the results
of multistate models. In particular, the workflow mainly consists of
estimating a discrete-time multistate model and then applying methods
for absorbing Markov chains. Discrete-time multistate models lend
themselves well for modelling trajectories captured in panel data.

## Disclaimer

This package is currently undergoing development and many functions are
experimental. The content of this repository will change in the future,
and functions and features might be added or removed without warning.

## Installation

You can install the development version of dtms like this:

``` r
library(devtools)
install_github("christiandudel/dtms")
```

## Example 1: Artificial data

This is a basic example using artificial data which is provided with the
package. The state space consists of two transient states (A, B), and
one absorbing state (X). The time scale goes from 0 to 20. Transition
probabilities do change depending on time.

The following code loads the package and the data set. The data set is
called ‘simpledata’.

``` r
## Load package 
library(dtms)

## Look at data
head(simpledata)
#>   id time state
#> 1  1    0     A
#> 2  1    1     B
#> 3  1    2     B
#> 4  1    3     A
#> 5  1    4     B
#> 6  1    5     A

## Number of units
simpledata$id |> unique() |> length()
#> [1] 993

## Number of observations
dim(simpledata)
#> [1] 12173     3
```

The data set is in long format and contains three variables. ‘id’ is an
unit identifier; ‘time’ contains the value of the time scale; and
‘state’ contains the state the unit occupied at a given time. In total,
there are 993 units, each of them contributing to the total of 12,173
observations.

To work with this data set, we first define a basic discrete-time
multistate model. This is done with the function ‘dtms’ and requires us
to specify the names of the transient states, the names of the absorbing
states, and the possible values of the time scale:

``` r
## Define model: Absorbing and transient states, time scale
simple <- dtms(transient=c("A","B"),
               absorbing="X",
               timescale=0:19)
```

The resulting object of class ‘dtms’ can be passed to other functions of
the package, as shown below.

In a second step, we transform the data from long format to what we call
transition format. This can be done using the function ‘dtms_format’. In
this example, we need to specify the name of the object containing the
data, and in addition a ‘dtms’ object as created above. Moreover, we
need to specify which variables contain the unit identifier, which
variable contains the values of the time scale, and which variable
contains the information on the state:

``` r
## Reshape to transition format
estdata <- dtms_format(data=simpledata,
                       dtms=simple,
                       idvar="id",
                       timevar="time",
                       statevar="state")
#> Kept original name for time 
#> Kept original name for id

## Look at reshaped data
head(estdata)
#>   id time from to
#> 1  1    0    A  B
#> 2  1    1    B  B
#> 3  1    2    B  A
#> 4  1    3    A  B
#> 5  1    4    B  A
#> 6  1    5    A  B
```

While in long format each row contains information on the currently
occupied state, in transition format also the next state is shown in a
variable. If there is no observation at time t+1, the next state is NA.
The names of the variables in the resulting data set are by default
chosen such that they match defaults of other functions of the package.

Depending on the original data, there can be missing values in
transition format data due to several reasons. The function ‘dtms_clean’
provides a convenient way to remove such rows in the data, as well as
other potentially problematic or unwanted rows. It returns a cleaned
data set and prints a brief overview of the dropped rows to the console:

``` r
## Missing values?
estdata$to |> table(useNA="always")
#> 
#>    A    B    X <NA> 
#> 3629 6668  616 1260

## Clean
estdata <- dtms_clean(data=estdata,
                      dtms=simple)
#> Dropping  0  rows not in time range
#> Dropping  1260  rows starting or ending in NA
#> Dropping  0  rows starting in absorbing state
```

In this example, 1,260 transitions were dropped because they end in a
missing value. No observations were dropped because they are out of the
time range specified with the ‘dtms’ object, and no observations were
dropped because they start in an absorbing state.

A brief overview of the data is provided when using the function
‘summary’:

``` r
## Summary of data
summary(estdata)
#>   from to COUNT       PROP       PROB
#> 1    A  A   384 0.03518739 0.09741248
#> 2    A  B  3398 0.31137176 0.86199899
#> 3    A  X   160 0.01466141 0.04058853
#> 4    B  A  3245 0.29735178 0.46549993
#> 5    B  B  3270 0.29964263 0.46908621
#> 6    B  X   456 0.04178503 0.06541386
```

This shows for all possible transitions the absolute number each
transition is observed (e.g., there are 160 transitions from A to X);
the proportion of each transition relative to all transitions (e.g., a
bit more than 1% of all observed transitions are from A to X); and raw
transition probabilities (e.g., the probability of transitioning to X
starting in A is around 4%).

Some more information on the data is provided by the function
‘dtms_censoring’. It can be used in different ways, but a basic version
shows an overview of the number of units with left censoring, the number
of units with gaps in their series of observations, and the number of
units with right censoring:

``` r
dtms_censoring(data=estdata,
               dtms=simple)
#> Units with left censoring:  240 
#> Units with gaps:  193 
#> Units with right censoring:  334
```

To estimate the transition probabilities of the multistate model, the
function ‘dtms_fit’ is used. In this simple example, it is enough to
specify the name of the object with the transition data:

``` r
## Fit model 
fit <- dtms_fit(data=estdata)
```

To predict transition probabilities and to arrange them in a matrix, the
functions ‘dtms_transitions’ and ‘dtms_matrix’ are used. The function
‘dtms_transitions’ needs a ‘dtms’ object as well as a fitted model,
while the function ‘dtms_matrix’ requires a ‘dtms’ object and predicted
probabilities:

``` r
## Predict probabilities
probs    <- dtms_transitions(dtms=simple,
                             model = fit)

## Get transition matrix 
Tp <- dtms_matrix(dtms=simple,
                  probs=probs)
```

In more complex examples, the previous functions would need more
information. For instance, on which covariates to include in the
estimation step, and which covariate values to use in the prediction
step. To get an overview of the transition probabilities, the function
summary can be used:

``` r
## Summary of probabilities
summary(probs)
#>   from to    MIN MINtime    MAX MAXtime MEDIAN   MEAN
#> 1    A  A 0.0827      18 0.0994       3 0.0982 0.0956
#> 3    A  B 0.7095      18 0.8950       0 0.8640 0.8420
#> 5    A  X 0.0058       0 0.2078      18 0.0378 0.0624
#> 2    B  A 0.3530      18 0.4884       0 0.4688 0.4512
#> 4    B  B 0.3449      18 0.5018       0 0.4698 0.4527
#> 6    B  X 0.0098       0 0.3021      18 0.0614 0.0961
```

For all combinations of starting and receiving state, this shows the
lowest transition probability and at what value of the time scale it
occurs. It also shows the same for the highest transition probability,
and in addition it shows the median and the mean of all transition
probabilities between two states.

Another useful way to look at the transition probabilies is to plot
them. To make this easy, the package provides the function dtms_simplify
which can be applied to an object created with dtms_transitions to make
it easier to plot. For instance, using ggplot2, a simple plot could look
like this:

``` r
## Simple plot
library(ggplot2)
probs |>  dtms_simplify() |> 
          ggplot(aes(x=time,y=P,color=to)) + 
          geom_line() + 
          facet_wrap(~from)
```

<img src="man/figures/README-example1-probsplot-1.png" width="100%" /> A
simpler way is available which builds on base-R and does not require
ggplot2. However, this creates less nice figures and is mainly intended
as a very quick way of checking results:

``` r
## Simple base plot
plot(probs)
```

<img src="man/figures/README-example1-baseplot-1.png" width="100%" />

Before we generate more results, we calculate the starting distribution
of the states; i.e., the distribution of states at the first value of
the time scale.

``` r
## Get starting distribution 
S <- dtms_start(dtms=simple,
                data=estdata)
```

This step is not necessary, but its result can be passed to several of
the functions used for calculating results, providing additional
information.

Most functions used to calculate results need a transition matrix and a
‘dtms’ object, and potentially further arguments. The two examples below
calculate the expected time spent in a state (dtms_expectancy) and the
lifetime risk of ever reaching a state (dtms_risk). In the first case,
the starting distribution of states is passed to the function; this is
optional. In the second case, one or several states need to be specified
for which the lifetime risk os be calculated:

``` r
## State expectancies 
dtms_expectancy(dtms=simple,
                matrix=Tp,
                start_distr=S)
#>                  A        B    TOTAL
#> start:A_0 5.006818 8.700854 13.70767
#> start:B_0 4.774829 8.896551 13.67138
#> AVERAGE   4.891804 8.797876 13.68968

## Lifetime risk 
dtms_risk(dtms=simple,
          matrix=Tp,
          risk="A")
#>      A_0      B_0 
#> 1.000000 0.974725
```

The results of the call of ‘dtms_expectancy’ show the starting states in
rows and the states in which the time is spent as columns. For instance,
around 5.05 time units are spent in state A when starting in state A at
time 0. The last column shows the total time until absorbtion. The last
row is shown because the starting distribution was specified. It shows
the average time spent in a state irrespective of the starting state.
That is, on average 4.93 time units are spent in state A, 8.88 time
units are spent in state B, for a total of 13.81 time units.

The result of the call of ‘dtms_risk’ above is the lifetime risk of ever
reaching state A depending on the starting state. Obviously, when
starting in state A at time 0, this risk amounts to 1. When starting in
state B, the risk is also very high and around 97%.

The function calls below are all similar in that they provide full
distributions as a result. Specifically, ‘dtms_visits’ calculates the
distribution of the time spent in a state; the mean over this
distribution is equal to the state expectancy as provided by
‘dtms_expectancy’, and the one minus the proportion of 0 time units
spent in a state is equal to the lifetime risk provided by ‘dtms_risk’.
‘dtms_first’ calculates the distribution of the waiting time until a
given state is reached for the first time, conditional on ever reaching
this state. ‘dtms_last’ calculates the distribution of the waiting time
until a state is left for the last time; i.e., there is no return back
to this state.

``` r
## Distribution of visits
dtms_visits(dtms=simple,
            matrix=Tp,
            risk="A",
            start_distr=S)
#>                         0        0.5          1        1.5          2
#> start:A_0      0.00000000 0.03344800 0.00000000 0.05798264 0.00000000
#> start:B_0      0.02527496 0.00000000 0.05016242 0.00000000 0.08200508
#> AVERAGE        0.01253069 0.01686533 0.02486925 0.02923632 0.04065604
#> AVERAGE(COND.) 0.02527496 0.00000000 0.05016242 0.00000000 0.08200508
#>                          2.5             3        3.5          4        4.5
#> start:A_0       9.272739e-02 -2.775558e-17 0.13278564 0.00000000 0.16612892
#> start:B_0      -2.775558e-17  1.210458e-01 0.00000000 0.15724305 0.00000000
#> AVERAGE         4.675550e-02  6.001144e-02 0.06695389 0.07795712 0.08376641
#> AVERAGE(COND.) -2.775558e-17  1.210458e-01 0.00000000 0.15724305 0.00000000
#>                         5        5.5          6        6.5          7
#> start:A_0      0.00000000 0.17786454 0.00000000 0.15780490 0.00000000
#> start:B_0      0.17606542 0.00000000 0.16547248 0.00000000 0.12398646
#> AVERAGE        0.08728877 0.08968381 0.08203706 0.07956923 0.06146934
#> AVERAGE(COND.) 0.17606542 0.00000000 0.16547248 0.00000000 0.12398646
#>                       7.5          8        8.5          9         9.5
#> start:A_0      0.10881912 0.00000000 0.05290557 0.00000000 0.016290535
#> start:B_0      0.00000000 0.06803989 0.00000000 0.02464253 0.000000000
#> AVERAGE        0.05486936 0.03373245 0.02667633 0.01221714 0.008214101
#> AVERAGE(COND.) 0.00000000 0.06803989 0.00000000 0.02464253 0.000000000
#>                         10        10.5           11         11.5           12
#> start:A_0      0.000000000 0.002926048 0.0000000000 0.0002983776 0.000000e+00
#> start:B_0      0.005353188 0.000000000 0.0006604016 0.0000000000 4.634940e-05
#> AVERAGE        0.002653975 0.001475388 0.0003274104 0.0001504495 2.297886e-05
#> AVERAGE(COND.) 0.005353188 0.000000000 0.0006604016 0.0000000000 4.634940e-05
#>                        12.5           13         13.5           14         14.5
#> start:A_0      1.767521e-05 0.000000e+00 6.331159e-07 0.000000e+00 1.418449e-08
#> start:B_0      0.000000e+00 1.915403e-06 0.000000e+00 4.832799e-08 0.000000e+00
#> AVERAGE        8.912289e-06 9.496080e-07 3.192331e-07 2.395979e-08 7.152182e-09
#> AVERAGE(COND.) 0.000000e+00 1.915403e-06 0.000000e+00 4.832799e-08 0.000000e+00
#>                          15         15.5           16         16.5           17
#> start:A_0      0.000000e+00 2.023649e-10 0.000000e+00 1.831979e-12 0.000000e+00
#> start:B_0      7.612809e-10 0.000000e+00 7.491119e-12 0.000000e+00 4.485301e-14
#> AVERAGE        3.774238e-10 1.020375e-10 3.713907e-12 9.237303e-13 2.223699e-14
#> AVERAGE(COND.) 7.612809e-10 0.000000e+00 7.491119e-12 0.000000e+00 4.485301e-14
#>                        17.5 18 18.5 19         19.5            20         20.5
#> start:A_0      1.021405e-14  0    0  0 2.220446e-16 -2.220446e-16 2.220446e-16
#> start:B_0      0.000000e+00  0    0  0 0.000000e+00  0.000000e+00 0.000000e+00
#> AVERAGE        5.150184e-15  0    0  0 1.119605e-16 -1.119605e-16 1.119605e-16
#> AVERAGE(COND.) 0.000000e+00  0    0  0 0.000000e+00  0.000000e+00 0.000000e+00
#> attr(,"class")
#> [1] "dtms_distr" "matrix"

## Distribution of waiting time to first visit
dtms_first(dtms=simple,
           matrix=Tp,
           risk="A",
           start_distr=S)
#>                        0       0.5       1.5        2.5        3.5        4.5
#> start:A_0      1.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> start:B_0      0.0000000 0.5010276 0.2512180 0.12543490 0.06232800 0.03079645
#> AVERAGE        0.5106239 0.2451909 0.1229401 0.06138484 0.03050183 0.01507104
#> AVERAGE(COND.) 0.0000000 0.5010276 0.2512180 0.12543490 0.06232800 0.03079645
#>                        5.5         6.5         7.5          8.5          9.5
#> start:A_0      0.000000000 0.000000000 0.000000000 0.0000000000 0.0000000000
#> start:B_0      0.015116437 0.007362339 0.003552826 0.0016957606 0.0007988565
#> AVERAGE        0.007397623 0.003602953 0.001738668 0.0008298648 0.0003909413
#> AVERAGE(COND.) 0.015116437 0.007362339 0.003552826 0.0016957606 0.0007988565
#>                        10.5         11.5         12.5         13.5         14.5
#> start:A_0      0.0000000000 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> start:B_0      0.0003704930 1.686444e-04 7.506966e-05 3.253792e-05 1.366336e-05
#> AVERAGE        0.0001813104 8.253052e-05 3.673730e-05 1.592328e-05 6.686524e-06
#> AVERAGE(COND.) 0.0003704930 1.686444e-04 7.506966e-05 3.253792e-05 1.366336e-05
#>                        15.5         16.5         17.5 TOTAL(RESCALED)
#> start:A_0      0.000000e+00 0.000000e+00 0.000000e+00               1
#> start:B_0      5.526197e-06 2.138358e-06 7.856347e-07               1
#> AVERAGE        2.704389e-06 1.046462e-06 3.844709e-07               1
#> AVERAGE(COND.) 5.526197e-06 2.138358e-06 7.856347e-07               1
#> attr(,"class")
#> [1] "dtms_distr" "matrix"

## Distribution of waiting time to last exit
dtms_last(dtms=simple,
          matrix=Tp,
          risk="A",
          start_distr=S,
          rescale=T,
          total=F)
#>                       0.5         1.5        2.5        3.5        4.5
#> start:A_0      0.03344862 0.004045972 0.02216475 0.01841978 0.02582569
#> start:B_0      0.00000000 0.020434867 0.01493956 0.02251011 0.02469741
#> AVERAGE        0.01707967 0.012066305 0.01862892 0.02042149 0.02527354
#> AVERAGE(COND.) 0.00000000 0.020434867 0.01493956 0.02251011 0.02469741
#>                       5.5        6.5        7.5        8.5        9.5
#> start:A_0      0.02880081 0.03469428 0.04000662 0.04623384 0.05244868
#> start:B_0      0.03027824 0.03512111 0.04110841 0.04722560 0.05370060
#> AVERAGE        0.02952383 0.03490316 0.04054581 0.04671918 0.05306134
#> AVERAGE(COND.) 0.03027824 0.03512111 0.04110841 0.04722560 0.05370060
#>                      10.5       11.5       12.5       13.5       14.5
#> start:A_0      0.05872871 0.06459506 0.06978833 0.07408259 0.07770617
#> start:B_0      0.06007381 0.06609947 0.07140288 0.07580110 0.07950678
#> AVERAGE        0.05938697 0.06533128 0.07057845 0.07492359 0.07858735
#> AVERAGE(COND.) 0.06007381 0.06609947 0.07140288 0.07580110 0.07950678
#>                      15.5       16.5      17.5       18.5
#> start:A_0      0.08178381 0.08894421 0.1006557 0.07762634
#> start:B_0      0.08367974 0.09100577 0.1029889 0.07942565
#> AVERAGE        0.08271163 0.08995309 0.1017975 0.07850688
#> AVERAGE(COND.) 0.08367974 0.09100577 0.1029889 0.07942565
#> attr(,"class")
#> [1] "dtms_distr" "matrix"
```

The output from these functions tends to be difficult to read, and often
results on the distribution are used to calculate other statistics. A
small set of such statistics can be generated using the function
‘summary’:

``` r
## Distribution of visits
example <- dtms_visits(dtms=simple,
                       matrix=Tp,
                       risk="A",
                       start_distr=S)
summary(example)
#>                    MEAN VARIANCE       SD MEDIAN      RISK0
#> start:A_0      5.006818 4.467193 2.113573    5.5 0.00000000
#> start:B_0      4.774829 4.591303 2.142733    5.0 0.02527496
#> AVERAGE        4.891804 4.542177 2.131238    5.0 0.01253069
#> AVERAGE(COND.) 4.774829 4.591303 2.142733    5.0 0.02527496
```

In the example above this returns, respectively, the average lifetime
spent in state A; the variance of the lifetime spent in state A; the
standard deviation of the lifetime spent in state A; the median of the
lifetime spent in state A; and the probability of spending zero lifetime
in state A. Depending on which distribution this function is applied to,
some entries might not be defined. For instance:

``` r
## Distribution of waiting time to last exit
example2 <- dtms_last(dtms=simple,
                      matrix=Tp,
                      risk="A",
                      start_distr=S,
                      rescale=T,
                      total=F)
summary(example2)
#>                    MEAN VARIANCE       SD MEDIAN RISK0
#> start:A_0      12.12560 23.55401 4.853247   13.5    NA
#> start:B_0      12.40256 20.71913 4.551827   13.5    NA
#> AVERAGE        12.26114 22.18586 4.710186   13.5    NA
#> AVERAGE(COND.) 12.40256 20.71913 4.551827   13.5    NA
```

In this case, the distribution is conditional on ever experiencing the
exit from state A, such that the waiting time until exit always has to
be above 0.

## Example 2: Working trajectories during late working life

Here we provide an example using simulated data based on the Health and
Retirement Study (HRS). The simulations are are conducted using
transition probabilities estimated from the HRS and published by Dudel &
Myrskylä (2017) who studied working trajectories in late working life
and old age. These transition probabilities are used to simulate
artificial but realistic trajectories. There are three transient states
(employed, inactive or unemployed, retired) and one absorbing state
(dead). The time scale represents age and ranges from 50 to 99, as the
focus is on older individuals. Note that the actual HRS data is
collected every two years and while the simulated data is annual. The
data set also contains each individual’s gender, and the transition
probabilities underlying the simulated trajectories differ between men
and women.

The workflow is similar to the previous example. First, a ‘dtms’ model
is defined using the function \`dtms’. Second, the data is brought into
transition format and cleaned. Third, transition probabilities are
estimated and put into a transition matrix. In this example,
probabilities are estimated and predicted using time-constant and
time-varying covariates. Finally, the transition matrix is used to
calculate state expectancies and similar measures.

``` r
## Load package
library(dtms)

## Define model: Absorbing and transient states, time scale
hrs <- dtms(transient=c("Employed","Inactive","Retired"),
            absorbing="Dead",
            timescale=50:99)

## Quick look at data
head(hrsdata)
#>   ID Gender Age State
#> 1  1      1  50  <NA>
#> 2  1      1  51  <NA>
#> 3  1      1  52  <NA>
#> 4  1      1  53  <NA>
#> 5  1      1  54  <NA>
#> 6  1      1  55  <NA>

## Reshape
estdata <- dtms_format(data=hrsdata,
                       dtms=hrs,
                       idvar="ID",
                       timevar="Age",
                       statevar="State")

## Drop dead-to-dead transitions etc
estdata <- dtms_clean(data=estdata,
                      dtms=hrs)
#> Dropping  0  rows not in time range
#> Dropping  98287  rows starting or ending in NA
#> Dropping  51935  rows starting in absorbing state

## Overview
summary(estdata)
#>        from       to COUNT        PROP        PROB
#> 1  Employed     Dead   306 0.003066808 0.008699855
#> 2  Employed Employed 30623 0.306911343 0.870639411
#> 3  Employed Inactive  2066 0.020705967 0.058738237
#> 4  Employed  Retired  2178 0.021828459 0.061922497
#> 5  Inactive     Dead   197 0.001974383 0.013936050
#> 6  Inactive Employed  1404 0.014071238 0.099320883
#> 7  Inactive Inactive 10635 0.106586622 0.752334465
#> 8  Inactive  Retired  1900 0.019042274 0.134408602
#> 9   Retired     Dead  2602 0.026077893 0.051556401
#> 10  Retired Employed   838 0.008398645 0.016604252
#> 11  Retired Inactive   606 0.006073483 0.012007371
#> 12  Retired  Retired 46423 0.465262884 0.919831976
dtms_censoring(data=estdata,
               dtms=hrs)
#> Units with left censoring:  2036 
#> Units with gaps:  1720 
#> Units with right censoring:  1323

## Add age squared
estdata$time2 <- estdata$time^2
  
## Fit model
fit <- dtms_fit(data=estdata,
                controls=c("Gender","time2"))

## Transition probabilities by gender
  
# Men
probs_m <- dtms_transitions(dtms=hrs,
                             model = fit,
                             constant = list(Gender=0),
                             varying = list(time2 = (50:98)^2))
  
# Women
probs_w <- dtms_transitions(dtms=hrs,
                            model = fit,
                            constant = list(Gender=1),
                            varying = list(time2 = (50:98)^2))

# Overview
summary(probs_m)
#>        from       to    MIN MINtime    MAX MAXtime MEDIAN   MEAN
#> 1  Employed     Dead 0.0013      50 0.4725      98 0.0288 0.0926
#> 4  Employed Employed 0.3990      98 0.9456      50 0.7862 0.7447
#> 7  Employed Inactive 0.0000      98 0.0704      58 0.0104 0.0252
#> 10 Employed  Retired 0.0075      50 0.2443      85 0.1524 0.1375
#> 2  Inactive     Dead 0.0064      50 0.7701      98 0.1080 0.2055
#> 5  Inactive Employed 0.0238      98 0.1668      50 0.1008 0.0879
#> 8  Inactive Inactive 0.0000      98 0.8173      54 0.1402 0.3186
#> 11 Inactive  Retired 0.0356      50 0.7089      79 0.3904 0.3880
#> 3   Retired     Dead 0.0231      58 0.4588      98 0.0364 0.0930
#> 6   Retired Employed 0.0047      98 0.2132      50 0.0121 0.0365
#> 9   Retired Inactive 0.0000      98 0.1627      50 0.0025 0.0324
#> 12  Retired  Retired 0.5365      98 0.9492      73 0.8855 0.8381
summary(probs_w)
#>        from       to    MIN MINtime    MAX MAXtime MEDIAN   MEAN
#> 1  Employed     Dead 0.0009      50 0.3908      98 0.0208 0.0718
#> 4  Employed Employed 0.4542      98 0.9311      50 0.7815 0.7462
#> 7  Employed Inactive 0.0000      98 0.0922      58 0.0138 0.0332
#> 10 Employed  Retired 0.0078      50 0.2633      86 0.1700 0.1488
#> 2  Inactive     Dead 0.0036      50 0.6979      98 0.0743 0.1649
#> 5  Inactive Employed 0.0297      98 0.1312      50 0.0822 0.0791
#> 8  Inactive Inactive 0.0000      98 0.8565      54 0.1780 0.3471
#> 11 Inactive  Retired 0.0296      50 0.7504      80 0.4309 0.4090
#> 3   Retired     Dead 0.0156      57 0.3678      98 0.0253 0.0687
#> 6   Retired Employed 0.0052      98 0.1966      50 0.0115 0.0341
#> 9   Retired Inactive 0.0000      98 0.2013      50 0.0032 0.0404
#> 12  Retired  Retired 0.5857      50 0.9600      74 0.9001 0.8567

# Plotting, men as example
probs_m |>  dtms_simplify() |> 
            ggplot(aes(x=time,y=P,color=to)) + 
            geom_line() + 
            facet_wrap(~from)
```

<img src="man/figures/README-example2-1.png" width="100%" />

``` r
 
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
#>                    Employed Inactive  Retired    TOTAL
#> start:Employed_50 13.307334 3.074605 13.53758 29.91952
#> start:Inactive_50  8.782989 6.496935 13.75116 29.03108
#> start:Retired_50   8.821763 3.905134 15.19067 27.91757
#> AVERAGE           12.445981 3.589228 13.65526 29.69047
    
dtms_expectancy(dtms=hrs,
                matrix=Tw,
                start_distr=Sw)
#>                    Employed Inactive  Retired    TOTAL
#> start:Employed_50 12.066485 4.386659 16.53253 32.98567
#> start:Inactive_50  7.445512 8.199122 16.77896 32.42359
#> start:Retired_50   7.836675 5.486704 18.24103 31.56440
#> AVERAGE           10.450573 5.607548 16.68838 32.74651

## A variant: ignoring retirement as a starting state (shown only for men)
limited <- c("Employed","Inactive")

Smwr <- dtms_start(dtms=hrs,
                   data=estdata,
                   start_state=limited,
                   variables=list(Gender=0))

dtms_expectancy(dtms=hrs,
                matrix=Tm,
                start_state=limited,
                start_distr=Smwr)
#>                    Employed Inactive  Retired    TOTAL
#> start:Employed_50 13.307334 3.074605 13.53758 29.91952
#> start:Inactive_50  8.782989 6.496935 13.75116 29.03108
#> AVERAGE           12.650574 3.571394 13.56859 29.79056

## Lifetime risk of reaching retirement
dtms_risk(dtms=hrs,
          matrix=Tm,
          risk="Retired",
          start_distr=Sm)
#>    Employed_50    Inactive_50     Retired_50        AVERAGE AVERAGE(COND.) 
#>      0.8828701      0.8806115      1.0000000      0.8888186      0.8825422
  
dtms_risk(dtms=hrs,
          matrix=Tw,
          risk="Retired",
          start_distr=Sw)
#>    Employed_50    Inactive_50     Retired_50        AVERAGE AVERAGE(COND.) 
#>      0.9172407      0.9159547      1.0000000      0.9207352      0.9168268
  
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
#>                       MEAN  VARIANCE        SD MEDIAN     RISK0
#> start:Employed_50 13.53758  97.35708  9.866969   13.0 0.1171299
#> start:Inactive_50 13.75116 103.70848 10.183736   13.0 0.1193885
#> start:Retired_50  15.19067 112.59117 10.610899   14.5 0.0000000
#> AVERAGE           13.65526  99.18227  9.959030   13.0 0.1111814
#> AVERAGE(COND.)    13.56859  98.28472  9.913865   13.0 0.1174578
summary(visitsw)
#>                       MEAN VARIANCE       SD MEDIAN      RISK0
#> start:Employed_50 17.44977 112.1093 10.58817     18 0.08275934
#> start:Inactive_50 17.69491 117.8990 10.85813     18 0.08404534
#> start:Retired_50  18.74103 124.6118 11.16296     19 0.00000000
#> AVERAGE           17.58562 114.5507 10.70284     18 0.07926479
#> AVERAGE(COND.)    17.52865 113.9856 10.67640     18 0.08317318
  
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
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                       MEAN VARIANCE       SD MEDIAN      RISK0
#> start:Employed_50 14.25887 42.52401 6.521043   14.5 0.00000000
#> start:Inactive_50 12.31587 50.90513 7.134783   12.5 0.00000000
#> start:Retired_50   0.00000  0.00000 0.000000    0.0 1.00000000
#> AVERAGE           13.13712 52.58727 7.251708   13.5 0.06011927
#> AVERAGE(COND.)    13.97744 44.20558 6.648728   13.5 0.00000000
summary(firstw)
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                       MEAN VARIANCE       SD MEDIAN      RISK0
#> start:Employed_50 14.10717 40.00302 6.324794   14.5 0.00000000
#> start:Inactive_50 12.54709 46.49749 6.818907   12.5 0.00000000
#> start:Retired_50   0.00000  0.00000 0.000000    0.0 1.00000000
#> AVERAGE           12.91123 49.41190 7.029360   13.5 0.05103632
#> AVERAGE(COND.)    13.60562 42.62186 6.528542   13.5 0.00000000

## Last exit
  
# Leaving employment to any state
last1m <- dtms_last(dtms=hrs,
                    matrix=Tm,
                    risk="Employed",
                    start_distr=Sm)  
  
last1w <- dtms_last(dtms=hrs,
                    matrix=Tw,
                    risk="Employed",
                    start_distr=Sw) 

summary(last1m)
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                       MEAN VARIANCE       SD MEDIAN RISK0
#> start:Employed_50 16.50265 76.98676 8.774210   15.5    NA
#> start:Inactive_50 18.02302 68.14259 8.254853   17.5    NA
#> start:Retired_50  17.83146 69.33252 8.326615   17.5    NA
#> AVERAGE           16.73797 75.91238 8.712771   16.5    NA
#> AVERAGE(COND.)    17.97027 68.47754 8.275116   17.5    NA
summary(last1w)
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                       MEAN VARIANCE       SD MEDIAN RISK0
#> start:Employed_50 16.15218 87.76963 9.368545   15.5    NA
#> start:Inactive_50 18.31738 77.06422 8.778623   17.5    NA
#> start:Retired_50  17.94411 78.93510 8.884543   17.5    NA
#> AVERAGE           16.78741 85.57486 9.250668   16.5    NA
#> AVERAGE(COND.)    18.26738 77.33095 8.793802   17.5    NA
  
# Leaving employment for retirement
last2m <- dtms_last(dtms=hrs,
                    matrix=Tm,
                    risk="Employed",
                    risk_to="Retired",
                    start_distr=Sm)  
  
last2w <- dtms_last(dtms=hrs,
                    matrix=Tw,
                    risk="Employed",
                    risk_to="Retired",
                    start_distr=Sw)  

summary(last2m)
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                       MEAN VARIANCE       SD MEDIAN RISK0
#> start:Employed_50 18.74988 64.64429 8.040167   18.5    NA
#> start:Inactive_50 19.72542 56.78054 7.535286   19.5    NA
#> start:Retired_50  19.60618 57.79592 7.602363   19.5    NA
#> AVERAGE           18.90783 63.49802 7.968565   18.5    NA
#> AVERAGE(COND.)    19.69276 57.06149 7.553906   19.5    NA
summary(last2w)
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                       MEAN VARIANCE       SD MEDIAN RISK0
#> start:Employed_50 19.33660 73.65218 8.582085   19.5    NA
#> start:Inactive_50 20.62023 63.26114 7.953687   20.5    NA
#> start:Retired_50  20.40228 64.97334 8.060604   20.5    NA
#> AVERAGE           19.73733 70.75006 8.411305   19.5    NA
#> AVERAGE(COND.)    20.59143 63.49278 7.968235   20.5    NA
```
