
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
               timescale=0:20)
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
# Fit model 
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

To estimate the transition probabilities of the multistate model, the
function ‘dtms_fit’ is used. In this simple example, it is enough to
specify the name of the object with the transition data:

``` r
# Fit model 
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
#> 1    A  A 0.0791      19 0.0994       3 0.0979 0.0948
#> 3    A  B 0.6763      19 0.8950       0 0.8601 0.8337
#> 5    A  X 0.0058       0 0.2447      19 0.0420 0.0715
#> 2    B  A 0.3302      19 0.4884       0 0.4658 0.4451
#> 4    B  B 0.3217      19 0.5018       0 0.4661 0.4462
#> 6    B  X 0.0098       0 0.3481      19 0.0681 0.1087
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
plot(probs)
```

<img src="man/figures/README-example1-baseplot-1.png" width="100%" />

Before we generate several results, we calculate the starting
distribution of the states; i.e., the distribution of states at the
first value of the time scale.

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
#> start:A_0 5.049179 8.779119 13.82830
#> start:B_0 4.817076 8.974606 13.79168
#> AVERAGE   4.934108 8.876037 13.81014

## Lifetime risk 
dtms_risk(dtms=simple,
          matrix=Tp,
          risk="A")
#>       A_0       B_0 
#> 1.0000000 0.9747251
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
#>                         0           0.5          1        1.5          2
#> start:A_0      0.00000000  3.344785e-02 0.00000000 0.05797786 0.00000000
#> start:B_0      0.02527488 -3.469447e-18 0.05015961 0.00000000 0.08196620
#> AVERAGE        0.01253064  1.686525e-02 0.02486786 0.02923391 0.04063676
#> AVERAGE(COND.) 0.02527488 -3.469447e-18 0.05015961 0.00000000 0.08196620
#>                       2.5          3        3.5          4        4.5
#> start:A_0      0.09266499 0.00000000 0.13234577 0.00000000 0.16433225
#> start:B_0      0.00000000 0.12075167 0.00000000 0.15592903 0.00000000
#> AVERAGE        0.04672404 0.05986562 0.06673209 0.07730566 0.08286049
#> AVERAGE(COND.) 0.00000000 0.12075167 0.00000000 0.15592903 0.00000000
#>                         5        5.5          6        6.5          7
#> start:A_0      0.00000000 0.17367290 0.00000000 0.15298371 0.00000000
#> start:B_0      0.17259920 0.00000000 0.16056822 0.00000000 0.12189236
#> AVERAGE        0.08557031 0.08757028 0.07960565 0.07713827 0.06043114
#> AVERAGE(COND.) 0.17259920 0.00000000 0.16056822 0.00000000 0.12189236
#>                       7.5          8        8.5          9        9.5
#> start:A_0      0.10829957 0.00000000 0.05765539 0.00000000 0.02107182
#> start:B_0      0.00000000 0.07149131 0.00000000 0.02983741 0.00000000
#> AVERAGE        0.05460739 0.03544358 0.02907131 0.01479263 0.01062494
#> AVERAGE(COND.) 0.00000000 0.07149131 0.00000000 0.02983741 0.00000000
#>                         10        10.5           11         11.5           12
#> start:A_0      0.000000000 0.004833560 0.0000000000 0.0006585990 0.000000e+00
#> start:B_0      0.008078237 0.000000000 0.0013183451 0.0000000000 1.260989e-04
#> AVERAGE        0.004004985 0.002437203 0.0006536021 0.0003320823 6.251665e-05
#> AVERAGE(COND.) 0.008078237 0.000000000 0.0013183451 0.0000000000 1.260989e-04
#>                        12.5           13         13.5           14         14.5
#> start:A_0      5.305246e-05 0.000000e+00 2.597229e-06 0.000000e+00 7.993303e-08
#> start:B_0      0.000000e+00 7.176230e-06 0.000000e+00 2.508942e-07 0.000000e+00
#> AVERAGE        2.675040e-05 3.557793e-06 1.309589e-06 1.243870e-07 4.030426e-08
#> AVERAGE(COND.) 0.000000e+00 7.176230e-06 0.000000e+00 2.508942e-07 0.000000e+00
#>                          15         15.5           16         16.5           17
#> start:A_0      0.000000e+00 1.586369e-09 0.000000e+00 2.052092e-11 2.220446e-16
#> start:B_0      5.542817e-09 0.000000e+00 7.849632e-11 2.220446e-16 7.083223e-13
#> AVERAGE        2.747988e-09 7.998875e-10 3.891649e-11 1.034728e-11 3.512802e-13
#> AVERAGE(COND.) 5.542817e-09 0.000000e+00 7.849632e-11 2.220446e-16 7.083223e-13
#>                        17.5           18         18.5           19
#> start:A_0      1.711964e-13 0.000000e+00 6.661338e-16 2.220446e-16
#> start:B_0      2.220446e-16 3.552714e-15 2.220446e-16 2.220446e-16
#> AVERAGE        8.643164e-14 1.761345e-15 4.459656e-16 2.220446e-16
#> AVERAGE(COND.) 2.220446e-16 3.552714e-15 2.220446e-16 2.220446e-16
#>                         19.5            20 20.5           21          21.5
#> start:A_0       2.220446e-16 -2.220446e-16    0 2.220446e-16  0.000000e+00
#> start:B_0      -2.220446e-16  2.220446e-16    0 0.000000e+00 -2.220446e-16
#> AVERAGE         1.876433e-18 -1.876433e-18    0 1.119605e-16 -1.100841e-16
#> AVERAGE(COND.) -2.220446e-16  2.220446e-16    0 0.000000e+00 -2.220446e-16
#> attr(,"class")
#> [1] "dtms_distr" "matrix"

## Distribution of waiting time to first visit
dtms_first(dtms=simple,
           matrix=Tp,
           risk="A",
           start_distr=S)
#>                        0       0.5       1.5        2.5        3.5        4.5
#> start:A_0      1.0000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> start:B_0      0.0000000 0.5010275 0.2512179 0.12543486 0.06232798 0.03079644
#> AVERAGE        0.5106238 0.2451909 0.1229401 0.06138483 0.03050183 0.01507104
#> AVERAGE(COND.) 0.0000000 0.5010275 0.2512179 0.12543486 0.06232798 0.03079644
#>                        5.5         6.5         7.5          8.5          9.5
#> start:A_0      0.000000000 0.000000000 0.000000000 0.0000000000 0.0000000000
#> start:B_0      0.015116433 0.007362337 0.003552825 0.0016957602 0.0007988563
#> AVERAGE        0.007397622 0.003602952 0.001738668 0.0008298646 0.0003909412
#> AVERAGE(COND.) 0.015116433 0.007362337 0.003552825 0.0016957602 0.0007988563
#>                        10.5         11.5         12.5         13.5         14.5
#> start:A_0      0.0000000000 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> start:B_0      0.0003704929 1.686443e-04 7.506964e-05 3.253791e-05 1.366336e-05
#> AVERAGE        0.0001813104 8.253051e-05 3.673729e-05 1.592328e-05 6.686523e-06
#> AVERAGE(COND.) 0.0003704929 1.686443e-04 7.506964e-05 3.253791e-05 1.366336e-05
#>                        15.5         16.5         17.5         18.5
#> start:A_0      0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> start:B_0      5.526195e-06 2.138358e-06 7.856345e-07 2.717496e-07
#> AVERAGE        2.704388e-06 1.046461e-06 3.844708e-07 1.329878e-07
#> AVERAGE(COND.) 5.526195e-06 2.138358e-06 7.856345e-07 2.717496e-07
#>                TOTAL(RESCALED)
#> start:A_0                    1
#> start:B_0                    1
#> AVERAGE                      1
#> AVERAGE(COND.)               1
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
#> start:A_0      0.03344805 0.004045859 0.02216373 0.01841838 0.02582240
#> start:B_0      0.00000000 0.020434290 0.01493887 0.02250838 0.02469425
#> AVERAGE        0.01707937 0.012065966 0.01862806 0.02041993 0.02527031
#> AVERAGE(COND.) 0.00000000 0.020434290 0.01493887 0.02250838 0.02469425
#>                       5.5        6.5        7.5        8.5        9.5
#> start:A_0      0.02879462 0.03468157 0.03998133 0.04618269 0.05234543
#> start:B_0      0.03027172 0.03510823 0.04108241 0.04717334 0.05359487
#> AVERAGE        0.02951748 0.03489037 0.04052017 0.04666749 0.05295687
#> AVERAGE(COND.) 0.03027172 0.03510823 0.04108241 0.04717334 0.05359487
#>                      10.5       11.5       12.5       13.5       14.5
#> start:A_0      0.05851903 0.06416783 0.06891338 0.07228196 0.07398394
#> start:B_0      0.05985931 0.06566226 0.07050766 0.07395867 0.07569828
#> AVERAGE        0.05917493 0.06489917 0.06969359 0.07310250 0.07482290
#> AVERAGE(COND.) 0.05985931 0.06566226 0.07050766 0.07395867 0.07569828
#>                      15.5       16.5       17.5       18.5       19.5
#> start:A_0      0.07409558 0.07348171 0.07420418 0.07762634 0.05684198
#> start:B_0      0.07581325 0.07518486 0.07592420 0.07942562 0.05815952
#> AVERAGE        0.07493617 0.07431519 0.07504592 0.07850686 0.05748675
#> AVERAGE(COND.) 0.07581325 0.07518486 0.07592420 0.07942562 0.05815952
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
#> start:A_0      5.049179 4.677853 2.162835    5.5 0.00000000
#> start:B_0      4.817076 4.802630 2.191490    5.0 0.02527488
#> AVERAGE        4.934108 4.753181 2.180179    5.0 0.01253064
#> AVERAGE(COND.) 4.817076 4.802630 2.191490    5.0 0.02527488
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
#> start:A_0      12.29874 25.43944 5.043752   13.5    NA
#> start:B_0      12.57971 22.54946 4.748628   13.5    NA
#> AVERAGE        12.43624 24.04488 4.903558   13.5    NA
#> AVERAGE(COND.) 12.57971 22.54946 4.748628   13.5    NA
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
