
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
step.

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
#>                         0          0.5            1          1.5            2
#> start:A_0      0.00000000 3.344785e-02 1.524098e-07 5.797771e-02 4.926242e-06
#> start:B_0      0.02527488 8.546268e-08 5.015952e-02 2.894042e-06 8.196331e-02
#> AVERAGE        0.01253064 1.686530e-02 2.486790e-02 2.923527e-02 4.063781e-02
#> AVERAGE(COND.) 0.02527488 8.546268e-08 5.015952e-02 2.894042e-06 8.196331e-02
#>                         2.5            3         3.5            4         4.5
#> start:A_0      9.266006e-02 6.732588e-05 0.132278444 0.0005072002 0.163825048
#> start:B_0      4.177893e-05 1.207099e-01 0.000335906 0.1555931217 0.001649927
#> AVERAGE        4.674227e-02 5.987885e-02 0.066864678 0.0773948683 0.083422734
#> AVERAGE(COND.) 4.177893e-05 1.207099e-01 0.000335906 0.1555931217 0.001649927
#>                          5         5.5           6        6.5          7
#> start:A_0      0.002303868 0.171369034 0.006495506 0.14648821 0.01131669
#> start:B_0      0.170949275 0.005116145 0.155452071 0.01002041 0.11187195
#> AVERAGE        0.085913985 0.088945066 0.080344395 0.07883093 0.06116944
#> AVERAGE(COND.) 0.170949275 0.005116145 0.155452071 0.01002041 0.11187195
#>                       7.5          8         8.5           9         9.5
#> start:A_0      0.09698287 0.01183625 0.045819144 0.007086424 0.013985394
#> start:B_0      0.01211451 0.05937680 0.008663089 0.021174322 0.003468210
#> AVERAGE        0.05490729 0.03540565 0.027398114 0.014070847 0.008771241
#> AVERAGE(COND.) 0.01211451 0.05937680 0.008663089 0.021174322 0.003468210
#>                         10         10.5           11         11.5           12
#> start:A_0      0.002305142 0.0025284183 0.0003976300 2.609691e-04 3.740852e-05
#> start:B_0      0.004610027 0.0007431612 0.0005751838 8.521778e-05 4.088115e-05
#> AVERAGE        0.003447845 0.0016433331 0.0004856567 1.738360e-04 3.913016e-05
#> AVERAGE(COND.) 0.004610027 0.0007431612 0.0005751838 8.521778e-05 4.088115e-05
#>                        12.5           13         13.5           14         14.5
#> start:A_0      1.564394e-05 2.031264e-06 5.659645e-07 6.715140e-08 1.278164e-08
#> start:B_0      5.468247e-06 1.707983e-06 2.074194e-07 4.347478e-08 4.853210e-09
#> AVERAGE        1.059909e-05 1.870990e-06 3.882069e-07 5.541313e-08 8.850925e-09
#> AVERAGE(COND.) 5.468247e-06 1.707983e-06 2.074194e-07 4.347478e-08 4.853210e-09
#>                          15         15.5           16         16.5           17
#> start:A_0      1.402855e-09 1.835138e-10 1.885092e-11 1.669886e-12 1.622036e-13
#> start:B_0      6.896080e-10 7.167300e-11 6.823320e-12 6.676881e-13 4.085621e-14
#> AVERAGE        1.049245e-09 1.280659e-10 1.288794e-11 1.173022e-12 1.020426e-13
#> AVERAGE(COND.) 6.896080e-10 7.167300e-11 6.823320e-12 6.676881e-13 4.085621e-14
#>                        17.5           18         18.5 19          19.5
#> start:A_0      9.325873e-15 7.771561e-16 2.220446e-16  0 -2.220446e-16
#> start:B_0      3.774758e-15 2.220446e-16 0.000000e+00  0  0.000000e+00
#> AVERAGE        6.573771e-15 5.019459e-16 1.119605e-16  0 -1.119605e-16
#> AVERAGE(COND.) 3.774758e-15 2.220446e-16 0.000000e+00  0  0.000000e+00
#>                          20 20.5
#> start:A_0      2.220446e-16    0
#> start:B_0      0.000000e+00    0
#> AVERAGE        1.119605e-16    0
#> AVERAGE(COND.) 0.000000e+00    0

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
```

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
dtms_visits(dtms=hrs,
            matrix=Tm,
            risk="Retired",
            start_distr=Sm)
#>                           0          0.5            1          1.5            2
#> start:Employed_50 0.1171299 9.360193e-08 3.338830e-02 1.507503e-07 3.165376e-02
#> start:Inactive_50 0.1193885 7.303538e-08 3.551620e-02 1.215724e-07 3.265360e-02
#> start:Retired_50  0.0000000 6.864444e-02 3.210211e-08 5.126978e-02 6.982272e-08
#> AVERAGE           0.1111814 3.668109e-03 3.189657e-02 2.739745e-03 3.009973e-02
#> AVERAGE(COND.)    0.1174578 9.061647e-08 3.369719e-02 1.465148e-07 3.179890e-02
#>                            2.5            3          3.5            4
#> start:Employed_50 2.366803e-07 3.056120e-02 3.662404e-07 2.990759e-02
#> start:Inactive_50 1.936576e-07 3.085437e-02 3.018077e-07 2.972953e-02
#> start:Retired_50  4.169992e-02 1.211402e-07 3.612819e-02 1.955214e-07
#> AVERAGE           2.228458e-03 2.896845e-02 1.930852e-03 2.828502e-02
#> AVERAGE(COND.)    2.304351e-07 3.060376e-02 3.568873e-07 2.988174e-02
#>                            4.5            5          5.5            6
#> start:Employed_50 5.609928e-07 2.956740e-02 8.518947e-07 2.945739e-02
#> start:Inactive_50 4.640912e-07 2.905712e-02 7.064232e-07 2.870169e-02
#> start:Retired_50  3.274916e-02 3.057844e-07 3.065171e-02 4.700002e-07
#> AVERAGE           1.750473e-03 2.791737e-02 1.638664e-03 2.777952e-02
#> AVERAGE(COND.)    5.469264e-07 2.949333e-02 8.307779e-07 2.934769e-02
#>                            6.5            7          7.5            8
#> start:Employed_50 1.282735e-06 2.951813e-02 1.914391e-06 2.970389e-02
#> start:Inactive_50 1.065493e-06 2.857463e-02 1.592399e-06 2.861362e-02
#> start:Retired_50  2.934796e-02 7.137462e-07 2.856100e-02 1.072877e-06
#> AVERAGE           1.569396e-03 2.781123e-02 1.527928e-03 2.796691e-02
#> AVERAGE(COND.)    1.251200e-06 2.938117e-02 1.867651e-06 2.954563e-02
#>                            8.5            9          9.5           10
#> start:Employed_50 2.829862e-06 2.997680e-02 4.139878e-06 3.030336e-02
#> start:Inactive_50 2.356977e-06 2.877136e-02 3.452664e-06 2.900925e-02
#> start:Retired_50  2.812408e-02 1.596881e-06 2.793025e-02 2.352787e-06
#> AVERAGE           1.505427e-03 2.820944e-02 1.496280e-03 2.850641e-02
#> AVERAGE(COND.)    2.761217e-06 2.980181e-02 4.040121e-06 3.011551e-02
#>                           10.5           11         11.5           12
#> start:Employed_50 5.988626e-06 3.065242e-02 8.558830e-06 3.099388e-02
#> start:Inactive_50 5.001558e-06 2.929362e-02 7.159033e-06 2.959354e-02
#> start:Retired_50  2.790611e-02 3.429448e-06 2.799753e-02 4.941850e-06
#> AVERAGE           1.496699e-03 2.882798e-02 1.503961e-03 2.914557e-02
#> AVERAGE(COND.)    5.845342e-06 3.045517e-02 8.355634e-06 3.079061e-02
#>                           12.5           13         13.5           14
#> start:Employed_50 1.207506e-05 3.129812e-02 1.680379e-05 3.153575e-02
#> start:Inactive_50 1.011701e-05 2.987959e-02 1.410466e-05 3.012316e-02
#> start:Retired_50  2.816163e-02 7.034839e-06 2.836215e-02 9.885443e-06
#> AVERAGE           1.515981e-03 2.943116e-02 1.531070e-03 2.965707e-02
#> AVERAGE(COND.)    1.179082e-05 3.109220e-02 1.641198e-05 3.133070e-02
#>                           14.5           15         15.5           16
#> start:Employed_50 2.304866e-05 3.167777e-02 3.113930e-05 3.169594e-02
#> start:Inactive_50 1.938517e-05 3.029629e-02 2.624719e-05 3.037179e-02
#> start:Retired_50  2.856679e-02 1.370273e-05 2.874573e-02 1.872408e-05
#> AVERAGE           1.547783e-03 2.979597e-02 1.564835e-03 2.982132e-02
#> AVERAGE(COND.)    2.251687e-05 3.147723e-02 3.042916e-05 3.150372e-02
#>                           16.5           17         17.5           18
#> start:Employed_50 4.141275e-05 3.156346e-02 5.418761e-05 3.125578e-02
#> start:Inactive_50 3.499010e-05 3.032359e-02 4.590292e-05 3.012743e-02
#> start:Retired_50  2.887086e-02 2.520683e-05 2.891557e-02 3.341472e-05
#> AVERAGE           1.581035e-03 2.970784e-02 1.595261e-03 2.943237e-02
#> AVERAGE(COND.)    4.048043e-05 3.138348e-02 5.298500e-05 3.109199e-02
#>                           18.5           19         19.5           20
#> start:Employed_50 6.973246e-05 3.075164e-02 0.0000882321 3.003417e-02
#> start:Inactive_50 5.923803e-05 2.976156e-02 0.0000751825 2.920773e-02
#> start:Retired_50  2.885486e-02 4.359932e-05 0.0286657524 5.597767e-05
#> AVERAGE           1.606427e-03 2.897471e-02 0.0016134821 2.831873e-02
#> AVERAGE(COND.)    6.820908e-05 3.060792e-02 0.0000863378 2.991420e-02
#>                           20.5           21         21.5           22
#> start:Employed_50 1.097570e-04 2.909211e-02 0.0001342433 2.792103e-02
#> start:Inactive_50 9.383137e-05 2.845221e-02 0.0001151686 2.748687e-02
#> start:Retired_50  2.832791e-02 7.070885e-05 0.0278243957 8.787306e-05
#> AVERAGE           1.615409e-03 2.745343e-02 0.0016112488 2.637411e-02
#> AVERAGE(COND.)    1.074453e-04 2.899923e-02 0.0001314744 2.785801e-02
#>                           22.5           23         23.5           24
#> start:Employed_50 0.0001614903 0.0265244386 0.0001911870 0.0249146833
#> start:Inactive_50 0.0001390629 0.0263102548 0.0001652859 0.0249285104
#> start:Retired_50  0.0271425453 0.0001074585 0.0262748747 0.0001293621
#> AVERAGE           0.0016001444 0.0250834144 0.0015814129 0.0235921767
#> AVERAGE(COND.)    0.0001582347 0.0264933474 0.0001874272 0.0249166904
#>                           24.5           25         25.5          26
#> start:Employed_50 0.0002229740 0.0231135170 0.0002565347 0.021152069
#> start:Inactive_50 0.0001935629 0.0233560943 0.0002236538 0.021616041
#> start:Retired_50  0.0252199567 0.0001534109 0.0239831607 0.000179407
#> AVERAGE           0.0015546494 0.0219199724 0.0015198517 0.020095144
#> AVERAGE(COND.)    0.0002187046 0.0231487299 0.0002517616 0.021219420
#>                           26.5           27         27.5           28
#> start:Employed_50 0.0002916499 0.0190701510 0.0003280609 0.0169148368
#> start:Inactive_50 0.0002554175 0.0197396930 0.0002887361 0.0177658030
#> start:Retired_50  0.0225771417 0.0002071741 0.0210219544 0.0002365445
#> AVERAGE           0.0014774992 0.0181542039 0.0014284381 0.0161405570
#> AVERAGE(COND.)    0.0002863904 0.0191673426 0.0003223525 0.0170383642
#>                           28.5          29         29.5           30
#> start:Employed_50 0.0003650050 0.014738393 0.0004005601 0.0125956646
#> start:Inactive_50 0.0003231639 0.015739049 0.0003573606 0.0137080097
#> start:Retired_50  0.0193447051 0.000267196 0.0175787087 0.0002983306
#> AVERAGE           0.0013734383 0.014102618 0.0013125406 0.0120913965
#> AVERAGE(COND.)    0.0003589313 0.014883650 0.0003942892 0.0127571341
#>                           30.5           31         31.5           32
#> start:Employed_50 0.0004312945 0.0105409077 0.0004527207 0.0086240981
#> start:Inactive_50 0.0003886853 0.0117225912 0.0004133962 0.0098308902
#> start:Retired_50  0.0157621593 0.0003283527 0.0139363120 0.0003548112
#> AVERAGE           0.0012446463 0.0101575674 0.0011678145 0.0083480466
#> AVERAGE(COND.)    0.0004251093 0.0107124424 0.0004470123 0.0087992776
#>                           32.5           33         33.5           34
#> start:Employed_50 0.0004605277 0.0068870027 0.0004519943 0.0053597069
#> start:Inactive_50 0.0004275437 0.0080756777 0.0004281519 0.0064910287
#> start:Retired_50  0.0121431750 0.0003747658 0.0104228006 0.0003854727
#> AVERAGE           0.0010802592 0.0067023500 0.0009815094 0.0052493569
#> AVERAGE(COND.)    0.0004557397 0.0070595523 0.0004485333 0.0055239310
#>                           34.5           35         35.5           36
#> start:Employed_50 0.0004268750 0.0040584157 0.0003874349 0.0029850758
#> start:Inactive_50 0.0004141124 0.0050998018 0.0003864399 0.0039125245
#> start:Retired_50  0.0088104574 0.0003850981 0.0073341282 0.0003731726
#> AVERAGE           0.0008730990 0.0040052228 0.0007584956 0.0029729442
#> AVERAGE(COND.)    0.0004250224 0.0042095846 0.0003872905 0.0031197055
#>                           36.5           37         37.5           38
#> start:Employed_50 0.0003377729 0.0021288407 0.0002828183 0.0014689282
#> start:Inactive_50 0.0003479105 0.0029278554 0.0003023389 0.0021343862
#> start:Retired_50  0.0060127935 0.0003506578 0.0048558014 0.0003196728
#> AVERAGE           0.0006424113 0.0021436116 0.0005298584 0.0014989546
#> AVERAGE(COND.)    0.0003392445 0.0022448267 0.0002856519 0.0015655270
#>                           38.5           39         39.5           40
#> start:Employed_50 0.0002273697 0.0009781790 0.0001753973 0.0006266301
#> start:Inactive_50 0.0002537898 0.0015132871 0.0002059304 0.0010412397
#> start:Retired_50  0.0038633684 0.0002830249 0.0030280256 0.0002437041
#> AVERAGE           0.0004252899 0.0010145596 0.0003320232 0.0006631377
#> AVERAGE(COND.)    0.0002312049 0.0010558560 0.0001798295 0.0006868154
#>                           40.5           41         41.5           42
#> start:Employed_50 0.0001296737 0.0003845869 9.170030e-05 0.0002248998
#> start:Inactive_50 0.0001616196 0.0006931823 1.227366e-04 0.0004445548
#> start:Retired_50  0.0023366742 0.0002044561 1.772870e-03 0.0001674948
#> AVERAGE           0.0002519946 0.0004173640 1.857983e-04 0.0002520139
#> AVERAGE(COND.)    0.0001343110 0.0004293830 9.620557e-05 0.0002567852
#>                           42.5           43         43.5           44
#> start:Employed_50 6.184728e-05 0.0001243628 3.961826e-05 6.430686e-05
#> start:Inactive_50 9.020356e-05 0.0002728925 6.413947e-05 1.587646e-04
#> start:Retired_50  1.318991e-03 0.0001343661 9.580105e-04 1.059437e-04
#> AVERAGE           1.329192e-04 0.0001453060 9.206198e-05 7.951065e-05
#> AVERAGE(COND.)    6.596351e-05 0.0001459236 4.317779e-05 7.801846e-05
#>                           44.5           45         45.5           46
#> start:Employed_50 2.395999e-05 3.055319e-05 1.355697e-05 1.293298e-05
#> start:Inactive_50 4.407919e-05 8.614471e-05 2.920498e-05 4.235844e-05
#> start:Retired_50  6.746922e-04 8.252341e-05 4.561282e-04 6.398295e-05
#> AVERAGE           6.149641e-05 4.096875e-05 3.935592e-05 1.970403e-05
#> AVERAGE(COND.)    2.688052e-05 3.862292e-05 1.582846e-05 1.720442e-05
#>                           46.5           47         47.5           48
#> start:Employed_50 7.071047e-06 4.573533e-06 3.306073e-06 1.118201e-06
#> start:Inactive_50 1.855026e-05 1.776190e-05 1.115119e-05 5.292543e-06
#> start:Retired_50  2.916346e-04 5.000727e-05 1.721469e-04 4.063792e-05
#> AVERAGE           2.385403e-05 8.813433e-06 1.340606e-05 3.803515e-06
#> AVERAGE(COND.)    8.737385e-06 6.487974e-06 4.444881e-06 1.724154e-06
#>                           48.5           49 49.5
#> start:Employed_50 1.294267e-06 0.000000e+00    0
#> start:Inactive_50 6.125878e-06 0.000000e+00    0
#> start:Retired_50  8.923535e-05 1.032859e-04    0
#> AVERAGE           6.657294e-06 5.519092e-06    0
#> AVERAGE(COND.)    1.995630e-06 0.000000e+00    0
  
dtms_visits(dtms=hrs,
            matrix=Tw,
            risk="Retired",
            start_distr=Sw,
            method="end")
#>                            0            1          2          3          4
#> start:Employed_50 0.08275934 2.005943e-07 0.02446359 0.02331318 0.02267940
#> start:Inactive_50 0.08404534 1.602553e-07 0.02580474 0.02396753 0.02289828
#> start:Retired_50  0.00000000 4.940366e-02 0.03703364 0.03034240 0.02655296
#> AVERAGE           0.07926479 2.321703e-03 0.02546557 0.02384416 0.02292855
#> AVERAGE(COND.)    0.08317318 1.876132e-07 0.02489517 0.02352375 0.02274984
#>                            5          6          7          8          9
#> start:Employed_50 0.02240828 0.02240784 0.02261977 0.02300464 0.02353366
#> start:Inactive_50 0.02233017 0.02211278 0.02215617 0.02240317 0.02281492
#> start:Retired_50  0.02435729 0.02310129 0.02243892 0.02217696 0.02220148
#> AVERAGE           0.02247591 0.02234994 0.02246910 0.02278129 0.02325064
#> AVERAGE(COND.)    0.02238314 0.02231289 0.02247059 0.02281109 0.02330237
#>                           10         11         12         13         14
#> start:Employed_50 0.02418394 0.02493545 0.02576902 0.02666502 0.02760228
#> start:Inactive_50 0.02336281 0.02402377 0.02477738 0.02560390 0.02648297
#> start:Retired_50  0.02244130 0.02284880 0.02338949 0.02403598 0.02476434
#> AVERAGE           0.02385023 0.02455780 0.02535309 0.02621605 0.02712566
#> AVERAGE(COND.)    0.02391970 0.02464207 0.02544991 0.02632355 0.02724209
#>                           15         16         17         18         19
#> start:Employed_50 0.02855747 0.02950453 0.03041446 0.03125527 0.03199221
#> start:Inactive_50 0.02739270 0.02830899 0.02920513 0.03005165 0.03081634
#> start:Retired_50  0.02555185 0.02637545 0.02721086 0.02803190 0.02881029
#> AVERAGE           0.02805902 0.02899084 0.02989305 0.03073468 0.03148207
#> AVERAGE(COND.)    0.02818265 0.02911980 0.03002530 0.03086794 0.03161381
#>                           20         21         22         23         24
#> start:Employed_50 0.03258832 0.03300536 0.03320515 0.03315138 0.03281209
#> start:Inactive_50 0.03146462 0.03196024 0.03226633 0.03234694 0.03216915
#> start:Retired_50  0.02951556 0.03011541 0.03057636 0.03086471 0.03094803
#> AVERAGE           0.03209931 0.03254905 0.03279370 0.03279722 0.03252732
#> AVERAGE(COND.)    0.03222671 0.03266905 0.03290304 0.03289251 0.03260519
#>                           25         26         27         28         29
#> start:Employed_50 0.03216257 0.03118845 0.02988779 0.02827096 0.02635877
#> start:Inactive_50 0.03170570 0.03093793 0.02985801 0.02846940 0.02678551
#> start:Retired_50  0.03079711 0.03038818 0.02970511 0.02874068 0.02749672
#> AVERAGE           0.03195830 0.03107402 0.02987007 0.02835389 0.02654311
#> AVERAGE(COND.)    0.03201555 0.03110783 0.02987820 0.02833481 0.02649609
#>                           30         31         32         33         34
#> start:Employed_50 0.02418157 0.02178195 0.01921955 0.01657327 0.01393676
#> start:Inactive_50 0.02482891 0.02263346 0.02024877 0.01774301 0.01520029
#> start:Retired_50  0.02598389 0.02422257 0.02224520 0.02009870 0.01784429
#> AVERAGE           0.02446479 0.02215778 0.01967737 0.01709767 0.01450788
#> AVERAGE(COND.)    0.02438989 0.02205597 0.01955075 0.01694969 0.01434336
#>                           35          36          37          38          39
#> start:Employed_50 0.01140757 0.009073858 0.007003207 0.005236357 0.003786377
#> start:Inactive_50 0.01271202 0.010365284 0.008231970 0.006361694 0.004779304
#> start:Retired_50  0.01555368 0.013302026 0.011159651 0.009184792 0.007418791
#> AVERAGE           0.01200245 0.009668595 0.007575356 0.005767013 0.004261576
#> AVERAGE(COND.)    0.01182734 0.009489438 0.007398622 0.005598490 0.004105900
#>                            40          41          42           43           44
#> start:Employed_50 0.002642174 0.001774584 0.001143319 0.0007035042 0.0004109874
#> start:Inactive_50 0.003486504 0.002466255 0.001688505 0.0011160234 0.0007094946
#> start:Retired_50  0.005884327 0.004586378 0.003515201 0.0026504212 0.0019654289
#> AVERAGE           0.003053463 0.002118833 0.001421972 0.0009215020 0.0005755775
#> AVERAGE(COND.)    0.002913879 0.001997163 0.001318760 0.0008362526 0.0005070468
#>                             45           46           47           48
#> start:Employed_50 0.0002260680 0.0001156258 5.385538e-05 2.194196e-05
#> start:Inactive_50 0.0004314232 0.0002486803 1.337895e-04 6.518272e-05
#> start:Retired_50  0.0014314397 0.0010207935 7.093414e-04 4.786617e-04
#> AVERAGE           0.0003456874 0.0001989653 1.091713e-04 5.666462e-05
#> AVERAGE(COND.)    0.0002921511 0.0001584427 7.957814e-05 3.585681e-05
#>                             49            50
#> start:Employed_50 7.060355e-06 -2.220446e-16
#> start:Inactive_50 2.670889e-05  0.000000e+00
#> start:Retired_50  5.276947e-04  0.000000e+00
#> AVERAGE           3.755122e-05 -1.435144e-16
#> AVERAGE(COND.)    1.338324e-05 -1.505908e-16
  
## First visit
dtms_first(dtms=hrs,
           matrix=Tm,
           risk="Retired",
           start_distr=Sm)  
#>                            0         0.5        1.5        2.5        3.5
#> start:Employed_50 0.00000000 0.008510551 0.01176588 0.01526444 0.01906744
#> start:Inactive_50 0.00000000 0.040384498 0.03705685 0.03536803 0.03492557
#> start:Retired_50  1.00000000 0.000000000 0.00000000 0.00000000 0.00000000
#> AVERAGE           0.06011927 0.012338089 0.01450153 0.01708357 0.02007998
#> AVERAGE(COND.)    0.00000000 0.013127292 0.01542912 0.01817632 0.02136439
#>                          4.5        5.5        6.5        7.5        8.5
#> start:Employed_50 0.02320676 0.02768571 0.03247537 0.03750739 0.04266425
#> start:Inactive_50 0.03546315 0.03678288 0.03871845 0.04110939 0.04378161
#> start:Retired_50  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> AVERAGE           0.02348012 0.02725971 0.03137288 0.03574283 0.04025142
#> AVERAGE(COND.)    0.02498202 0.02900338 0.03337964 0.03802912 0.04282609
#>                          9.5       10.5       11.5       12.5       13.5
#> start:Employed_50 0.04776870 0.05257583 0.05677318 0.05999664 0.06187033
#> start:Inactive_50 0.04653262 0.04912241 0.05127326 0.05268355 0.05306076
#> start:Retired_50  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> AVERAGE           0.04472860 0.04894487 0.05261128 0.05539412 0.05695143
#> AVERAGE(COND.)    0.04758966 0.05207562 0.05597655 0.05893739 0.06059432
#>                         14.5       15.5       16.5       17.5       18.5
#> start:Employed_50 0.06207455 0.06043520 0.05700952 0.05212563 0.04633450
#> start:Inactive_50 0.05217525 0.04992730 0.04640600 0.04190717 0.03688154
#> start:Retired_50  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> AVERAGE           0.05699502 0.05537138 0.05213863 0.04760078 0.04226202
#> AVERAGE(COND.)    0.06064070 0.05891320 0.05547366 0.05064555 0.04496530
#>                         19.5       20.5       21.5       22.5       23.5
#> start:Employed_50 0.04026816 0.03445914 0.02921569 0.02461588 0.02059841
#> start:Inactive_50 0.03181514 0.02709062 0.02290292 0.01927069 0.01611706
#> start:Retired_50  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> AVERAGE           0.03669651 0.03138436 0.02659987 0.02240832 0.01874998
#> AVERAGE(COND.)    0.03904379 0.03339186 0.02830133 0.02384166 0.01994931
#>                         24.5       25.5        26.5        27.5        28.5
#> start:Employed_50 0.01706864 0.01395775 0.011230991 0.008872823 0.006871546
#> start:Inactive_50 0.01335305 0.01091893 0.008785762 0.006941012 0.005375457
#> start:Retired_50  0.00000000 0.00000000 0.000000000 0.000000000 0.000000000
#> AVERAGE           0.01553666 0.01270493 0.010222908 0.008076407 0.006254762
#> AVERAGE(COND.)    0.01653046 0.01351760 0.010876815 0.008593012 0.006654847
#>                          29.5        30.5        31.5        32.5        33.5
#> start:Employed_50 0.005210436 0.003864658 0.002801740 0.001983932 0.001371271
#> start:Inactive_50 0.004076008 0.003023236 0.002191738 0.001551986 0.001072714
#> start:Retired_50  0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> AVERAGE           0.004742752 0.003517770 0.002550258 0.001805856 0.001248187
#> AVERAGE(COND.)    0.005046121 0.003742783 0.002713385 0.001921367 0.001328027
#>                           34.5         35.5         36.5         37.5
#> start:Employed_50 0.0009245231 0.0006075356 0.0003887592 0.0002419632
#> start:Inactive_50 0.0007232337 0.0004752615 0.0003041176 0.0001892824
#> start:Retired_50  0.0000000000 0.0000000000 0.0000000000 0.0000000000
#> AVERAGE           0.0008415388 0.0005530038 0.0003538645 0.0002202448
#> AVERAGE(COND.)    0.0008953676 0.0005883765 0.0003764994 0.0002343328
#>                           38.5         39.5         40.5         41.5
#> start:Employed_50 0.0001462739 8.573963e-05 4.862724e-05 2.661725e-05
#> start:Inactive_50 0.0001144268 6.707219e-05 3.804000e-05 2.082208e-05
#> start:Retired_50  0.0000000000 0.000000e+00 0.000000e+00 0.000000e+00
#> AVERAGE           0.0001331445 7.804372e-05 4.426251e-05 2.422811e-05
#> AVERAGE(COND.)    0.0001416610 8.303577e-05 4.709375e-05 2.577786e-05
#>                           42.5         43.5         44.5         45.5
#> start:Employed_50 1.401955e-05 7.080767e-06 3.415627e-06 1.566582e-06
#> start:Inactive_50 1.096718e-05 5.539126e-06 2.671969e-06 1.225502e-06
#> start:Retired_50  0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> AVERAGE           1.276117e-05 6.445204e-06 3.109044e-06 1.425967e-06
#> AVERAGE(COND.)    1.357743e-05 6.857470e-06 3.307913e-06 1.517179e-06
#>                           46.5         47.5 TOTAL(RESCALED)
#> start:Employed_50 6.797764e-07 2.775573e-07               1
#> start:Inactive_50 5.317739e-07 2.171269e-07               1
#> start:Retired_50  0.000000e+00 0.000000e+00               1
#> AVERAGE           6.187603e-07 2.526440e-07               1
#> AVERAGE(COND.)    6.583392e-07 2.688043e-07               1
  
dtms_first(dtms=hrs,
           matrix=Tw,
           risk="Retired",
           start_distr=Sw)  
#>                            0         0.5        1.5        2.5        3.5
#> start:Employed_50 0.00000000 0.008543757 0.01175132 0.01520147 0.01894215
#> start:Inactive_50 0.00000000 0.032365881 0.03116009 0.03099414 0.03168189
#> start:Retired_50  1.00000000 0.000000000 0.00000000 0.00000000 0.00000000
#> AVERAGE           0.05103632 0.015375496 0.01707291 0.01924375 0.02186212
#> AVERAGE(COND.)    0.00000000 0.016202408 0.01799111 0.02027871 0.02303789
#>                          4.5        5.5        6.5        7.5        8.5
#> start:Employed_50 0.02300192 0.02739042 0.03209511 0.03707431 0.04224656
#> start:Inactive_50 0.03309104 0.03511678 0.03766292 0.04062628 0.04388199
#> start:Retired_50  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> AVERAGE           0.02490603 0.02834971 0.03215575 0.03626582 0.04058940
#> AVERAGE(COND.)    0.02624550 0.02987439 0.03388513 0.03821624 0.04277234
#>                          9.5       10.5       11.5       12.5       13.5
#> start:Employed_50 0.04747706 0.05256289 0.05722147 0.06108921 0.06374106
#> start:Inactive_50 0.04726906 0.05057716 0.05353778 0.05582563 0.05707859
#> start:Retired_50  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> AVERAGE           0.04499054 0.04927446 0.05317726 0.05636560 0.05845533
#> AVERAGE(COND.)    0.04741019 0.05192449 0.05603719 0.05939701 0.06159912
#>                         14.5       15.5       16.5       17.5       18.5
#> start:Employed_50 0.06474218 0.06373837 0.06057602 0.05541516 0.04877345
#> start:Inactive_50 0.05694475 0.05515992 0.05164535 0.04659268 0.04048424
#> start:Retired_50  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> AVERAGE           0.05905909 0.05786824 0.05475983 0.04989536 0.04375531
#> AVERAGE(COND.)    0.06223536 0.06098046 0.05770487 0.05257879 0.04610852
#>                         19.5       20.5       21.5       22.5       23.5
#> start:Employed_50 0.04144050 0.03425935 0.02786711 0.02254152 0.01823599
#> start:Inactive_50 0.03400332 0.02784071 0.02248284 0.01810249 0.01460971
#> start:Retired_50  0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> AVERAGE           0.03705656 0.03055265 0.02480221 0.02003680 0.01619897
#> AVERAGE(COND.)    0.03904950 0.03219580 0.02613610 0.02111440 0.01707017
#>                         24.5        25.5        26.5        27.5        28.5
#> start:Employed_50 0.01474438 0.011858851 0.009438700 0.007404906 0.005711551
#> start:Inactive_50 0.01180073 0.009488268 0.007551315 0.005924116 0.004569379
#> start:Retired_50  0.00000000 0.000000000 0.000000000 0.000000000 0.000000000
#> AVERAGE           0.01309382 0.010530389 0.008381170 0.006575219 0.005071594
#> AVERAGE(COND.)    0.01379802 0.011096725 0.008831918 0.006928842 0.005344350
#>                          29.5        30.5        31.5        32.5         33.5
#> start:Employed_50 0.004324293 0.003210304 0.002335277 0.001663703 0.0011603683
#> start:Inactive_50 0.003459539 0.002568321 0.001868278 0.001331002 0.0009283226
#> start:Retired_50  0.000000000 0.000000000 0.000000000 0.000000000 0.0000000000
#> AVERAGE           0.003839773 0.002850602 0.002073619 0.001477291 0.0010303536
#> AVERAGE(COND.)    0.004046281 0.003003911 0.002185140 0.001556742 0.0010857672
#>                           34.5         35.5         36.5         37.5
#> start:Employed_50 0.0007920612 0.0005289541 0.0003454565 0.0002205181
#> start:Inactive_50 0.0006336680 0.0004231760 0.0002763735 0.0001764198
#> start:Retired_50  0.0000000000 0.0000000000 0.0000000000 0.0000000000
#> AVERAGE           0.0007033139 0.0004696868 0.0003067494 0.0001958099
#> AVERAGE(COND.)    0.0007411389 0.0004949471 0.0003232468 0.0002063408
#>                           38.5         39.5         40.5         41.5
#> start:Employed_50 0.0001374813 8.362788e-05 4.956649e-05 2.857676e-05
#> start:Inactive_50 0.0001099883 6.690432e-05 3.965439e-05 2.286210e-05
#> start:Retired_50  0.0000000000 0.000000e+00 0.000000e+00 0.000000e+00
#> AVERAGE           0.0001220771 7.425771e-05 4.401276e-05 2.537485e-05
#> AVERAGE(COND.)    0.0001286425 7.825137e-05 4.637982e-05 2.673953e-05
#>                           42.5         43.5         44.5         45.5
#> start:Employed_50 1.599180e-05 8.663770e-06 4.529845e-06 2.277371e-06
#> start:Inactive_50 1.279383e-05 6.931225e-06 3.623985e-06 1.821952e-06
#> start:Retired_50  0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> AVERAGE           1.419999e-05 7.693029e-06 4.022294e-06 2.022201e-06
#> AVERAGE(COND.)    1.496368e-05 8.106769e-06 4.238617e-06 2.130957e-06
#>                           46.5         47.5 TOTAL(RESCALED)
#> start:Employed_50 1.096292e-06 5.029110e-07               1
#> start:Inactive_50 8.770597e-07 4.023409e-07               1
#> start:Retired_50  0.000000e+00 0.000000e+00               1
#> AVERAGE           9.734565e-07 4.465618e-07               1
#> AVERAGE(COND.)    1.025810e-06 4.705784e-07               1

## Last exit
  
# Leaving employment to any state
dtms_last(dtms=hrs,
          matrix=Tm,
          risk="Employed",
          start_distr=Sm)  
#>                          0.5         1.5         2.5        3.5        4.5
#> start:Employed_50 0.01291930 0.015581788 0.018435838 0.02139260 0.02435957
#> start:Inactive_50 0.00000000 0.003371822 0.007140880 0.01120125 0.01543120
#> start:Retired_50  0.00000000 0.004412367 0.008671143 0.01289217 0.01708130
#> AVERAGE           0.01084782 0.013669976 0.016692358 0.01983317 0.02300084
#> AVERAGE(COND.)    0.00000000 0.003658321 0.007562215 0.01166682 0.01588553
#>                          5.5        6.5        7.5        8.5        9.5
#> start:Employed_50 0.02724904 0.02998559 0.03251109 0.03478660 0.03679097
#> start:Inactive_50 0.01970305 0.02389432 0.02789749 0.03162640 0.03501854
#> start:Retired_50  0.02118950 0.02514605 0.02887925 0.03232845 0.03544917
#> AVERAGE           0.02610474 0.02906418 0.03181469 0.03431089 0.03652579
#> AVERAGE(COND.)    0.02011232 0.02423897 0.02816780 0.03181970 0.03513711
#>                         10.5       11.5       12.5       13.5       14.5
#> start:Employed_50 0.03851654 0.03996294 0.04113005 0.04201154 0.04259039
#> start:Inactive_50 0.03803302 0.04064518 0.04283914 0.04460008 0.04590811
#> start:Retired_50  0.03821283 0.04060258 0.04260675 0.04421206 0.04539813
#> AVERAGE           0.03844695 0.04007045 0.04139383 0.04240945 0.04309984
#> AVERAGE(COND.)    0.03808253 0.04063345 0.04277516 0.04449324 0.04576770
#>                         15.5       16.5       17.5       18.5       19.5
#> start:Employed_50 0.04283776 0.04271613 0.04218708 0.04122217 0.03981376
#> start:Inactive_50 0.04673563 0.04704943 0.04681799 0.04602227 0.04466652
#> start:Retired_50  0.04613551 0.04638821 0.04612133 0.04531210 0.04396098
#> AVERAGE           0.04343625 0.04338174 0.04289885 0.04196047 0.04056070
#> AVERAGE(COND.)    0.04657039 0.04686737 0.04662618 0.04582674 0.04447226
#>                         20.5       21.5       22.5       23.5       24.5
#> start:Employed_50 0.03798195 0.03577477 0.03326223 0.03052757 0.02765939
#> start:Inactive_50 0.04278485 0.04044063 0.03771962 0.03472048 0.03154663
#> start:Retired_50  0.04209843 0.03978455 0.03710226 0.03414785 0.03102267
#> AVERAGE           0.03872174 0.03649393 0.03394967 0.03117458 0.02825954
#> AVERAGE(COND.)    0.04259585 0.04025999 0.03754964 0.03456282 0.03140237
#>                         25.5       26.5       27.5       28.5       29.5
#> start:Employed_50 0.02474592 0.02187065 0.01910816 0.01652012 0.01415214
#> start:Inactive_50 0.02830050 0.02507905 0.02196927 0.01904370 0.01635671
#> start:Retired_50  0.02782728 0.02465695 0.02159714 0.01871908 0.01607615
#> AVERAGE           0.02529497 0.02236645 0.01955048 0.01691042 0.01449324
#> AVERAGE(COND.)    0.02817020 0.02496283 0.02186681 0.01895432 0.01627946
#>                         30.5       31.5        32.5        33.5        34.5
#> start:Employed_50 0.01203229 0.01017145 0.008565468 0.007198436 0.006046463
#> start:Inactive_50 0.01394271 0.01181644 0.009975293 0.008402934 0.007073609
#> start:Retired_50  0.01370209 0.01161130 0.009801119 0.008255421 0.006948812
#> AVERAGE           0.01232798 0.01042615 0.008783831 0.007385053 0.006205646
#> AVERAGE(COND.)    0.01387646 0.01175996 0.009927337 0.008362319 0.007039248
#>                          35.5        36.5        37.5        38.5        39.5
#> start:Employed_50 0.005081374 0.004273862 0.003595900 0.003022358 0.002531929
#> start:Inactive_50 0.005956321 0.005018470 0.004228642 0.003558524 0.002984010
#> start:Retired_50  0.005850764 0.004929184 0.004153159 0.003494829 0.002930481
#> AVERAGE           0.005217003 0.004389311 0.003694022 0.003105515 0.002602052
#> AVERAGE(COND.)    0.005927257 0.004993886 0.004207859 0.003540986 0.002969272
#>                          40.5        41.5        42.5        43.5         44.5
#> start:Employed_50 0.002107505 0.001736165 0.001408898 0.001120128 0.0008670725
#> start:Inactive_50 0.002485693 0.002048894 0.001663383 0.001322859 0.0010242250
#> start:Retired_50  0.002441028 0.002012032 0.001633428 0.001299020 0.0010057589
#> AVERAGE           0.002166172 0.001784681 0.001448380 0.001151582 0.0008914552
#> AVERAGE(COND.)    0.002473395 0.002038745 0.001655136 0.001316295 0.0010191406
#>                           45.5         46.5         47.5         48.5
#> start:Employed_50 0.0006489230 0.0004658608 0.0003179800 0.0002042865
#> start:Inactive_50 0.0007666527 0.0005504360 0.0003757345 0.0002414024
#> start:Retired_50  0.0007528259 0.0005405064 0.0003689554 0.0002370465
#> AVERAGE           0.0006671894 0.0004789833 0.0003269411 0.0002100453
#> AVERAGE(COND.)    0.0007628457 0.0005477020 0.0003738680 0.0002402031
#>                   TOTAL(RESCALED)
#> start:Employed_50               1
#> start:Inactive_50               1
#> start:Retired_50                1
#> AVERAGE                         1
#> AVERAGE(COND.)                  1
  
dtms_last(dtms=hrs,
          matrix=Tw,
          risk="Employed",
          start_distr=Sw)  
#>                          0.5         1.5         2.5        3.5        4.5
#> start:Employed_50 0.01766783 0.020613619 0.023555315 0.02638148 0.02899653
#> start:Inactive_50 0.00000000 0.003708492 0.007807114 0.01212918 0.01650949
#> start:Retired_50  0.00000000 0.005504919 0.010474444 0.01510825 0.01945313
#> AVERAGE           0.01236192 0.015609013 0.018933181 0.02222113 0.02536489
#> AVERAGE(COND.)    0.00000000 0.003949094 0.008164359 0.01252817 0.01690375
#>                          5.5        6.5        7.5        8.5        9.5
#> start:Employed_50 0.03132930 0.03333782 0.03500950 0.03635722 0.03741231
#> start:Inactive_50 0.02079613 0.02486179 0.02861160 0.03198569 0.03495657
#> start:Retired_50  0.02348879 0.02717535 0.03047892 0.03338316 0.03589142
#> AVERAGE           0.02827434 0.03088540 0.03316323 0.03510059 0.03671242
#> AVERAGE(COND.)    0.02115677 0.02517166 0.02886170 0.03217285 0.03508178
#>                         10.5       11.5       12.5       13.5       14.5
#> start:Employed_50 0.03821572 0.03880879 0.03922461 0.03948108 0.03957618
#> start:Inactive_50 0.03752221 0.03969621 0.04149697 0.04293734 0.04401647
#> start:Retired_50  0.03802264 0.03980372 0.04126049 0.04240896 0.04324858
#> AVERAGE           0.03802758 0.03907962 0.03989752 0.04049779 0.04087878
#> AVERAGE(COND.)    0.03758924 0.03971061 0.04146529 0.04286657 0.04391362
#>                         15.5       16.5       17.5       18.5       19.5
#> start:Employed_50 0.03948656 0.03917027 0.03857452 0.03764826 0.03635742
#> start:Inactive_50 0.04471538 0.04499779 0.04481752 0.04413194 0.04291898
#> start:Retired_50  0.04375890 0.04390115 0.04362597 0.04288692 0.04165733
#> AVERAGE           0.04101838 0.04087625 0.04040146 0.03954533 0.03827720
#> AVERAGE(COND.)    0.04458727 0.04485092 0.04465793 0.04396519 0.04275000
#>                         20.5       21.5       22.5       23.5       24.5
#> start:Employed_50 0.03469816 0.03270274 0.03043490 0.02797735 0.02541750
#> start:Inactive_50 0.04119170 0.03900417 0.03644520 0.03362298 0.03064867
#> start:Retired_50  0.03994487 0.03779763 0.03529832 0.03254947 0.02965736
#> AVERAGE           0.03659812 0.03454662 0.03219376 0.02962964 0.02694862
#> AVERAGE(COND.)    0.04102471 0.03884258 0.03629159 0.03347920 0.03051590
#>                         25.5       26.5       27.5       28.5       29.5
#> start:Employed_50 0.02283752 0.02030970 0.01789480 0.01564100 0.01358284
#> start:Inactive_50 0.02762527 0.02464299 0.02177785 0.01909070 0.01662595
#> start:Retired_50  0.02672093 0.02382699 0.02104876 0.01844477 0.01605764
#> AVERAGE           0.02423898 0.02157823 0.01903161 0.01665102 0.01447387
#> AVERAGE(COND.)    0.02750415 0.02453370 0.02168020 0.01900419 0.01654984
#>                         30.5       31.5        32.5        33.5        34.5
#> start:Employed_50 0.01174068 0.01012123 0.008719335 0.007520512 0.006503987
#> start:Inactive_50 0.01441085 0.01245604 0.010757569 0.009299914 0.008059579
#> start:Retired_50  0.01391343 0.01202210 0.010379561 0.008970547 0.007772131
#> AVERAGE           0.01252256 0.01080496 0.009316243 0.008041645 0.006959593
#> AVERAGE(COND.)    0.01434423 0.01239792 0.010706941 0.009255800 0.008021080
#>                          35.5        36.5        37.5        38.5        39.5
#> start:Employed_50 0.005645577 0.004920159 0.004303547 0.003773729 0.003311576
#> start:Inactive_50 0.007008582 0.006117469 0.005357622 0.004702820 0.004130147
#> start:Retired_50  0.006757092 0.005896824 0.005163569 0.004531912 0.003979662
#> AVERAGE           0.006044792 0.005270855 0.004612296 0.004045875 0.003551352
#> AVERAGE(COND.)    0.006974899 0.006087917 0.005331632 0.004679929 0.004109992
#>                          40.5        41.5        42.5        43.5        44.5
#> start:Employed_50 0.002901139 0.002529706 0.002187719 0.001868660 0.001568931
#> start:Inactive_50 0.003620419 0.003158287 0.002732189 0.002334248 0.001960146
#> start:Retired_50  0.003488249 0.003042822 0.002632199 0.002248759 0.001888321
#> AVERAGE           0.003111834 0.002713834 0.002347210 0.002005045 0.001683530
#> AVERAGE(COND.)    0.003602717 0.003142822 0.002718797 0.002322799 0.001950526
#>                          45.5        46.5         47.5         48.5
#> start:Employed_50 0.001287683 0.001026522 0.0007889588 0.0005794964
#> start:Inactive_50 0.001608940 0.001282717 0.0009859112 0.0007241830
#> start:Retired_50  0.001549964 0.001235687 0.0009497582 0.0006976246
#> AVERAGE           0.001381789 0.001101570 0.0008466524 0.0006218797
#> AVERAGE(COND.)    0.001601041 0.001276418 0.0009810691 0.0007206259
#>                   TOTAL(RESCALED)
#> start:Employed_50               1
#> start:Inactive_50               1
#> start:Retired_50                1
#> AVERAGE                         1
#> AVERAGE(COND.)                  1
  
# Leaving employment for retirement
dtms_last(dtms=hrs,
          matrix=Tm,
          risk="Employed",
          risk_to="Retired",
          start_distr=Sm)  
#>                           0.5         1.5         2.5         3.5         4.5
#> start:Employed_50 0.006374811 0.007556942 0.008949602 0.010566320 0.012418675
#> start:Inactive_50 0.000000000 0.001548804 0.003283188 0.005239979 0.007450894
#> start:Retired_50  0.000000000 0.002041284 0.004015319 0.006074200 0.008306721
#> AVERAGE           0.005306920 0.006573072 0.008033973 0.009712344 0.011625755
#> AVERAGE(COND.)    0.000000000 0.001683701 0.003483728 0.005468483 0.007685316
#>                           5.5        6.5        7.5        8.5        9.5
#> start:Employed_50 0.014514853 0.01685850 0.01944769 0.02227364 0.02531878
#> start:Inactive_50 0.009940264 0.01272343 0.01580537 0.01917927 0.02282457
#> start:Retired_50  0.010766760 0.01348587 0.01647879 0.01974545 0.02327075
#> AVERAGE           0.013786454 0.01620079 0.01886844 0.02178126 0.02492143
#> AVERAGE(COND.)    0.010166653 0.01293227 0.01598983 0.01933435 0.02294679
#>                         10.5       11.5       12.5       13.5       14.5
#> start:Employed_50 0.02855360 0.03193213 0.03538631 0.03882020 0.04210632
#> start:Inactive_50 0.02670407 0.03075974 0.03490760 0.03903264 0.04298611
#> start:Retired_50  0.02702251 0.03094760 0.03496692 0.03897021 0.04281308
#> AVERAGE           0.02825838 0.03174436 0.03530884 0.03885292 0.04224576
#> AVERAGE(COND.)    0.02679130 0.03081120 0.03492385 0.03901554 0.04293872
#>                         15.5       16.5       17.5       18.5       19.5
#> start:Employed_50 0.04508681 0.04758253 0.04941167 0.05041617 0.05048984
#> start:Inactive_50 0.04658799 0.04963787 0.05193570 0.05331018 0.05364831
#> start:Retired_50  0.04631919 0.04929084 0.05152937 0.05286352 0.05317911
#> AVERAGE           0.04532595 0.04791092 0.04981585 0.05088047 0.05099740
#> AVERAGE(COND.)    0.04651436 0.04954281 0.05182440 0.05318784 0.05351979
#>                         20.5       21.5       22.5       23.5       24.5
#> start:Employed_50 0.04959917 0.04778926 0.04517363 0.04191371 0.03819621
#> start:Inactive_50 0.05291640 0.05116518 0.04851814 0.04514948 0.04126043
#> start:Retired_50  0.05244040 0.05069566 0.04806589 0.04472292 0.04086577
#> AVERAGE           0.05013302 0.04833324 0.04571314 0.04243619 0.03869141
#> AVERAGE(COND.)    0.05278601 0.05103658 0.04839426 0.04503264 0.04115233
#>                         25.5       26.5       27.5       28.5       29.5
#> start:Employed_50 0.03421361 0.03014881 0.02616306 0.02238720 0.01891658
#> start:Inactive_50 0.03705889 0.03274330 0.02848974 0.02444225 0.02070711
#> start:Retired_50  0.03670024 0.03242280 0.02820777 0.02419770 0.02049771
#> AVERAGE           0.03467379 0.03056872 0.02653988 0.02272024 0.01920691
#> AVERAGE(COND.)    0.03696065 0.03265551 0.02841250 0.02437527 0.02064976
#>                         30.5       31.5       32.5        33.5        34.5
#> start:Employed_50 0.01581002 0.01309268 0.01076171 0.008793454 0.007151018
#> start:Inactive_50 0.01735140 0.01440574 0.01187022 0.009721997 0.007923382
#> start:Retired_50  0.01717409 0.01425704 0.01174650 0.009619744 0.007839347
#> AVERAGE           0.01606009 0.01330582 0.01094173 0.008944309 0.007276546
#> AVERAGE(COND.)    0.01730283 0.01436501 0.01183633 0.009693988 0.007900364
#>                          35.5        36.5        37.5        38.5        39.5
#> start:Employed_50 0.005791071 0.004669362 0.003744616 0.002980886 0.002348590
#> start:Inactive_50 0.006429230 0.005192918 0.004170649 0.003324087 0.002621556
#> start:Retired_50  0.006360529 0.005137065 0.004125542 0.003287972 0.002592971
#> AVERAGE           0.005894822 0.004754504 0.003813914 0.003036721 0.002393005
#> AVERAGE(COND.)    0.006410412 0.005177619 0.004158294 0.003314195 0.002613726
#>                          40.5        41.5        42.5         43.5         44.5
#> start:Employed_50 0.001824526 0.001391178 0.001035606 0.0007481226 0.0005209633
#> start:Inactive_50 0.002038130 0.001554942 0.001158005 0.0008368000 0.0005828411
#> start:Retired_50  0.002015845 0.001537904 0.001145297 0.0008276064 0.0005764326
#> AVERAGE           0.001859285 0.001417829 0.001055527 0.0007625558 0.0005310349
#> AVERAGE(COND.)    0.002032026 0.001550275 0.001154524 0.0008342818 0.0005810858
#>                           45.5         46.5         47.5         48.5
#> start:Employed_50 0.0003471251 0.0002195238 0.0001305913 7.234773e-05
#> start:Inactive_50 0.0003884138 0.0002456605 0.0001461500 8.097110e-05
#> start:Retired_50  0.0003841407 0.0002429569 0.0001445411 8.007958e-05
#> AVERAGE           0.0003538456 0.0002237781 0.0001331238 7.375138e-05
#> AVERAGE(COND.)    0.0003872433 0.0002449199 0.0001457093 8.072690e-05
#>                   TOTAL(RESCALED)
#> start:Employed_50               1
#> start:Inactive_50               1
#> start:Retired_50                1
#> AVERAGE                         1
#> AVERAGE(COND.)                  1
  
dtms_last(dtms=hrs,
          matrix=Tw,
          risk="Employed",
          risk_to="Retired",
          start_distr=Sw)  
#>                           0.5         1.5         2.5         3.5         4.5
#> start:Employed_50 0.007325823 0.008464907 0.009760816 0.011217980 0.012841688
#> start:Inactive_50 0.000000000 0.001390216 0.002953280 0.004708294 0.006674623
#> start:Retired_50  0.000000000 0.002096664 0.004025669 0.005958534 0.007990528
#> AVERAGE           0.004986331 0.006235416 0.007632078 0.009191864 0.010927758
#> AVERAGE(COND.)    0.000000000 0.001483540 0.003094946 0.004873455 0.006848459
#>                           5.5        6.5        7.5        8.5        9.5
#> start:Employed_50 0.014637676 0.01661185 0.01876998 0.02111699 0.02365542
#> start:Inactive_50 0.008869943 0.01130915 0.01400351 0.01695953 0.02017724
#> start:Retired_50  0.010178692 0.01255931 0.01515610 0.01798369 0.02104828
#> AVERAGE           0.012850970 0.01497118 0.01729644 0.01983252 0.02258141
#> AVERAGE(COND.)    0.009042833 0.01147430 0.01415577 0.01709482 0.02029231
#>                         10.5       11.5       12.5       13.5       14.5
#> start:Employed_50 0.02638242 0.02928480 0.03233197 0.03546735 0.03859977
#> start:Inactive_50 0.02364710 0.02734501 0.03122530 0.03521208 0.03919069
#> start:Retired_50  0.02434584 0.02785773 0.03154407 0.03533518 0.03912305
#> AVERAGE           0.02553837 0.02868696 0.03199200 0.03539103 0.03878563
#> AVERAGE(COND.)    0.02373940 0.02741274 0.03126742 0.03522835 0.03918176
#>                         15.5       16.5       17.5       18.5       19.5
#> start:Employed_50 0.04159773 0.04429066 0.04648140 0.04797178 0.04859845
#> start:Inactive_50 0.04300253 0.04644764 0.04929959 0.05133468 0.05237158
#> start:Retired_50  0.04275596 0.04604065 0.04875664 0.05068457 0.05164530
#> AVERAGE           0.04203595 0.04496232 0.04735848 0.04901829 0.04977275
#> AVERAGE(COND.)    0.04296996 0.04639388 0.04922787 0.05124880 0.05227563
#>                         20.5       21.5       22.5       23.5       24.5
#> start:Employed_50 0.04826892 0.04698400 0.04483650 0.04198762 0.03863311
#> start:Inactive_50 0.05231036 0.05115568 0.04901365 0.04606466 0.04252608
#> start:Retired_50  0.05153854 0.05036634 0.04823072 0.04530735 0.04180896
#> AVERAGE           0.04952699 0.04828292 0.04613744 0.04325767 0.03984607
#> AVERAGE(COND.)    0.05220840 0.05105140 0.04891022 0.04596461 0.04243135
#>                         25.5       26.5       27.5       28.5       29.5
#> start:Employed_50 0.03497347 0.03119439 0.02745575 0.02388563 0.02057743
#> start:Inactive_50 0.03862009 0.03455280 0.03050270 0.02661406 0.02299344
#> start:Retired_50  0.03795345 0.03394315 0.02995318 0.02612495 0.02256277
#> AVERAGE           0.03610989 0.03224117 0.02840560 0.02473632 0.02133081
#> AVERAGE(COND.)    0.03853202 0.03447226 0.03043011 0.02654944 0.02293655
#>                         30.5       31.5       32.5       33.5        34.5
#> start:Employed_50 0.01758983 0.01494981 0.01265790 0.01069488 0.009028591
#> start:Inactive_50 0.01970947 0.01679574 0.01425638 0.01207325 0.010213380
#> start:Retired_50  0.01933359 0.01646996 0.01397550 0.01183197 0.010006688
#> AVERAGE           0.01825088 0.01552556 0.01315652 0.01112488 0.009398232
#> AVERAGE(COND.)    0.01965982 0.01675270 0.01421927 0.01204137 0.010186076
#>                          35.5        36.5        37.5        38.5        39.5
#> start:Employed_50 0.007620264 0.006429561 0.005418289 0.004552758 0.003805006
#> start:Inactive_50 0.008635922 0.007297780 0.006157785 0.005179395 0.004332146
#> start:Retired_50  0.008459243 0.007147106 0.006029696 0.005071019 0.004241083
#> AVERAGE           0.007937160 0.006700469 0.005649042 0.004748302 0.003969506
#> AVERAGE(COND.)    0.008612582 0.007277876 0.006140864 0.005165078 0.004320117
#>                          40.5        41.5        42.5        43.5        44.5
#> start:Employed_50 0.003153160 0.002581193 0.002078273 0.001637830 0.001256424
#> start:Inactive_50 0.003592139 0.002941839 0.002369402 0.001867680 0.001432972
#> start:Retired_50  0.003516372 0.002879631 0.002319208 0.001828064 0.001402550
#> AVERAGE           0.003290151 0.002693740 0.002169127 0.001709561 0.001311521
#> AVERAGE(COND.)    0.003582130 0.002933621 0.002362772 0.001862446 0.001428953
#>                           45.5         46.5         47.5         48.5
#> start:Employed_50 0.0009324810 0.0006649574 0.0004520926 0.0002904283
#> start:Inactive_50 0.0010636232 0.0007585305 0.0005157364 0.0003313244
#> start:Retired_50  0.0010410284 0.0007424103 0.0005047729 0.0003242799
#> AVERAGE           0.0009734079 0.0006941598 0.0004719547 0.0003031913
#> AVERAGE(COND.)    0.0010606383 0.0007564010 0.0005142881 0.0003303938
#>                   TOTAL(RESCALED)
#> start:Employed_50               1
#> start:Inactive_50               1
#> start:Retired_50                1
#> AVERAGE                         1
#> AVERAGE(COND.)                  1
```
