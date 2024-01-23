
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dtms

<!-- badges: start -->
<!-- badges: end -->

## Authors

Christian Dudel, <dudel@demogr.mpg.de>

Peng Li, <li@demogr.mpg.de>

## Overview

The package dtms is a user-friendly implementation of discrete-time
multistate models in R. It comes with many tools to analyze the results
of multistate models. In particular, the worklfow mainly consists of
estimating a discrete-time multistate model and then applying methods
for absorbing Markov chains. Discrete-time multistate models lend
themselves well for modelling trajectories captured in panel data.

## Disclaimer

This package is still undergoing development and many functions are
still experimental. The content of this repository will likely change in
the future, and functions and features might be added or removed without
warning.

## Installation

You can install the development version of dtms like this:

``` r
library(devtools)
install_github("christiandudel/dtms")
```

## Example 1: Artificial data

This is a very basic example using artificial data provided with the
package. The state space consists of two transient states (A, B), and
there is one absorbing state (X). The time scale goes from 0 to 20.
Transition probabilities do change depending on time.

``` r
## Load package
library(dtms)

## Define model: Absorbing and transient states, time scale
simple <- dtms(transient=c("A","B"),
               absorbing="X",
               timescale=0:20)

## Look at original data
head(simpledata)
#>   id time state
#> 1  2    0     A
#> 2  2    1     B
#> 3  2    2     A
#> 4  2    3     B
#> 5  2    4     A
#> 6  2    5     A

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
#> # A tibble: 6 × 4
#>      id  time from  to   
#>   <dbl> <dbl> <chr> <chr>
#> 1     2     0 A     B    
#> 2     2     1 B     A    
#> 3     2     2 A     B    
#> 4     2     3 B     A    
#> 5     2     4 A     A    
#> 6     2     5 A     B

## Clean
estdata <- dtms_clean(data=estdata,
                      dtms=simple)
#> Dropping  0  rows because not in time range
#> Dropping  1215  rows because gap, last obs, ...
#> Dropping  0  rows because starting in absorbing state

# Fit model 
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

## State expectancies 
dtms_expectancy(dtms=simple,
                matrix=Tp,
                start_distr=S)
#>                  A        B    TOTAL
#> start:A_0 5.015765 8.514421 13.53019
#> start:B_0 4.799698 8.711435 13.51113
#> AVERAGE   4.905684 8.614794 13.52048

## Lifetime risk 
dtms_risk(dtms=simple,
          matrix=Tp,
          risk="A")
#>       A_0       B_0 
#> 1.0000000 0.9846722

## Distribution of visits
dtms_visits(dtms=simple,
            matrix=Tp,
            risk="A",
            start_distr=S)
#>                          0        0.5          1        1.5          2
#> start:A_0      0.000000000 0.02139860 0.00000000 0.04556435 0.00000000
#> start:B_0      0.015327809 0.00000000 0.03727791 0.00000000 0.07359658
#> AVERAGE        0.007809115 0.01049658 0.01899212 0.02235051 0.03749552
#> AVERAGE(COND.) 0.015327809 0.00000000 0.03727791 0.00000000 0.07359658
#>                         2.5          3        3.5          4           4.5
#> start:A_0      8.732850e-02 0.00000000 0.14249972 0.00000000  1.894252e-01
#> start:B_0      2.775558e-17 0.12565518 0.00000000 0.17717692 -5.551115e-17
#> AVERAGE        4.283692e-02 0.06401801 0.06989986 0.09026698  9.291803e-02
#> AVERAGE(COND.) 2.775558e-17 0.12565518 0.00000000 0.17717692 -5.551115e-17
#>                            5       5.5             6           6.5          7
#> start:A_0      -1.110223e-16 0.1989227 -1.110223e-16  1.611155e-01 0.00000000
#> start:B_0       1.995772e-01 0.0000000  1.752243e-01 -1.110223e-16 0.11674654
#> AVERAGE         1.016793e-01 0.0975768  8.927219e-02  7.903140e-02 0.05947929
#> AVERAGE(COND.)  1.995772e-01 0.0000000  1.752243e-01 -1.110223e-16 0.11674654
#>                       7.5          8        8.5           9         9.5
#> start:A_0      0.09760683 0.00000000 0.04197567 0.000000000 0.011882873
#> start:B_0      0.00000000 0.05652125 0.00000000 0.018572979 0.000000000
#> AVERAGE        0.04787872 0.02879609 0.02059017 0.009462444 0.005828862
#> AVERAGE(COND.) 0.00000000 0.05652125 0.00000000 0.018572979 0.000000000
#>                         10        10.5           11         11.5           12
#> start:A_0      0.000000000 0.002056761 0.0000000000 0.0002100808 0.000000e+00
#> start:B_0      0.003823249 0.000000000 0.0004654104 0.0000000000 3.317334e-05
#> AVERAGE        0.001947845 0.001008895 0.0002371144 0.0001030502 1.690095e-05
#> AVERAGE(COND.) 0.003823249 0.000000000 0.0004654104 0.0000000000 3.317334e-05
#>                        12.5           13          13.5            14
#> start:A_0      1.278861e-05 0.000000e+00  4.781199e-07 -1.110223e-16
#> start:B_0      0.000000e+00 1.417303e-06 -1.110223e-16  3.736743e-08
#> AVERAGE        6.273150e-06 7.220785e-07  2.345304e-07  1.903772e-08
#> AVERAGE(COND.) 0.000000e+00 1.417303e-06 -1.110223e-16  3.736743e-08
#>                        14.5            15          15.5            16
#> start:A_0      1.128639e-08 -1.110223e-16  1.706568e-10 -1.110223e-16
#> start:B_0      0.000000e+00  6.191081e-10 -1.110223e-16  6.434520e-12
#> AVERAGE        5.536273e-09  3.154192e-10  8.371161e-11  3.278164e-12
#> AVERAGE(COND.) 0.000000e+00  6.191081e-10 -1.110223e-16  6.434520e-12
#>                         16.5            17         17.5            18
#> start:A_0       1.643796e-12 -3.330669e-16 9.992007e-15 -1.110223e-16
#> start:B_0      -2.220446e-16  4.063416e-14 0.000000e+00  2.220446e-16
#> AVERAGE         8.062122e-13  2.053866e-14 4.901342e-15  5.866652e-17
#> AVERAGE(COND.) -2.220446e-16  4.063416e-14 0.000000e+00  2.220446e-16
#>                        18.5           19         19.5            20
#> start:A_0      1.110223e-16 2.220446e-16 0.000000e+00  0.000000e+00
#> start:B_0      0.000000e+00 0.000000e+00 2.220446e-16 -2.220446e-16
#> AVERAGE        5.445936e-17 1.089187e-16 1.131259e-16 -1.131259e-16
#> AVERAGE(COND.) 0.000000e+00 0.000000e+00 2.220446e-16 -2.220446e-16
#>                         20.5
#> start:A_0      -2.220446e-16
#> start:B_0      -1.110223e-16
#> AVERAGE        -1.654817e-16
#> AVERAGE(COND.) -1.110223e-16

## First/last visit
dtms_first(dtms=simple,
           matrix=Tp,
           risk="A",
           start_distr=S)
#>                       0       0.5       1.5        2.5        3.5        4.5
#> start:A_0      1.000000 0.0000000 0.0000000 0.00000000 0.00000000 0.00000000
#> start:B_0      0.000000 0.5298969 0.2487133 0.11704870 0.05519050 0.02604746
#> AVERAGE        0.494387 0.2679227 0.1257527 0.05918134 0.02790503 0.01316994
#> AVERAGE(COND.) 0.000000 0.5298969 0.2487133 0.11704870 0.05519050 0.02604746
#>                        5.5         6.5         7.5         8.5          9.5
#> start:A_0      0.000000000 0.000000000 0.000000000 0.000000000 0.0000000000
#> start:B_0      0.012289033 0.005786386 0.002713402 0.001263733 0.0005825315
#> AVERAGE        0.006213494 0.002925672 0.001371931 0.000638960 0.0002945355
#> AVERAGE(COND.) 0.012289033 0.005786386 0.002713402 0.001263733 0.0005825315
#>                        10.5         11.5         12.5         13.5         14.5
#> start:A_0      0.0000000000 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> start:B_0      0.0002645973 1.177700e-04 5.100991e-05 2.131786e-05 8.508056e-06
#> AVERAGE        0.0001337838 5.954603e-05 2.579127e-05 1.077858e-05 4.301783e-06
#> AVERAGE(COND.) 0.0002645973 1.177700e-04 5.100991e-05 2.131786e-05 8.508056e-06
#>                        15.5         16.5         17.5         18.5
#> start:A_0      0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> start:B_0      3.203492e-06 1.122048e-06 3.598564e-07 1.038734e-07
#> AVERAGE        1.619727e-06 5.673221e-07 1.819480e-07 5.251973e-08
#> AVERAGE(COND.) 3.203492e-06 1.122048e-06 3.598564e-07 1.038734e-07
#>                TOTAL(RESCALED)
#> start:A_0                    1
#> start:B_0                    1
#> AVERAGE                      1
#> AVERAGE(COND.)               1

dtms_last(dtms=simple,
          matrix=Tp,
          risk="A",
          start_distr=S,
          rescale=T,
          total=F)
#>                       0.5         1.5        2.5        3.5        4.5
#> start:A_0      0.02139860 0.003165010 0.01716991 0.01512523 0.02287011
#> start:B_0      0.00000000 0.014828691 0.01124049 0.01858749 0.02150094
#> AVERAGE        0.01057919 0.009062318 0.01417192 0.01687579 0.02217784
#> AVERAGE(COND.) 0.00000000 0.014828691 0.01124049 0.01858749 0.02150094
#>                       5.5        6.5        7.5        8.5        9.5
#> start:A_0      0.02701001 0.03466032 0.04220125 0.05132105 0.06075695
#> start:B_0      0.02826202 0.03470815 0.04301741 0.05193800 0.06166635
#> AVERAGE        0.02764304 0.03468450 0.04261391 0.05163299 0.06121675
#> AVERAGE(COND.) 0.02826202 0.03470815 0.04301741 0.05193800 0.06166635
#>                      10.5       11.5       12.5       13.5       14.5
#> start:A_0      0.07033993 0.07894684 0.08551388 0.08873791 0.08758683
#> start:B_0      0.07130922 0.08007239 0.08671669 0.08999288 0.08882282
#> AVERAGE        0.07083002 0.07951593 0.08612203 0.08937244 0.08821176
#> AVERAGE(COND.) 0.07130922 0.08007239 0.08671669 0.08999288 0.08882282
#>                      15.5       16.5       17.5       18.5       19.5
#> start:A_0      0.08163513 0.07156045 0.05947112 0.04907794 0.03145157
#> start:B_0      0.08278814 0.07257081 0.06031091 0.04977093 0.03189568
#> AVERAGE        0.08221811 0.07207130 0.05989573 0.04942832 0.03167611
#> AVERAGE(COND.) 0.08278814 0.07257081 0.06031091 0.04977093 0.03189568
```

## Example 2: Simulated HRS data

Here we will provide a second example based on simulated data from the
Health and Retirement Study.

``` r
## Load package
library(dtms)

## Define model: Absorbing and transient states, time scale
hrs <- dtms(transient=c("Employed","Inactive","Retired"),
            absorbing="Dead",
            timescale=50:99)

## Quick look at data
head(hrsdata)
#> # A tibble: 6 × 4
#>      ID Gender   Age State
#>   <dbl>  <dbl> <dbl> <chr>
#> 1     1      1    50 <NA> 
#> 2     1      1    51 <NA> 
#> 3     1      1    52 <NA> 
#> 4     1      1    53 <NA> 
#> 5     1      1    54 <NA> 
#> 6     1      1    55 <NA>
```
