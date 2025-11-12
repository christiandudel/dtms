
# dtms - An R package for discrete-time multistate models

<!-- badges: start -->

[![R-CMD-check](https://github.com/christiandudel/dtms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/christiandudel/dtms/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Authors

Christian Dudel, <dudel@demogr.mpg.de>

## Disclaimer

This package is currently undergoing development and many functions are
experimental. The package comes with no warranty. The content of this
repository will change in the future, and functions and features might
be removed or changed without warning.

## Acknowledgements

We thank Peng Li, Alessandro Feraldi, Aapo Hiilamo, Daniel Schneider,
Donata Stonkute, Marcus Ebeling, and Angelo Lorenti for helpful
comments, suggestions, and code snippets. All errors remain our own.

## Citation

If you use this package in your work, please use the following citation
(or a variation):

Dudel, C. (2025). dtms: discrete-time multistate models in R. R package
version 0.3.5, available at <https://github.com/christiandudel/dtms>

## Overview

The package `dtms` implements discrete-time multistate models in R. It
comes with many tools to analyze the results of multistate models. The
workflow mainly consists of estimating a discrete-time multistate model
and then applying methods for absorbing Markov chains.

Currently, the following features are implemented:

- Data handling: functions for reshaping and aggregating data, cleaning
  data, editing states, generating indicators of duration and number of
  occurrences of a state, indicators of censoring, descriptive
  information on different types of censoring, and other general
  descriptive statistics.
- Estimation of transition probabilities: nonparametric estimation;
  semiparametric estimation
  ([VGAM](https://cran.r-project.org/web/packages/VGAM)), random effects
  and random intercepts
  ([mclogit](https://cran.r-project.org/web/packages/mclogit)), and
  neural networks
  ([nnet](https://cran.r-project.org/web/packages/nnet)); all possible
  for constrained and unconstrained/fully interacted models. Functions
  for descriptive statistics on transition probabilities and for
  plotting them are also available.
- Markov chain methods: survivorship function, (partial) state/life
  expectancy, (partial) lifetime risk, (partial) distribution of
  occupation time, (partial) distribution of waiting time to first
  visit, (partial) distribution of waiting time to last exit, (partial)
  distribution of waiting time to absorption, based on (partial)
  distributions variance/standard deviation and median of occupancy
  time/first visit/last exit/time to absorption, Markov chains with
  rewards.
- Inference: analytic standard errors and variance-covariance matrix for
  transition probabilities; simulated inference using the bootstrap and
  the block bootstrap for other quantities, supporting parallel
  computing.
- Other features: simulation of Markov chains using the package
  [markovchain](https://cran.r-project.org/web/packages/markovchain);
  survey weights (experimental); irregular time intervals
  (experimental).
- Examples: the package comes with two simulated data sets which are
  used for examples. These are described further below. The input data
  and code for the simulations is available at
  <https://github.com/christiandudel/dtms_data/>.

The documentation provided below does currently not cover all features.

## Content

Currently, the following topics are covered in this documentation

- Installation
- General workflow and basic principles
  - Model setup
  - Preparing and handling data
  - Estimating transition probabilities
  - Markov chain methods
- Example 1: artificial data
  - Data description
  - Model setup
  - Preparing and handling data
  - Estimating transition probabilities
  - Markov chain methods
- Example 2: simulated working trajectories
  - Data description
  - Analysis
  - Variance estimation
- Irregular intervals
- Combining dtms with other software
- Using dtms with secure data environments
- References

## Installation

You can install the development version of `dtms` like this:

``` r
install.packages("remotes")
remotes::install_github("christiandudel/dtms")
```

## General workflow and basic principles

The basic workflow consists of four main steps. First, the multistate
model is defined in a general way which describes the states included in
the model and its timescale (model setup). Second, the input data has to
be reshaped and cleaned. Third, a regression model is fitted to predict
transition probabilities. Fourth, Markov chain methods are applied to
calculate statistics to describe the model. These steps and the
corresponding functions are described below. Note that not all arguments
of each function are described, and the help files for the individual
functions usually contain useful additional information.

### Model setup

To use the `dtms` package, in a first step disrete-time multistate
models are defined in a abstract way using three components: the names
of the transient states; the names of the absorbing states; and the
values of the time scale. Moreover, there are two additional components
which the user not necessarily needs to specify: the step length of the
timescale, and a separator.

To define these components, the function `dtms()` is used. It has an
argument for each of the components, but only three are necessary: the
names of the transient states, the names of the absorbing states, and
the values of the time scale. The step length of the timescale is
implicitly defined by the values of the timescale, and the separator
uses a default value which users likely don’t want to change in a
majority of applications. In the first example provided further below,
the function `dtms()` is called like this:

``` r
## Load package
library(dtms)
## Define model: Absorbing and transient states, time scale
simple <- dtms(transient=c("A","B"),
               absorbing="X",
               timescale=0:19)
```

The arguments `transient` and `absorbing` take the names of the
transient and absorbing states, respectively, which are specified as
character vectors. In this case, there are two transient states called
`A` and `B`, and one absorbing state called `X`. Each model needs at
least one transient state and one absorbing state. The argument
`timescale` takes the values of the timescale which are specified with a
numeric vector. In this example, the time scale starts at 0 and stops at
19, with a step length of 1. The step length has to be consistent along
the timescale. For instance, a timescale with values 0, 1, 2, 4, 5, 7 is
not allowed. However, the step length does not need to equal 1; e.g., 0,
2, 4, 6, … would be fine. Moreover, using the argument `timestep`,
several values for the step length can be specified.

The separator is a character string used to construct what we call long
state names. Its default is `_`. Long state names consist of a
combination of the names of the transient states with values of the time
scale. They are used internally to handle that transition probabilities
might depend on values of the time scale. Long state names are never
constructed for absorbing states. For instance, if the transient states
are called `A` and `B`, the time scale can take on values 0, 1, and 2,
and the separator is `_`, then the following long state names will be
used: `A_0`, `A_1`, `A_2`, `B_0`, `B_1`, and `B_2`. Due to the temporal
ordering of states not all transitions between these states are
possible; e.g., it is not possible to transition from `A_2` to `B_0`.

### Preparing and handling data

The input data has to be panel data in long format. If your data is not
in this shape, there are many tools already available in R and its
extensions which allow you to reshape it. An example of data in long
format could look like this:

| idvar | timevar | statevar | X   | Y    |
|:------|:--------|:---------|:----|:-----|
| 1     | 0       | A        | 2   | 1020 |
| 1     | 1       | A        | 2   | 1025 |
| 1     | 2       | B        | 2   | 1015 |
| 1     | 3       | A        | 2   | 1000 |
| 2     | 0       | B        | 1   | 2300 |
| 2     | 1       | A        | 1   | 2321 |
| …     | …       | …        | …   | …    |

The first variable, `idvar`, contains a unit identifier. The first four
rows of the data belong to unit `1`. The variable `timevar` has the
values of the timescale. `statevar` shows the state each unit is
occupying at a given time. Ideally, the states are provided as character
strings; numeric values will also work, factors are, however, currently
not supported. `X` and `Z` are additional covariates.

The `dtms` package provides tools to reshape this data into what we call
transition format. For the example data shown above the transition
format looks like this:

| idvar | timevar | fromvar | tovar | X   | Y    |
|:------|:--------|:--------|:------|:----|:-----|
| 1     | 0       | A       | A     | 2   | 1020 |
| 1     | 1       | A       | B     | 2   | 1025 |
| 1     | 2       | B       | A     | 2   | 1015 |
| 2     | 0       | B       | A     | 1   | 2300 |
| …     | …       | …       | …     | …   | …    |

Each row shows for each unit (`idvar`) and given time (`timevar`) the
state currently occupied (`fromvar`) and the state the unit will
transition to at the next value of the time scale (`tovar`). For unit 1,
the last observation in long format is at time 3. However, this is the
final observation and there is no transition to another state after
this. This means that the last observed transition for unit 1 starts at
time 2.

To reshape data into transition format, the function `dtms_format()` can
be used. It is one of the many functions of the package which takes the
result of the function `dtms()` as one of its inputs. In the second
example provided below `dtms_format()` is used like this:

``` r
## Load package
library(dtms)
## Define model: Absorbing and transient states, time scale
work <- dtms(transient=c("Working","Non-working","Retired"),
             absorbing="Dead",
             timescale=50:99)
## Reshape
estdata <- dtms_format(data=workdata,
                       dtms=work,
                       idvar="ID",
                       timevar="Age",
                       statevar="State")
```

First, `dtms()`is called to define the model. The call of
`dtms_format()` specifies the data frame which contains the data in long
format (argument `data`), the definition of the multistate model
(argument `dtms`), the name of the variable containing the unit
identifier (argument `idvar`), the name of the variable containing the
values of the timescale (argument `timevar`), and the variable
indicating the states (argument `statevar`).

Note that at this stage, the states captured by the variable specified
with `statevar` do not necessarily need to match the states in the
`dtms` object. In particular, this means that if there are any states in
the input data which are not in the `dtms` object, there will be no
warning or similar at this stage. However, the function `dtms_clean()`
described below will remove states not included in the `dtms` object by
default.

The original data and the reshaped data can be compared like this:

``` r
workdata |> subset(ID==3) |> head()
#>     ID Gender Age       State
#> 101  3      1  50     Working
#> 102  3      1  51     Working
#> 103  3      1  52     Working
#> 104  3      1  53     Working
#> 105  3      1  54 Non-working
#> 106  3      1  55 Non-working
estdata |> subset(id==3) |> head()
#>     id Gender time        from          to
#> 101  3      1   50     Working     Working
#> 102  3      1   51     Working     Working
#> 103  3      1   52     Working     Working
#> 104  3      1   53     Working Non-working
#> 105  3      1   54 Non-working Non-working
#> 106  3      1   55 Non-working Non-working
```

Three things are important to note. First, if not specified otherwise,
`dtms_format()` changes variable names to default names. These default
names are also default values for other functions, meaning that certain
variable names do not need to be specified all the time, making the
workflow easier. Specifically, the variable with the unit identifier
gets the name `id`; the variable with the timescale gets the name
`time`; and the variables with the starting and receiving state get the
names `from` and `to`. Other names can of course be specified, but they
have to be used consistently.

Second, the object returned by `dtms_format()` technically is a standard
data frame. This comes with benefits and costs. On the upside this means
that data in transition format can easily be viewed and modified using
standard tools, making it very accessible to users. The main downside is
that it does not contain general information on the model or the data,
which could make the workflow slightly more convenient. We decided to
keep intermediate steps as accessible as possible.

Third, the data in long format contains an additional variable
(`Gender`). All variables which are not `idvar`, `timevar`, or
`statevar` do not need to be specified and are handled in the same way.
For any variables X, Y, Z, … the value at time $t$ is assigned to the
transition starting in time $t$.

A useful function for handling data in transition format is
`dtms_clean()`. It can be used to remove, for instance, transitions
starting and/or ending in a missing value. Depending on the data, these
can occur quite frequently. For instance, in the example used above:

``` r
estdata |> subset(id==1) |> head()
#>   id Gender time from   to
#> 1  1      1   50 <NA> <NA>
#> 2  1      1   51 <NA> <NA>
#> 3  1      1   52 <NA> <NA>
#> 4  1      1   53 <NA> <NA>
#> 5  1      1   54 <NA> <NA>
#> 6  1      1   55 <NA> <NA>
```

This means that for the unit with ID 1, at least the first couple of
transitions only contain missing values. This is because in the original
data these values are missing:

``` r
workdata |> subset(ID==1) |> head()
#>   ID Gender Age State
#> 1  1      1  50  <NA>
#> 2  1      1  51  <NA>
#> 3  1      1  52  <NA>
#> 4  1      1  53  <NA>
#> 5  1      1  54  <NA>
#> 6  1      1  55  <NA>
```

Technically, this is because `dtms_format()` takes each row in the input
data and adds the state at the next value of the time scale to that row.
If this next value is NA, or if it is not included in the data, the
resulting value of \`to’ will be NA.

Such missing values and a few other things can be removed as follows:

``` r
estdata <- dtms_clean(data=estdata,
                      dtms=work)
#> Dropping  0  rows not in state space
#> Dropping  0  rows not in time range
#> Dropping  81903  rows starting or ending in NA
#> Dropping  68319  rows starting in absorbing state
```

`dtms_clean()` by default removes transitions starting or ending in a
missing value; transitions starting in an absorbing state; transitions
which start at a time not covered by the time scale; and transitions
which start or end in a state which is not contained in the state space
in the `dtms` object. This can be changed through some of the arguments
of `dtms_clean()`.

### Estimating transition probabilities

Getting transition probabilities ready requires three steps. First,
estimating a regression model for the transition probabilities. Second,
predicting transition probabilities using this model. Third, putting the
transition probabilities into a transition matrix.

Estimating a regression model is done using `dtms_fit()`. Currently,
this function builds on several other packages to allow for, among other
things, semiparametric estimation and random effects. In the first
example below, we show a very basic call of this function:

``` r
## Fit model 
fit <- dtms_fit(data=estdata)
```

This will estimate a model which only uses the starting state as a
control. Covariates (including time) can be included in several ways. A
convenient way for including covariates is to use the argument
`controls`. For example, if in a data set there are two variables named
`Z1` and `Z2`, and time is captured with a variable called `time`, then
these variables could be included as follows:

``` r
fit <- dtms_fit(data=somedata,
                controls=c("time","Z1","Z2"))
```

It is also possible to specify a formula, like is common practice for
most regression functions. This way you have to specify the names of the
state variables:

``` r
fit <- dtms_fit(data=somedata,
                formula=to~from+time+Z1+Z2)
```

What class of object is generated by `dtms_fit()` depends on the package
which is used for estimation. This is controlled by the argument
`package`, with `VGAM` as default and then using the function `vgam()`.
This means that the resulting objects of the short examples above have
class `vgam` and behave accordingly, in particular with methods for
functions like `summary()` or `coef()`. In addition to `VGAM`, currently
`mclogit` and `nnet` are supported. If the package is set to `mclogit`
the function `mblogit()` is used. If it is `nnet` the function is
`multinom()`. Arguments of these functions can be passed via
`dtms_fit()`. For instance, this example includes random effects by unit
ID using the `mclogit` package and the argument `random` of `mblogit()`:

``` r
fit <- dtms_fit(data=somedata,
                covariates=c("time","Z1","Z2"),
                package="mclogit",
                random=~1|id)
```

Once a regression model has been estimated, it can be used to predict
transition probabilities. If the model includes covariates, the user
needs to specify covariate values which are used in the prediction. If
no covariates are included, then no values are needed. For each
time-constant covariate one value needs to be specified, and for
time-varying variables a value at each value of the timescale has to be
provided minus the last value. For instance, if the timescale has values
0, 1, 2, and 3, and there is a time-varying covariate X, then values for
X at 0, 1, and 2 need to be specified. Alternatively, predictor values
for all values of the time scale can be specified. In this case, the
last predictor value will be dropped. For a time-constant covariate Y
only one values is necessary, which is then used at all values of the
time scale.

To predict transition probabilities, the function `dtms_transitions()`
is used. It has three main arguments: `model`, `dtms`, and `controls`.
The first argument is used to specify the object which contains the
regression model. The argument `dtms` takes the definition of the
multistate model as generated by `dtms()`. The argument `controls` is
only needed if there are covariates in the model. It take as lists with
named entries, where the names need to correspond to the names of the
covariates in the model. For instance, in a model with a timescale with
values 0, 1, 2, and 3, and a time-constant covariate Y and a
time-varying covariate X, the call of `dtms_transitions()` could look as
follows:

``` r
probs <- dtms_transitions(model=fit,
                          dtms=example,
                          controls=list(time=0:3,
                                        Y=1,
                                        X=c(500,501,500,503))) 
```

The result of calling `dtms_transitions()` is a data frame with
predicted transition probabilities, where each row contains the
transition probability for one specific transition, as indicated by the
variables in the data. This object will often be only an intermediate
step. However, the package provides several tools to look at the
transition probabilities (see the examples below), and putting the
probabilities in a data frame makes them easily accessible to the user.

Transition probabilities in a data frame have to be converted to a
transition matrix before they can be used further. For this, we use the
function `dtms_matrix()`. In most cases, this function will only require
two arguments. First, a data frame with transition probabilities as
created with `dtms_transitions()` and passed to the argument `probs`;
and, second, a `dtms` object as created with `dtms()`. For instance, in
the first example below, this function is used as follows:

``` r
Tp <- dtms_matrix(probs=probs,
                  dtms=simple)
```

This transition matrix can then be used to apply Markov chain methods.
The matrix itself often will not be of major interest, and it is easier
to look at transition probabilities in the data frame generated with
`dtms_transitions()`.

### Markov chain methods

The dtms package provides several functions which implement Markov chain
methods and which can be applied to a transition matrix generated with
`dtms_matrix()`. Most of these functions require at least two arguments.
First, a transition matrix; and, second, a `dtms` object. For instance,
to calculate the lifetime spent in the different states, the function
`dtms_expectancy()`. This and other functions are demonstrated in more
detail in the two examples below.

## Example 1: Artificial data

### Data description

This is a basic example using artificial data which is provided with the
package. The state space consists of two transient states (A, B), and
one absorbing state (X). The time scale goes from 0 to 20. Transition
probabilities do change depending on time, as we will see below.

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
#> 4  1    3     B
#> 5  1    4     A
#> 6  1    5     B

## States
simpledata$state |> unique()
#> [1] "A" "B" "X"

## Number of units
simpledata$id |> unique() |> length()
#> [1] 993

## Number of observations
dim(simpledata)
#> [1] 12179     3
```

The data set is in long format and contains three variables. ‘id’ is an
unit identifier; ‘time’ contains the value of the time scale; and
‘state’ contains the state the unit occupied at a given time. In total,
there are 993 units, each of them contributing to the total of 12,173
observations.

### Model setup

To work with this data set, we first define a basic discrete-time
multistate model. This is done with the function `dtms()` and requires
us to specify the names of the transient states, the names of the
absorbing states, and the possible values of the time scale:

``` r
## Define model: Absorbing and transient states, time scale
simple <- dtms(transient=c("A","B"),
               absorbing="X",
               timescale=0:19)
```

The resulting object of class `dtms` will be passed to other functions
of the package.

### Data handling

In a second step, we transform the data from long format to transition
format. This can be done using the function `dtms_format()`. In this
example, we need to specify the name of the object containing the data,
and in addition a `dtms` object as created above. Moreover, we need to
specify which variables contain the unit identifier, which variable
contains the values of the time scale, and which variable contains the
information on the state:

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
#> 3  1    2    B  B
#> 4  1    3    B  A
#> 5  1    4    A  B
#> 6  1    5    B  A
```

While in long format each row contains information on the currently
occupied state, in transition format also the next state is shown in a
variable. If there is no observation at time t+1, the next state is NA.
The names of the variables in the resulting data set are by default
chosen such that they match defaults of other functions of the package.

Depending on the original data, there can be missing values in
transition format data due to several reasons. The function
`dtms_clean()` provides a convenient way to remove such rows in the
data, as well as other potentially problematic or unwanted rows. It
returns a cleaned data set and prints a brief overview of the dropped
rows to the console:

``` r
## Missing values?
estdata$to |> table(useNA="always")
#> 
#>    A    B    X <NA> 
#> 3625 6664  627 1263

## Clean
estdata <- dtms_clean(data=estdata,
                      dtms=simple)
#> Dropping  0  rows not in state space
#> Dropping  0  rows not in time range
#> Dropping  1263  rows starting or ending in NA
#> Dropping  0  rows starting in absorbing state
```

In this example, 1,260 transitions were dropped because they end in a
missing value. No observations were dropped because the states are not
covered by the state space; no observations are dropped because they are
out of the time range specified with the ‘dtms’ object; and no
observations were dropped because they start in an absorbing state.

A brief overview of the data is provided when using the function
`summary()`:

``` r
## Summary of data
summary(estdata)
#>   from to COUNT       PROP       PROB
#> 1    A  A   372 0.03407842 0.09451220
#> 2    A  B  3413 0.31266032 0.86712398
#> 3    A  X   151 0.01383291 0.03836382
#> 4    B  A  3253 0.29800293 0.46604585
#> 5    B  B  3251 0.29781971 0.46575931
#> 6    B  X   476 0.04360572 0.06819484
```

This shows for all possible transitions the absolute number each
transition is observed (e.g., there are 160 transitions from A to X);
the proportion of each transition relative to all transitions (e.g., a
bit more than 1% of all observed transitions are from A to X); and raw
transition probabilities (e.g., the probability of transitioning to X
starting in A is around 4%).

Some more information on the data is provided by the function
`dtms_censoring()`. It can be used in different ways, but a basic
version shows an overview of the number of units with left censoring,
the number of units with gaps in their series of observations, and the
number of units with right censoring:

``` r
dtms_censoring(data=estdata,
               dtms=simple)
#> Units with left censoring:  242 
#> Units with gaps:  196 
#> Units with right censoring:  325
```

### Estimating transition probabilities

To estimate the transition probabilities of the multistate model, the
function `dtms_fit()` is used. In this simple example, it is sufficient
to specify the name of the object with the transition data and the time
scale as a control variable:

``` r
## Fit model 
fit <- dtms_fit(data=estdata,
                controls="time")
```

To predict transition probabilities and to arrange them in a matrix, the
functions `dtms_transitions()` and `dtms_matrix()` are used. The
function `dtms_transitions()` needs a ‘dtms’ object as well as a fitted
model and values for the control variable, while the function
`dtms_matrix()` requires a `dtms` object and predicted probabilities:

``` r
## Predict probabilities
probs    <- dtms_transitions(dtms=simple,
                             model = fit,
                             controls=list(time=simple$timescale))

## Get transition matrix 
Tp <- dtms_matrix(dtms=simple,
                  probs=probs)
```

In more complex examples, the previous functions would need more
information. For instance, on which covariates to include in the
estimation step, and which covariate values to use in the prediction
step.

To get an overview of the transition probabilities, the function summary
can be used:

``` r
## Summary of probabilities
summary(probs)
#>   from to    MIN MINtime    MAX MAXtime MEDIAN   MEAN
#> 1    A  A 0.0775      18 0.0982       0 0.0946 0.0922
#> 3    A  B 0.7179      18 0.8964       0 0.8695 0.8473
#> 5    A  X 0.0054       0 0.2045      18 0.0359 0.0605
#> 2    B  A 0.3394      18 0.4967       0 0.4682 0.4497
#> 4    B  B 0.3422      18 0.4936       0 0.4686 0.4500
#> 6    B  X 0.0096       0 0.3184      18 0.0632 0.1003
```

For all combinations of starting and receiving state, this shows the
lowest transition probability and at what value of the time scale it
occurs. It also shows the same for the highest transition probability,
and in addition it shows the median and the mean of all transition
probabilities between two states.

Another useful way to look at the transition probabilities is to plot
them. To make this easy, the package provides the function
`dtms_simplify()` which can be applied to an object created with
`dtms_transitions()` to make it easier to plot. For instance, using
ggplot2, a simple plot could look like this:

``` r
## Simple plot
library(ggplot2)
probs |>  dtms_simplify() |> 
          ggplot(aes(x=time,y=P,color=to)) + 
          geom_line() + 
          facet_wrap(~from)
```

<img src="man/figures/README-example1-probsplot-1.png" width="100%" />
An even simpler way is available which builds on base-R and does not
require ggplot2. However, this creates less nice figures and is mainly
intended as a very quick way of checking results:

``` r
## Simple base plot
plot(probs,dtms=simple)
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
the functions used for calculating Markov chain methods, providing
additional information.

### Markov chain methods

Most functions used to calculate Markov chain methods need a transition
matrix and a `dtms` object, and potentially further arguments. The two
examples below calculate the expected time spent in a state
(`dtms_expectancy()`) and the lifetime risk of ever reaching a state
(`dtms_risk()`). In the first case, the starting distribution of states
is passed to the function; this is optional. In the second case, one or
several states need to be specified for which the lifetime risk will be
calculated:

``` r
## State expectancies 
dtms_expectancy(dtms=simple,
                matrix=Tp,
                start_distr=S)
#>                  A        B    TOTAL
#> start:A_0 4.989071 8.686304 13.67538
#> start:B_0 4.761130 8.875429 13.63656
#> AVERAGE   4.876064 8.780067 13.65613

## Lifetime risk 
dtms_risk(dtms=simple,
          matrix=Tp,
          risk="A")
#>       A_0       B_0 
#> 1.0000000 0.9754457
```

The results of the call of `dtms_expectancy()` show the starting states
in rows and the states in which the time is spent as columns. For
instance, around 5.05 time units are spent in state A when starting in
state A at time 0. The last column shows the total time until
absorption. The last row is shown because the starting distribution was
specified. It shows the average time spent in a state irrespective of
the starting state. That is, on average 4.93 time units are spent in
state A, and 8.88 time units are spent in state B, for a total of 13.81
time units.

The result of the call of `dtms_risk()` above is the lifetime risk of
ever reaching state A depending on the starting state. Obviously, when
starting in state A at time 0, this risk amounts to 1. When starting in
state B, the risk is also very high and around 97%. Note that for
consistent estimation the Markov assumption has to hold. `dtms_risk()`
can also be used such that the Markov assumption is not necessary, but
this requires additional data editing and use of the function
`dtms_forward()`. In particular, the latter function can be used to
create an absorbing set which contains the state of interest plus all
absorbing states. This has to be done in the data editing process, and
all following editing and estimation steps have to be repeated:

``` r
riskdata <- dtms_forward(data=estdata,
                         state="A",
                         dtms=simple)

riskfit <- dtms_fit(data=riskdata,
                controls="time",
                package="mclogit")
#> 
#> Iteration 1 - deviance = 6381.024 - criterion = 0.7477243
#> Iteration 2 - deviance = 5222.033 - criterion = 0.2219381
#> Iteration 3 - deviance = 4899.832 - criterion = 0.06575636
#> Iteration 4 - deviance = 4804.59 - criterion = 0.01982266
#> Iteration 5 - deviance = 4771.139 - criterion = 0.007010997
#> Iteration 6 - deviance = 4758.893 - criterion = 0.00257315
#> Iteration 7 - deviance = 4754.393 - criterion = 0.0009464568
#> Iteration 8 - deviance = 4752.738 - criterion = 0.000348183
#> Iteration 9 - deviance = 4752.13 - criterion = 0.0001280896
#> Iteration 10 - deviance = 4751.906 - criterion = 4.712157e-05
#> Iteration 11 - deviance = 4751.823 - criterion = 1.733506e-05
#> Iteration 12 - deviance = 4751.793 - criterion = 6.377213e-06
#> Iteration 13 - deviance = 4751.782 - criterion = 2.346046e-06
#> Iteration 14 - deviance = 4751.778 - criterion = 8.63062e-07
#> Iteration 15 - deviance = 4751.776 - criterion = 3.175028e-07
#> Iteration 16 - deviance = 4751.776 - criterion = 1.168027e-07
#> Iteration 17 - deviance = 4751.776 - criterion = 4.296933e-08
#> Iteration 18 - deviance = 4751.776 - criterion = 1.580753e-08
#> Iteration 19 - deviance = 4751.776 - criterion = 5.815266e-09
#> converged

riskprobs <- dtms_transitions(dtms=simple,
                          controls=list(time=simple$timescale),
                          model = riskfit)

riskTp <- dtms_matrix(dtms=simple,
                  probs=riskprobs)

dtms_risk(dtms=simple,
          matrix=riskTp,
          risk="A")
#>       A_0       B_0 
#> 1.0000000 0.9742659
```

In this example, the results of the naive use of `dtms_risk()` and the
more elaborate approach using `dtms_forward()` lead to very similar
results.

It is also possible to calculate state expectancies conditional on
values of the time scale. For this, a single transient state has to be
specified for which we want to know how long units spent in it:

``` r
dtms_expectancy(dtms=simple,
                risk="A",
                matrix=Tp)
#>                  0        1        2        3        4        5        6
#> start:A_0 4.989071 4.671601 4.360438 4.056333 3.760074 3.472470 3.194327
#> start:B_0 4.761130 4.441133 4.127109 3.819781 3.519910 3.228284 2.945695
#>                  7        8        9       10       11       12       13
#> start:A_0 2.926429 2.669490 2.424128 2.190751 1.969608 1.760326 1.562651
#> start:B_0 2.672918 2.410677 2.159596 1.920167 1.692603 1.476925 1.272314
#>                 14        15        16        17        18  19
#> start:A_0 1.373716 1.1943763 1.0091560 0.8407368 0.5775494 0.5
#> start:B_0 1.078217 0.8902818 0.7100656 0.5136192 0.3394079 0.0
```

In the example above, we look at the time spent in state A. The results
by row now indicate the remaining life expectancy in state A starting
from the state in the row at a given time. For instance, if a unit is in
state A at time 5, an additional 3.51 time units will be spent in state
A (as seen in row named “A”, column named “5”).

The function calls below are all similar in that they provide full
distributions as a result. Specifically, `dtms_visits()` calculates the
distribution of the time spent in a state; the mean over this
distribution is equal to the state expectancy as provided by
`dtms_expectancy()`, and one minus the proportion of 0 time units spent
in a state is equal to the lifetime risk provided by `dtms_risk()`.
`dtms_first()` calculates the distribution of the waiting time until a
given state is reached for the first time, conditional on ever reaching
this state. `dtms_last()` calculates the distribution of the waiting
time until a state is left for the last time; i.e., there is no return
back to this state.

``` r
## Distribution of visits
dtms_visits(dtms=simple,
            matrix=Tp,
            risk="A",
            start_distr=S)
#>                         0           0.5             1        1.5          2
#> start:A_0      0.00000000  3.243529e-02 -6.938894e-18 0.05723796 0.00000000
#> start:B_0      0.02455429 -3.469447e-18  4.915136e-02 0.00000000 0.08175976
#> AVERAGE        0.01217339  1.635470e-02  2.436800e-02 0.02886083 0.04053442
#> AVERAGE(COND.) 0.02455429 -3.469447e-18  4.915136e-02 0.00000000 0.08175976
#>                       2.5          3           3.5             4           4.5
#> start:A_0      0.09305499 0.00000000  1.348112e-01 -1.110223e-16  1.692934e-01
#> start:B_0      0.00000000 0.12233817 -5.551115e-17  1.599971e-01 -5.551115e-17
#> AVERAGE        0.04692069 0.06065216  6.797523e-02  7.932250e-02  8.536202e-02
#> AVERAGE(COND.) 0.00000000 0.12233817 -5.551115e-17  1.599971e-01 -5.551115e-17
#>                         5        5.5             6           6.5             7
#> start:A_0      0.00000000 0.18020372 -1.110223e-16  1.576207e-01 -1.110223e-16
#> start:B_0      0.17875845 0.00000000  1.661985e-01 -1.110223e-16  1.224907e-01
#> AVERAGE        0.08862391 0.09086329  8.239698e-02  7.947634e-02  6.072777e-02
#> AVERAGE(COND.) 0.17875845 0.00000000  1.661985e-01 -1.110223e-16  1.224907e-01
#>                      7.5          8        8.5          9         9.5
#> start:A_0      0.1066097 0.00000000 0.05067043 0.00000000 0.015169236
#> start:B_0      0.0000000 0.06591306 0.00000000 0.02331358 0.000000000
#> AVERAGE        0.0537553 0.03267802 0.02554931 0.01155828 0.007648713
#> AVERAGE(COND.) 0.0000000 0.06591306 0.00000000 0.02331358 0.000000000
#>                         10        10.5           11         11.5           12
#> start:A_0      0.000000000 0.002623795 0.0000000000 0.0002548551 0.000000e+00
#> start:B_0      0.004905569 0.000000000 0.0005794685 0.0000000000 3.852707e-05
#> AVERAGE        0.002432057 0.001322984 0.0002872858 0.0001285044 1.910074e-05
#> AVERAGE(COND.) 0.004905569 0.000000000 0.0005794685 0.0000000000 3.852707e-05
#>                        12.5           13         13.5           14         14.5
#> start:A_0      1.425497e-05 0.000000e+00 4.792622e-07 0.000000e+00 1.003948e-08
#> start:B_0      0.000000e+00 1.496713e-06 0.000000e+00 3.532144e-08 0.000000e+00
#> AVERAGE        7.187719e-06 7.420324e-07 2.416562e-07 1.751147e-08 5.062159e-09
#> AVERAGE(COND.) 0.000000e+00 1.496713e-06 0.000000e+00 3.532144e-08 0.000000e+00
#>                          15          15.5            16          16.5
#> start:A_0      0.000000e+00  1.335809e-10 -2.220446e-16  1.125877e-12
#> start:B_0      5.187080e-10 -1.110223e-16  4.748202e-12 -1.110223e-16
#> AVERAGE        2.571623e-10  6.735483e-11  2.353926e-12  5.676408e-13
#> AVERAGE(COND.) 5.187080e-10 -1.110223e-16  4.748202e-12 -1.110223e-16
#>                          17         17.5 18 18.5 19 19.5 20          20.5
#> start:A_0      0.000000e+00 5.884182e-15  0    0  0    0  0 -4.440892e-16
#> start:B_0      2.642331e-14 0.000000e+00  0    0  0    0  0  0.000000e+00
#> AVERAGE        1.310001e-14 2.966954e-15  0    0  0    0  0 -2.239210e-16
#> AVERAGE(COND.) 2.642331e-14 0.000000e+00  0    0  0    0  0  0.000000e+00
#> attr(,"class")
#> [1] "dtms_distr" "matrix"

## Distribution of waiting time to first visit
dtms_first(dtms=simple,
           matrix=Tp,
           risk="A",
           start_distr=S)
#>                        0       0.5       1.5       2.5        3.5        4.5
#> start:A_0      1.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.00000000
#> start:B_0      0.0000000 0.5092282 0.2506884 0.1231091 0.06026808 0.02938768
#> AVERAGE        0.5104391 0.2492982 0.1227272 0.0602694 0.02950490 0.01438706
#> AVERAGE(COND.) 0.0000000 0.5092282 0.2506884 0.1231091 0.06026808 0.02938768
#>                        5.5         6.5         7.5          8.5          9.5
#> start:A_0      0.000000000 0.000000000 0.000000000 0.0000000000 0.0000000000
#> start:B_0      0.014258726 0.006875277 0.003289492 0.0015587836 0.0007299212
#> AVERAGE        0.006980515 0.003365867 0.001610407 0.0007631195 0.0003573409
#> AVERAGE(COND.) 0.014258726 0.006875277 0.003289492 0.0015587836 0.0007299212
#>                        10.5         11.5         12.5         13.5         14.5
#> start:A_0      0.0000000000 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00
#> start:B_0      0.0003368277 1.526681e-04 6.769946e-05 2.923404e-05 1.222575e-05
#> AVERAGE        0.0001648977 7.474034e-05 3.314301e-05 1.431184e-05 5.985251e-06
#> AVERAGE(COND.) 0.0003368277 1.526681e-04 6.769946e-05 2.923404e-05 1.222575e-05
#>                        15.5         16.5        17.5         18.5
#> start:A_0      0.000000e+00 0.000000e+00 0.00000e+00 0.000000e+00
#> start:B_0      4.920140e-06 1.891560e-06 6.88995e-07 2.356016e-07
#> AVERAGE        2.408708e-06 9.260336e-07 3.37305e-07 1.153413e-07
#> AVERAGE(COND.) 4.920140e-06 1.891560e-06 6.88995e-07 2.356016e-07
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
#> start:A_0      0.03243583 0.003913531 0.02214835 0.01825495 0.02599464
#> start:B_0      0.00000000 0.020285002 0.01467573 0.02253729 0.02469664
#> AVERAGE        0.01655652 0.011928361 0.01849005 0.02035141 0.02535919
#> AVERAGE(COND.) 0.00000000 0.020285002 0.01467573 0.02253729 0.02469664
#>                       5.5        6.5        7.5        8.5        9.5
#> start:A_0      0.02898964 0.03517085 0.04067050 0.04718290 0.05364461
#> start:B_0      0.03053533 0.03552299 0.04177723 0.04813699 0.05488052
#> AVERAGE        0.02974635 0.03534324 0.04121231 0.04764999 0.05424966
#> AVERAGE(COND.) 0.03053533 0.03552299 0.04177723 0.04813699 0.05488052
#>                      10.5       11.5       12.5       13.5       14.5
#> start:A_0      0.06015945 0.06617171 0.07137684 0.07546877 0.07858666
#> start:B_0      0.06147681 0.06765129 0.07295944 0.07714784 0.08033266
#> AVERAGE        0.06080437 0.06689605 0.07215162 0.07629078 0.07944143
#> AVERAGE(COND.) 0.06147681 0.06765129 0.07295944 0.07714784 0.08033266
#>                      15.5       16.5       17.5       18.5
#> start:A_0      0.08172420 0.08733604 0.09696570 0.07380484
#> start:B_0      0.08354095 0.08927708 0.09912097 0.07544525
#> AVERAGE        0.08261361 0.08828629 0.09802083 0.07460792
#> AVERAGE(COND.) 0.08354095 0.08927708 0.09912097 0.07544525
#> attr(,"class")
#> [1] "dtms_distr" "matrix"
```

The output from these functions tends to be difficult to read, and often
results on the distribution are used to calculate other statistics. A
set of such statistics can be generated using the function `summary()`:

``` r
## Distribution of visits
example <- dtms_visits(dtms=simple,
                       matrix=Tp,
                       risk="A",
                       start_distr=S)
summary(example)
#>                    MEAN VARIANCE       SD MEDIAN      RISK0
#> start:A_0      4.989071 4.369975 2.090449    5.5 0.00000000
#> start:B_0      4.761130 4.496290 2.120446    5.0 0.02455429
#> AVERAGE        4.876064 4.445587 2.108456    5.0 0.01217339
#> AVERAGE(COND.) 4.761130 4.496290 2.120446    5.0 0.02455429
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
#> start:A_0      12.07534 23.17287 4.813821   12.5    NA
#> start:B_0      12.33995 20.47255 4.524660   13.5    NA
#> AVERAGE        12.20489 21.86840 4.676366   12.5    NA
#> AVERAGE(COND.) 12.33995 20.47255 4.524660   13.5    NA
```

In this case, the distribution is conditional on ever experiencing the
exit from state A, such that the waiting time until exit always has to
be above 0.

With respect to inference, the package currently provides analytic
standard errors for transition probabilities, and simulated inference
for all other statistics. Standard errors for transition probabilities
are provided by default by the function `dtms_transitions()`, and it can
also return confidence intervals. Two methods for simulated inference
are implemented: the (simple) bootstrap and the block bootstrap. Both
methods are provided by the function `dtms_boot()`. How to use
`dtms_transitions()` and `dtms_boot()` for inference is demonstrated in
the second example below.

## Example 2: Simulated working trajectories

### Data description

Here we provide an example using simulated working trajectories. The
simulations are are conducted using transition probabilities estimated
from the US Health and Retirement Study (HRS) and published by Dudel &
Myrskylä (2017) who studied working trajectories in late working life
and old age. These transition probabilities are used to simulate
artificial but realistic trajectories. There are three transient states
(working, non-working, retired) and one absorbing state (dead). The time
scale represents age and ranges from 50 to 99, as the focus is on older
individuals. Note that the actual HRS data is collected every two years
and while the simulated data is annual. The data set also contains each
individual’s gender, and the transition probabilities underlying the
simulated trajectories differ between men and women.

### Analysis

The workflow is similar to the previous example. First, a ‘dtms’ model
is defined using the function \`dtms’. Second, the data is brought into
transition format and cleaned. Third, transition probabilities are
estimated and put into a transition matrix. In this example,
probabilities are estimated and predicted using time-constant and
time-varying covariates, and the probabilities are plotted together with
confidence intervals. Finally, the transition matrix is used to
calculate state expectancies and similar measures.

``` r
## Load packages
library(dtms)
library(ggplot2)

## Define model: Absorbing and transient states, time scale
work <- dtms(transient=c("Working","Non-working","Retired"),
             absorbing="Dead",
             timescale=50:99)

## Quick look at data
head(workdata)
#>   ID Gender Age State
#> 1  1      1  50  <NA>
#> 2  1      1  51  <NA>
#> 3  1      1  52  <NA>
#> 4  1      1  53  <NA>
#> 5  1      1  54  <NA>
#> 6  1      1  55  <NA>

## Reshape
estdata <- dtms_format(data=workdata,
                       dtms=work,
                       idvar="ID",
                       timevar="Age",
                       statevar="State")

## Drop dead-to-dead transitions etc
estdata <- dtms_clean(data=estdata,
                      dtms=work)
#> Dropping  0  rows not in state space
#> Dropping  0  rows not in time range
#> Dropping  81903  rows starting or ending in NA
#> Dropping  68319  rows starting in absorbing state

## Overview
summary(estdata)
#>           from          to COUNT        PROP        PROB
#> 1  Non-working        Dead   197 0.001974383 0.013936050
#> 2  Non-working Non-working 10635 0.106586622 0.752334465
#> 3  Non-working     Retired  1900 0.019042274 0.134408602
#> 4  Non-working     Working  1404 0.014071238 0.099320883
#> 5      Retired        Dead  2602 0.026077893 0.051556401
#> 6      Retired Non-working   606 0.006073483 0.012007371
#> 7      Retired     Retired 46423 0.465262884 0.919831976
#> 8      Retired     Working   838 0.008398645 0.016604252
#> 9      Working        Dead   306 0.003066808 0.008699855
#> 10     Working Non-working  2066 0.020705967 0.058738237
#> 11     Working     Retired  2178 0.021828459 0.061922497
#> 12     Working     Working 30623 0.306911343 0.870639411

## Basic censoring
dtms_censoring(data=estdata,
               dtms=work)
#> Units with left censoring:  2036 
#> Units with gaps:  1720 
#> Units with right censoring:  1323

## More advanced censoring example
estdata <- dtms_censoring(data=estdata,
                          dtms=work,
                          add=T,
                          addtype="obs")
#> Units with left censoring:  2036 
#> Units with gaps:  1720 
#> Units with right censoring:  1323

estdata |>
  subset(subset=to!="Dead",select=c(RIGHT,to)) |>
  table() |>
  prop.table(margin=1)
#>        to
#> RIGHT   Non-working    Retired    Working
#>   FALSE  0.13846880 0.51950708 0.34202412
#>   TRUE   0.07860922 0.73015873 0.19123205

## Add age squared
estdata$time2 <- estdata$time^2
  
## Fit model
fit <- dtms_fit(data=estdata,
                controls=c("Gender","time","time2"))

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
                                            time2 =(50:98)^2),
                            CI=TRUE)

# Overview
summary(probs_m)
#>           from          to    MIN MINtime    MAX MAXtime MEDIAN   MEAN
#> 1  Non-working        Dead 0.0064      50 0.7701      98 0.1080 0.2055
#> 4  Non-working Non-working 0.0000      98 0.8173      54 0.1402 0.3186
#> 7  Non-working     Retired 0.0356      50 0.7089      79 0.3904 0.3880
#> 10 Non-working     Working 0.0238      98 0.1668      50 0.1008 0.0879
#> 2      Retired        Dead 0.0231      58 0.4588      98 0.0364 0.0930
#> 5      Retired Non-working 0.0000      98 0.1627      50 0.0025 0.0324
#> 8      Retired     Retired 0.5365      98 0.9492      73 0.8855 0.8381
#> 11     Retired     Working 0.0047      98 0.2132      50 0.0121 0.0365
#> 3      Working        Dead 0.0013      50 0.4725      98 0.0288 0.0926
#> 6      Working Non-working 0.0000      98 0.0704      58 0.0104 0.0252
#> 9      Working     Retired 0.0075      50 0.2443      85 0.1524 0.1375
#> 12     Working     Working 0.3990      98 0.9456      50 0.7862 0.7447
summary(probs_w)
#>           from          to    MIN MINtime    MAX MAXtime MEDIAN   MEAN
#> 1  Non-working        Dead 0.0036      50 0.6979      98 0.0743 0.1649
#> 4  Non-working Non-working 0.0000      98 0.8565      54 0.1780 0.3471
#> 7  Non-working     Retired 0.0296      50 0.7504      80 0.4309 0.4090
#> 10 Non-working     Working 0.0297      98 0.1312      50 0.0822 0.0791
#> 2      Retired        Dead 0.0156      57 0.3678      98 0.0253 0.0687
#> 5      Retired Non-working 0.0000      98 0.2013      50 0.0032 0.0404
#> 8      Retired     Retired 0.5857      50 0.9600      74 0.9001 0.8567
#> 11     Retired     Working 0.0052      98 0.1966      50 0.0115 0.0341
#> 3      Working        Dead 0.0009      50 0.3908      98 0.0208 0.0718
#> 6      Working Non-working 0.0000      98 0.0922      58 0.0138 0.0332
#> 9      Working     Retired 0.0078      50 0.2633      86 0.1700 0.1488
#> 12     Working     Working 0.4542      98 0.9311      50 0.7815 0.7462

# Plotting, men as example
probs_m |>  dtms_simplify() |> 
            ggplot(aes(x=time,y=P,color=to)) + 
            geom_ribbon(aes(ymin = CIlow, ymax = CIup,fill=to),alpha=0.5) +
            geom_line() + 
            facet_wrap(~from)
```

<img src="man/figures/README-example2-1.png" width="100%" />

``` r
 
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
#>                        Working Non-working  Retired    TOTAL
#> start:Working_50     13.307334    3.074605 13.53758 29.91952
#> start:Non-working_50  8.782989    6.496935 13.75116 29.03108
#> start:Retired_50      8.821763    3.905134 15.19067 27.91757
#> AVERAGE              12.445981    3.589228 13.65526 29.69047
    
dtms_expectancy(dtms=work,
                matrix=Tw,
                start_distr=Sw)
#>                        Working Non-working  Retired    TOTAL
#> start:Working_50     12.066485    4.386659 16.53253 32.98567
#> start:Non-working_50  7.445512    8.199122 16.77896 32.42359
#> start:Retired_50      7.836675    5.486704 18.24103 31.56440
#> AVERAGE              10.450573    5.607548 16.68838 32.74651

## Variant: ignoring retirement as a starting state (shown only for men)
limited <- c("Working","Non-working")

Smwr <- dtms_start(dtms=work,
                   data=estdata,
                   start_state=limited,
                   variables=list(Gender=0))

dtms_expectancy(dtms=work,
                matrix=Tm,
                start_state=limited,
                start_distr=Smwr)
#>                        Working Non-working  Retired    TOTAL
#> start:Working_50     13.307334    3.074605 13.53758 29.91952
#> start:Non-working_50  8.782989    6.496935 13.75116 29.03108
#> AVERAGE              12.650574    3.571394 13.56859 29.79056

## Lifetime risk of reaching retirement
dtms_risk(dtms=work,
          matrix=Tm,
          risk="Retired",
          start_distr=Sm)
#>     Working_50 Non-working_50     Retired_50        AVERAGE AVERAGE(COND.) 
#>      0.8828701      0.8806115      1.0000000      0.8888186      0.8825422
  
dtms_risk(dtms=work,
          matrix=Tw,
          risk="Retired",
          start_distr=Sw)
#>     Working_50 Non-working_50     Retired_50        AVERAGE AVERAGE(COND.) 
#>      0.9172407      0.9159547      1.0000000      0.9207352      0.9168268
  
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
#>                          MEAN  VARIANCE        SD MEDIAN     RISK0
#> start:Working_50     13.53758  97.35708  9.866969   13.0 0.1171299
#> start:Non-working_50 13.75116 103.70848 10.183736   13.0 0.1193885
#> start:Retired_50     15.19067 112.59117 10.610899   14.5 0.0000000
#> AVERAGE              13.65526  99.18227  9.959030   13.0 0.1111814
#> AVERAGE(COND.)       13.56859  98.28472  9.913865   13.0 0.1174578
summary(visitsw)
#>                          MEAN VARIANCE       SD MEDIAN      RISK0
#> start:Working_50     17.44977 112.1093 10.58817     18 0.08275934
#> start:Non-working_50 17.69491 117.8990 10.85813     18 0.08404534
#> start:Retired_50     18.74103 124.6118 11.16296     19 0.00000000
#> AVERAGE              17.58562 114.5507 10.70284     18 0.07926479
#> AVERAGE(COND.)       17.52865 113.9856 10.67640     18 0.08317318
  
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
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                          MEAN VARIANCE       SD MEDIAN      RISK0
#> start:Working_50     14.25887 42.52413 6.521053   14.5 0.00000000
#> start:Non-working_50 12.31587 50.90524 7.134791   12.5 0.00000000
#> start:Retired_50      0.00000  0.00000 0.000000    0.0 1.00000000
#> AVERAGE              13.13713 52.58739 7.251716   13.5 0.06011926
#> AVERAGE(COND.)       13.97744 44.20570 6.648737   13.5 0.00000000
summary(firstw)
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                          MEAN VARIANCE       SD MEDIAN      RISK0
#> start:Working_50     14.10718 40.00327 6.324814   14.5 0.00000000
#> start:Non-working_50 12.54709 46.49771 6.818923   12.5 0.00000000
#> start:Retired_50      0.00000  0.00000 0.000000    0.0 1.00000000
#> AVERAGE              12.91124 49.41214 7.029377   13.5 0.05103631
#> AVERAGE(COND.)       13.60562 42.62210 6.528560   13.5 0.00000000

## Last exit
  
# Leaving work to any state
last1m <- dtms_last(dtms=work,
                    matrix=Tm,
                    risk="Working",
                    start_distr=Sm)  
  
last1w <- dtms_last(dtms=work,
                    matrix=Tw,
                    risk="Working",
                    start_distr=Sw) 

summary(last1m)
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                          MEAN VARIANCE       SD MEDIAN RISK0
#> start:Working_50     16.50265 76.98676 8.774210   15.5    NA
#> start:Non-working_50 18.02302 68.14259 8.254853   17.5    NA
#> start:Retired_50     17.83146 69.33252 8.326615   17.5    NA
#> AVERAGE              16.73797 75.91238 8.712771   16.5    NA
#> AVERAGE(COND.)       17.97027 68.47754 8.275116   17.5    NA
summary(last1w)
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                          MEAN VARIANCE       SD MEDIAN RISK0
#> start:Working_50     16.15218 87.76963 9.368545   15.5    NA
#> start:Non-working_50 18.31738 77.06422 8.778623   17.5    NA
#> start:Retired_50     17.94411 78.93510 8.884543   17.5    NA
#> AVERAGE              16.78741 85.57486 9.250668   16.5    NA
#> AVERAGE(COND.)       18.26738 77.33095 8.793802   17.5    NA
  
# Leaving work for retirement
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
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                          MEAN VARIANCE       SD MEDIAN RISK0
#> start:Working_50     18.74988 64.64429 8.040167   18.5    NA
#> start:Non-working_50 19.72542 56.78054 7.535286   19.5    NA
#> start:Retired_50     19.60618 57.79592 7.602363   19.5    NA
#> AVERAGE              18.90783 63.49802 7.968565   18.5    NA
#> AVERAGE(COND.)       19.69276 57.06149 7.553906   19.5    NA
summary(last2w)
#> Warning in dtms_distr_summary(distr = object, ...): NAs introduced by coercion
#>                          MEAN VARIANCE       SD MEDIAN RISK0
#> start:Working_50     19.33660 73.65218 8.582085   19.5    NA
#> start:Non-working_50 20.62023 63.26114 7.953687   20.5    NA
#> start:Retired_50     20.40228 64.97334 8.060604   20.5    NA
#> AVERAGE              19.73733 70.75006 8.411305   19.5    NA
#> AVERAGE(COND.)       20.59143 63.49278 7.968235   20.5    NA
```

As already noted in the first example, consistent estimation of the
lifetime risk of reaching a state requires a different setup:

``` r
riskdata <- dtms_forward(data=workdata,
                         state="Retired",
                         dtms=work,
                         idvar="ID",
                         timevar="Age",
                         statevar="State")

riskdata <- dtms_format(data=riskdata,
                       dtms=work,
                       idvar="ID",
                       timevar="Age",
                       statevar="State")

riskdata <- dtms_clean(data=riskdata,
                       dtms=work)
#> Dropping  0  rows not in state space
#> Dropping  0  rows not in time range
#> Dropping  59719  rows starting or ending in NA
#> Dropping  68319  rows starting in absorbing state

riskdata$time2 <- riskdata$time^2

riskfit <- dtms_fit(data=riskdata,
                controls=c("Gender","time","time2"),
                package="mclogit")
#> 
#> Iteration 1 - deviance = 85313.14 - criterion = 0.9818182
#> Iteration 2 - deviance = 73020.05 - criterion = 0.168352
#> Iteration 3 - deviance = 69050.83 - criterion = 0.05748246
#> Iteration 4 - deviance = 67662.88 - criterion = 0.02051277
#> Iteration 5 - deviance = 67164.12 - criterion = 0.007425986
#> Iteration 6 - deviance = 66917.11 - criterion = 0.003691294
#> Iteration 7 - deviance = 66805.61 - criterion = 0.001668999
#> Iteration 8 - deviance = 66778.81 - criterion = 0.0004012381
#> Iteration 9 - deviance = 66769.97 - criterion = 0.0001324614
#> Iteration 10 - deviance = 66766.72 - criterion = 4.865712e-05
#> Iteration 11 - deviance = 66765.53 - criterion = 1.789762e-05
#> Iteration 12 - deviance = 66765.09 - criterion = 6.583857e-06
#> Iteration 13 - deviance = 66764.93 - criterion = 2.422023e-06
#> Iteration 14 - deviance = 66764.87 - criterion = 8.910069e-07
#> Iteration 15 - deviance = 66764.84 - criterion = 3.277824e-07
#> Iteration 16 - deviance = 66764.84 - criterion = 1.205843e-07
#> Iteration 17 - deviance = 66764.83 - criterion = 4.436046e-08
#> Iteration 18 - deviance = 66764.83 - criterion = 1.63193e-08
#> Iteration 19 - deviance = 66764.83 - criterion = 6.003535e-09
#> converged

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
#>     Working_50 Non-working_50     Retired_50 
#>      0.8861793      0.8951606      1.0000000
```

### Variance estimation

To use bootstrap methods, the function `dtms_boot()` is called, and
results can be conveniently viewed using `summary()`. The function
`dtms_boot()` needs data in transition format (argument `data`) and a
`dtms` object. The argument `method` is used to choose the bootstrap
method; here, we use the block bootstrap. In case the block bootstrap is
used the argument `idvar` needs to be specified, which takes the name of
the variable with the unit identifier. The argument `rep` sets the
number of bootstrap replications, and the argument `parallel` can be set
to `TRUE` to enable parallel processing using the packages
[foreach](https://cran.r-project.org/web/packages/foreach) and
[doParallel](https://cran.r-project.org/web/packages/doParallel).

Further required is the argument `fun`. This is a function which should
have two arguments, one called `data` and one called `dtms`. These are
used to pass the corresponding arguments from `dtms_boot()`. Other than
this the function can contain anything the user is interested in. In the
example above, the function is called `bootfun`. It estimates transition
probabilities, puts them into a transition matrix, and then calculates
state expectancies. Each bootstrap replication of the data is passed to
the function specified by `fun`, and the results are saved in a list
with as many entries as there are replications. The format of each entry
of the list obviously depends on the definition of `fun`.

The result of calling `summary(bootresults)` is by default a list with
two entries which together provide the bootstrap percentiles, i.e., the
bootstrap confidence interval. The entries have the structure defined by
`fun`. In this example, the upper half of each entry are the state
expectancies for men, while the lower half are the state expectancies
for women. For instance, the 95% confidence interval for the average
lifetime spent working ranges from 12.17 years to 12.67 years for men.
The returned percentiles can be controlled with the arguments of the
summary method, `dtms_boot_summary()`.

``` r
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
                         rep=50,
                         method="block",
                         parallel=TRUE)

summary(bootresults)
#> $`2.5%`
#>                        Working Non-working  Retired    TOTAL
#> start:Working_50     13.027111    2.917580 13.18544 29.46857
#> start:Non-working_50  8.483431    6.269246 13.26474 28.34093
#> start:Retired_50      8.406719    3.715340 14.63665 27.07383
#> AVERAGE              12.169058    3.450685 13.26579 29.17974
#> start:Working_50     11.765995    4.233728 16.10671 32.59493
#> start:Non-working_50  7.105853    7.979484 16.30210 32.07353
#> start:Retired_50      7.411987    5.210598 17.68557 30.98511
#> AVERAGE              10.082227    5.429166 16.23917 32.35578
#> 
#> $`97.5%`
#>                        Working Non-working  Retired    TOTAL
#> start:Working_50     13.653372    3.223444 14.08771 30.40601
#> start:Non-working_50  9.150639    6.750686 14.30660 29.58953
#> start:Retired_50      9.239235    4.131469 15.62659 28.49331
#> AVERAGE              12.820168    3.759309 14.21257 30.17806
#> start:Working_50     12.281845    4.616786 16.93534 33.42463
#> start:Non-working_50  7.687136    8.457582 17.20335 32.84186
#> start:Retired_50      8.163383    5.806967 18.72468 32.06217
#> AVERAGE              10.708743    5.793324 17.09537 33.17589
```

## Using dtms with irregular intervals

Longitudinal data is often not observed in regular intervals, but only
irregular. For instance, respondents in a biannual longitudinal survey
might not be interviewed exactly every two years, but sometimes only one
and a half years might have passed, or two and a half years, or
something in between. This applies, for instance, to the U.S. Health and
Retirement Study (HRS).

A complete application showing how to handle irregular intervals with
`dtms` using HRS data can be found online:
<https://github.com/christiandudel/hrs_hwle> The basic steps are as
follows, using the HRS code as an example. First, the `dtms`object used
for transforming the data using `dtms_format()` needs to include not
one, but several values for the steplength:

``` r
hrsdtms <- dtms(transient=c("working","not-working"),
                absorbing="dead",
                timescale=seq(50,98,1),
                timestep=1:3)
```

Here, values 1, 2, and 3 are used. Using this object with
`dtms_format()` will keep all transitions which have an interval of 1,
2, or 3 time units. If in another application observations are every
half year with a maximum interval of 3 and a half years, the argument
could be set to `seq(0.5,3.5,by=0.5)`. Moreover the values of the
timescale are set to increase by units of 1 from 50 to 98, as
respondents in the HRS cover the population aged 50 and older, and age
is measured in completed years; 98 is used as the highest age in which
transitions can start. Already note that this `dtms` object will not be
used to predict transition probabilities.

Second, using this `dtms`object with `dtms_format()` should use the
argument `steplength=TRUE`, like:

``` r
estdata <- dtms_format(data=hrsdata,
                       dtms=hrsdtms,
                       steplength=TRUE)
```

In the resulting data set, there will be a variable called `steplength`
which for each transition indicates the interval. The name of this
variable can be controlled using the argument `stepvar`.

Third, the interval should be included as a predictor when estimating
the transition probabilities:

``` r
fit <- dtms_fit(data=somedata,
                controls=c("time","steplength"))
```

The above example will include the interval width as a linear predictor.
Whether this is meaningful, or some other specification should be used,
will depend on the specific application.

Fourth, a new `dtms` object is created which is used for predicting
transition probabilities:

``` r
hrspredict <- dtms(transient=c("working","not-working"),
                     absorbing="dead",
                     timescale=seq(50,98,2))
```

In contrast to the previous `dtms` object, this uses only one fixed
interval length. Which fixed value is used is a choice of the user, and
will likely depend on the distribution of the intervals; if, for
instance, most observations are two years apart, and only a few one or
three years, then using a fixed interval of two years will be most
supported by the data. This `dtms` object is next combined with
`dtms_transitions`:

``` r
dtms_transitions(dtms=hrspredict,
                 model=fit,
                 controls=list(time=seq(50,98,2),
                               steplength=2))
```

This predicts transition probabilities transitioning from age 50 to age
52, from age 52 to age 54, and so on up to transitioning from 98 to 100
and then dying; thus it also uses a fixed interval of 2, which is
explicit not only in the values of the time scale used for prediction,
but also the value for the variable `steplength`. The resulting
transition probabilities can be used like any other transition
probabilities.

## Combining dtms with other software

`dtms` can easily be used without doing the full workflow in `dtms`
itself. For instance, transition probabilities could be calculated with
a different software and `dtms` could then only used to calculate state
expectancies or similar indicators.

To do this, two things need to be kept in mind. The first is that
objects from the `dtms` package have a specific structure, and data or
results from other software need to have the correct structure. For
instance, the function `dtms_transitions()` outputs transition
probabilities in a certain way. Users who want to import transition
probabilities should arrange these probabilities the same way.

Second, most objects generated by `dtms` have specific object classes.
For instance, objects generated by `dtms_transitions()` and
`dtms_nonparametric()` have two classes: `dtms_probs` and `data.frame`.
User-generated objects need to have the right class(es) to work
correctly with `dtms`.

## Using dtms with secure data environments

In many secure data environments, installing `dtms` will not be
possible. There are two major ways to still use `dtms`. The first is to
export a tabulation of transitions and covariates from the secure
environment, similar to what `dtms_aggregate()` does provide; i.e., for
each transition, such as from some state A to some other state B, there
is a count of transitions by the values of the time scale, and
potentially further covariates. These counts can then be used with the
`weights` argument of functions such as `dtms_fit()`.

Sometimes it is possible to import existing code to a secure data
environment. In such a case, one of the two files in the folder
`combined` on GitHub can be used. The file `all.R` is (almost) all code
from `dtms`. Sourcing it should provide a lot of the functionality of
the package; dependencies are still required, though. The file
`selected.R` only includes a subset of the functions and removes most of
the documentation, but is much smaller than `all.R`. Again dependencies
might be required.

## References

Papers using `dtms` for substantive questions:

- Moretti, M., Korhonen, K., van Raalte, A., Riffe, T., Martikainen, P.
  (2025): Evolution of widowhood lifespan and its gender and educational
  inequalities in Finland over three decades. Demography 65(2):
  1635-1660. <https://doi.org/10.1215/00703370-12269717>

- Hiilamo, A., Pitkänen, J., Moretti, M., Martikainen, P., Myrskylä, M.
  (2025): Children’s out-of-home care in Finland, 1993–2020: lifetime
  risks, expectancies, exit routes, and number of placements for
  synthetic cohorts. Child Abuse & Neglect 169(1): 107626.
  <https://doi.org/10.1016/j.chiabu.2025.107626>

- Hiilamo, A., Hermansen, Å. (2025): Financial strain in Norway: the
  lifetime risk of and expected time spent in payment problems. MPIDR
  Working Paper WP-2025-006.
  <https://dx.doi.org/10.4054/MPIDR-WP-2025-006>

- Feraldi, A., Dudel, C. (2025): Smoking and the length of working life:
  an examination using the U.S. Health and Retirement Study. MPIDR
  Working Paper WP-2025-017.
  <https://dx.doi.org/10.4054/MPIDR-WP-2025-017>

- Lam, A. A., Keenan, K., Kulu, H., Myrskylä, M. (2024): Working longer
  despite poorer health? Inequalities in working and health expectancies
  at older ages in South Korea. MPIDR Working Paper WP-2024-022.
  <https://dx.doi.org/10.4054/MPIDR-WP-2024-022>

Methodological papers:

- Schneider, D. C. (2023): Statistical inference for discrete-time
  multistate models: asymptotic covariance matrices, partial age ranges,
  and group contrasts. MPIDR Working Paper WP-2023-041.
  <https://dx.doi.org/10.4054/MPIDR-WP-2023-041>

- Dudel, C. (2021): Expanding the Markov chain tool box: distributions
  of occupation times and waiting times. Sociological Methods & Research
  50: 401-428. <https://doi.org/10.1177/0049124118782541>

- Dudel, C., Myrskylä, M. (2020): Estimating the number and length of
  episodes in disability using a Markov chain approach. Population
  Health Metrics 18: 15. <https://doi.org/10.1186/s12963-020-00217-0>
