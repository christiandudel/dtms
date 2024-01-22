
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

## Example

This is a placeholder for an example. It currently only shows how to
define an ‘dtms’ object

``` r
library(dtms)
dtms(transient=c("A","B"),
     absorbing="X",
     time=1:10)
#> $transient
#> [1] "A" "B"
#> 
#> $absorbing
#> [1] "X"
#> 
#> $time
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $timestep
#> [1] 1
#> 
#> attr(,"class")
#> [1] "dtms"
## basic example code
```
