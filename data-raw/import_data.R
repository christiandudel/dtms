### Package
library(devtools)

### Import simple example, currently only exists locally
load("U:/Documents/__CODE__/Markov Chain Simulation/Data/simple.rda")
simpledata <- as.data.frame(simdata)
usethis::use_data(simpledata,overwrite=TRUE)

### Import simulated HRS data, currently only exists locally
# based on: https://doi.org/10.1007/s13524-017-0619-6
load("U:/Documents/__CODE__/Markov Chain Simulation/Data/hrs.rda")
hrsdata <- as.data.frame(simdata)
usethis::use_data(hrsdata,overwrite=TRUE)
