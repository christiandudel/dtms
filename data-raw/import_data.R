# Sources/code on Github: https://www.github.com/christiandudel/dtms_data
### Package
library(devtools)

### Import simple example
load(url("https://github.com/christiandudel/dtms_data/raw/main/Output/simple.rda"))
simpledata <- as.data.frame(simdata)
usethis::use_data(simpledata,overwrite=TRUE)

### Import simulated working data
# based on: https://doi.org/10.1007/s13524-017-0619-6
load(url("https://github.com/christiandudel/dtms_data/raw/main/Output/work.rda"))
workdata <- as.data.frame(simdata)
usethis::use_data(workdata,overwrite=TRUE)
