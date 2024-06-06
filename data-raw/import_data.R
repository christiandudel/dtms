# Sources/code on Github: https://www.github.com/christiandudel/dtms_data
### Package
library(devtools)

### Import simple example, currently only exists locally
load(url("https://github.com/christiandudel/dtms_data/raw/main/Output/simple.rda"))
simpledata <- as.data.frame(simdata)
usethis::use_data(simpledata,overwrite=TRUE)

### Import simulated HRS data, currently only exists locally
# based on: https://doi.org/10.1007/s13524-017-0619-6
load(url("https://github.com/christiandudel/dtms_data/raw/main/Output/hrs.rda"))
hrsdata <- as.data.frame(simdata)
usethis::use_data(hrsdata,overwrite=TRUE)
