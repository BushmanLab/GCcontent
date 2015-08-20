[![Travis-CI Build Status](https://travis-ci.org/BushmanLab/GCcontent.svg?branch=master)](https://travis-ci.org/BushmanLab/GCcontent)

[![codecov.io](http://codecov.io/github/BushmanLab/GCcontent/coverage.svg?branch=master)](http://codecov.io/github/BushmanLab/GCcontent?branch=master)


# GCcontent

calculate GC percentage for reference genomes with multiple windows

# Dependency

IntSiteRetriver(only for testing and timing) should be installed: 
    https://github.com/anatolydryga/intSiteRetriever

# Testing 

Run in the R console:

```bash
library(testthat)
devtools::test()
```

# Timing 
```
source('timing_GCcontent.R')
time_GC()
```
