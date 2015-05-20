# GCcontent
calculate GC percentage for reference genomes with multiple windows

# Dependency

IntSiteRetriver(only for testing and timing): 
    https://github.com/anatolydryga/intSiteRetriever

has to be cloned:

```bash
git clone https://github.com/anatolydryga/intSiteRetriever

```


# Testing 

Run in the R console:

```bash
library(testthat)
source('GCcontent.R')
test_dir(".")
```

# Timing 
```
source('timing_GCcontent.R')
time_GC()
```
