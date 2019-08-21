# loading packages
library(vdiffr) # for plot tests
library(testthat) # for basic tests
library(modleR)

# remove test dir to start from scratch
unlink("tests/01_test", recursive = TRUE)

#runs all tests on testthat dir
test_check("modleR")

