# test file for do_any function

## goal: test output of function

test_that("do_many produces model and stats file", {
  my_dir <- "../01_test/"
  my_mod <- do_any(species_name=sp,
                   predictors=example_vars,
                   models_dir=my_dir,
                   algo = "bioclim")
  
})