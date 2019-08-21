context("basic do_any test")

test_that("do_any produces model and stats file", {
  my_dir <- "../01_test/"
  my_mod <- do_any(species_name=sp,
                   predictors=example_vars,
                   models_dir=my_dir,
                   algo = "bioclim")
  
})