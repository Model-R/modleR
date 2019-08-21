# test file for setup_sdmdata function

## goal: test output of function

context("basic setup test")

test_that("testing that seed works", {
  sp <- names(coordenadas)[1]
  sp_coord <- coordenadas[[1]]
  my_seed <- 42
  my_dir <- "../01_test/"
  my_setup <- run_setup(my_seed)
  expect_message(run_setup(my_seed), 
                 "same metadata, no need to run data partition")
  expect_message(run_setup(123), 
                 "saving metadata")
})

test_that("all outputs were generated", {
  setup_dir <- paste0(my_dir, sp, "/present/data_setup/")
  # does it have two metadata and sdmdata txt files?
  expect_length(list.files(path = setup_dir, 
                                 pattern="metadata.txt"), 1)
  expect_length(list.files(path = setup_dir, 
                                 pattern="sdmdata.txt"), 1)
  # does it have sdmdata png file?
  expect_equal(length(list.files(path = setup_dir, 
                                 pattern="sdmdata_.*png")), 1)
  
})

test_that("setup has numeric values", {
  lapply(apply(my_setup, 2, is.numeric), expect_true)
})
