# test file for setup_sdmdata function

## goal: test output of function

context("basic setup test")

# creating objects to be used on tests
my_seed <- 42
my_dir <- "../01_test/"
sp <- names(coordenadas)[1]
sp_coord <- coordenadas[[1]]
# helper function
run_setup <- function(my_seed, ...){
  setup_sdmdata(species_name=sp,
                occurrences=sp_coord,
                example_vars,
                seed=my_seed,
                boot_n = 4,
                models_dir = my_dir,
                ...)
}
my_setup <- run_setup(my_seed)

# the test itself
test_that("testing that seed works", {
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
  expect_length(list.files(path = setup_dir, 
                                 pattern="sdmdata_.*png"), 1)
  
})

test_that("setup has numeric values", {
  lapply(apply(my_setup, 2, is.numeric), expect_true)
})
