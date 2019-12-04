context("basic do_many test")

# objects to be used on tests
my_dir <- "../tmp_test/"
sp <- names(example_occs)[1]
sp_coord <- example_occs[[1]]
part <- 4
mod_dir <- paste0(my_dir, sp, "/present/partitions/")
algos <- c("bioclim", "glm", "rf", "svme", "brt")
n.algos <- length(algos)

test_that("do_many produces model and stats file", {
  my_mods <- do_many(species_name = sp,
                     predictors = example_vars,
                     models_dir = my_dir,
                     bioclim = TRUE,
                     glm = TRUE,
                     rf = TRUE,
                     svme = TRUE,
                     brt = TRUE)
  # does it have two txt files (matrix and evaluate) ?
  expect_length(list.files(path = mod_dir,
                           pattern = "confusion_matrices_.*csv"),
                part * n.algos)
  expect_length(list.files(path = mod_dir,
                           pattern = "evaluate.*csv"),
                part * n.algos)
  expect_length(list.files(path = mod_dir,
                           pattern = "eval_mod.*csv"),
                part * n.algos)
  # does it have tif file (model)?
  expect_length(list.files(path = mod_dir,
                           pattern = ".*.tif"),
                part * n.algos)
})
