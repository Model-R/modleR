context("basic ensemble test")

my_dir <- "../tmp_test"
sp <- names(example_occs)[1]
sp_coord <- example_occs[[1]]
ens_dir <- paste0(my_dir, "/", sp, "/present/ensemble/")
#lists all the options from which
ensemble_options <- c("best", "average", "weighted_average", "median",
                      "frequency", "consensus", "pca")
#creates a random combination of the options- ANY ONE SHOULD WORK
which <- sample(ensemble_options, sample(1:length(ensemble_options), 1))

test_that("ensemble output", {
  my_ens <- ensemble_model(species_name = sp,
                           occurrences = sp_coord,
                           performance_metric = "pROC",
                           dismo_threshold = "spec_sens",
                           which_ensemble = which,
                           which_final = "raw_mean",
                           models_dir = my_dir,
                           overwrite = TRUE,
                           uncertainty = TRUE,
                           consensus_level = 0.5,
                           )
  #the object is a raster
  expect_s4_class(my_ens, "RasterStack")
  # does it have png and tif files
  expect_length(list.files(path = ens_dir,
                            pattern = paste0("ensemble_", which)), 1)
  expect_length(list.files(path = ens_dir,
                           pattern = "uncertainty.tif"), 1)
})

# delete test dir at the end
teardown(unlink("../tmp_test", recursive = TRUE))
