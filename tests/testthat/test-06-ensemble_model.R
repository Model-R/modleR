context("basic ensemble test")

my_dir <- "../tmp_test/"
sp <- names(example_occs)[1]
sp_coord <- example_occs[[1]]
ens_dir <- paste0(my_dir, sp, "/present/ensemble/")
final_dir <- paste0(my_dir, sp, "/present/final_models/")

test_that("output is a rasterStack", {
  my_ens <- ensemble_model(species_name = sp,
                           occurrences = sp_coord,
                           models_dir = my_dir,
                           final_dir = final_dir,
                           overwrite = TRUE,
                           uncertainty = TRUE)
  # does it have png and tif files
  # expect_length(list.files(path = ens_dir,
  #                          pattern = ".*png"), 4)
  expect_s4_class(my_ens, "raster")
  #expect_length(list.files(path = ens_dir,
   #                         pattern = "ensemble_average"), 1)
  #expect_length(list.files(path = ens_dir,
   #                         pattern = "uncertainty.tif"), 1)
})

# delete test dir at the end
teardown(unlink("../tmp_test", recursive = TRUE))
