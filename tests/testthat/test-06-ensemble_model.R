context("basic ensemble test")

my_dir <- "../tmp_test/"
sp <- names(coordenadas)[1]
sp_coord <- coordenadas[[1]]
ens_dir <- paste0(my_dir, sp, "/present/ensemble/")

test_that("output is generated", {
  my_ens <- ensemble_model(species_name=sp,
                           occurrences = sp_coord,
                           models_dir = my_dir,
                           which_final = c("bin_consensus"),
                           overwrite = TRUE)
  # does it have png and tif files
  expect_length(list.files(path = ens_dir,
                           pattern=".*png"), 4)
  expect_length(list.files(path = ens_dir,
                           pattern=".*.tif"), 4)
})

# delete test dir at the end
teardown(unlink("../tmp_test", recursive = TRUE))