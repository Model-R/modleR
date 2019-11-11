context("basic final_model test")

my_dir <- "../tmp_test/"
sp <- names(example_occs)[1]
n.algos <- 2
final_dir <- paste0(my_dir, sp, "/present/final_models/")

test_that("final_model generates joint model per algorithm", {
        sp_final <- final_model(species_name = sp,
                                models_dir = my_dir,
                                #algorithms = "bioclim",
                                select_partitions = TRUE,
                                select_par = "TSSmax",
                                select_par_val = 0,
                                which_models = c("bin_consensus"),
                                consensus_level = 0.5,
                                overwrite = TRUE)
 # does it have one csv file (final statistics) ?
 expect_length(list.files(path = final_dir,
                          pattern = ".*final_statistics.csv"), 1)
 # does it have the model png and tif file
 expect_length(list.files(path = final_dir,
                          pattern = ".*png"), n.algos)
 expect_length(list.files(path = final_dir,
                          pattern = ".*.tif"), n.algos)
})
