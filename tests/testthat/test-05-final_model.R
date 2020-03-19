context("basic final_model test")

my_dir <- "../tmp_test/"
sp <- names(example_occs)[1]
algos <- c("bioclim", "maxnet", "svme")
n.algos <- length(algos)
final_folder <- paste0(my_dir, sp, "/present/final_models/")
#lists all the options from final
#final_options <- c("raw_mean", "raw_mean_th", "raw_mean_cut", "bin_mean",
                  # "bin_consensus")
#creates a random combination of the options- ANY ONE SHOULD WORK
#which <- sample(final_options, sample(1:length(final_options), 1))
which <- "raw_mean"
test_that("final_model generates joint model per algorithm", {
        sp_final <- final_model(species_name = sp,
                                models_dir = my_dir,
                                which_models = which,
                                consensus_level = 0.5,
                                overwrite = TRUE)
 # does it have one csv file (final statistics) ?
 expect_length(list.files(path = final_folder,
                          pattern = ".*final_statistics.csv"), 1)
 # does it have one csv file (mean statistics) ?
 expect_length(list.files(path = final_folder,
                          pattern = ".*mean_statistics.csv"), 1)
 # does it have the model png and tif file? there must be one per algorithm PER
 #WHICH option
 expect_length(list.files(path = final_folder,
                          pattern = ".*png"), n.algos * length(which))
 expect_length(list.files(path = final_folder,
                          pattern = ".*.tif"), n.algos * length(which))
})
