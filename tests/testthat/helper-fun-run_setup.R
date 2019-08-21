run_setup <- function(my_seed, ...){
  setup_sdmdata(species_name=sp, 
                occurrences=sp_coord, 
                example_vars, 
                seed=my_seed, 
                boot_n = 4,
                models_dir = my_dir, 
                partition_type = "crossvalidation",
                cv_partitions = 5,
                cv_n = 1,
                ...)
}