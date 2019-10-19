#' Prepares the dataset to perform ENM
#'
#' This function takes the occurrence points files and the predictor layers and
#' executes data cleaning, data partitioning, pseudo-absence point sampling and
#' variable selection according to their correlation. It saves the metadata and
#' sdmdata files into the hard disk.
#'
#' @inheritParams create_buffer
#' @param species_name A character string with the species name. Because species
#' name will be used as a directory name, avoid non-ASCII characters, spaces and
#' punctuation marks.
#' Recommendation is to adopt "Genus_species" format. See names in
#' \code{\link{example_occs}} as an example
#' @param occurrences A data frame with occurrence data. Data must have at least
#'  columns with latitude and longitude values of species occurrences.
#' See \code{\link{example_occs}} as an example
#' @param lon The name of the longitude column. Defaults to "lon"
#' @param lat The name of the latitude column. Defaults to "lat"
#' @param predictors A Raster or RasterStack object with the environmental
#' raster layers
#' @param seed Random number generator for reproducibility purposes. Used for
#'  sampling pseudoabsences
#' @param real_absences User-defined absence points
#' @param geo_filt Logical, delete occurrences that are too close to each other?
#'  See \insertCite{varela_environmental_2014;textual}{modleR}
#' @param geo_filt_dist The distance of the geographic filter in the unit of the
#' predictor raster, see \insertCite{varela_environmental_2014;textual}{modleR}
#' @param select_variables Logical. Whether a variable selection should be performed. It excludes highly correlated environmental
#'  variables. If TRUE, \code{cutoff} and \code{percent} parameters must be specified
#' @param models_dir Folder path to save the output files. Defaults to
#' "\code{./models}"
#' @param plot_sdmdata Logical, whether png files will be written
#' @param n_back Number of pseudoabsence points. Default is 1,000
#' @param partition_type Character. Type of data partitioning scheme, either
#' "\code{bootstrap}" or k-fold "\code{crossvalidation}". If set to bootstrap, \code{boot_proportion} and \code{boot_n} must be specified. If set to crossvalidation, \code{cv_n} and \code{cv_partitions} must be specified
#' @param boot_proportion Numerical 0 to 1, proportion of points to be sampled
#' for bootstrap
#' @param boot_n Number of bootstrap runs
#' @param cv_partitions Number of partitions in the crossvalidation
#' @param cv_n Number of crossvalidation runs
#' @param clean_dupl Logical. If TRUE, removes points with the same longitude and
#' latitude
#' @param clean_nas Logical. If TRUE, removes points that are outside the bounds	#' of the raster
#' @param clean_uni Logical. If TRUE, selects only one point per pixel
#' @param cutoff Cutoff value of correlation between variables to exclude
#' environmental layers
#' Default is to exclude environmental variables with correlation > 0.8
#' @param percent percentage of the raster values to be sampled to calculate the
#'  correlation
#' @param ... Othre parameters from \code{\link{create_buffer}}
#' @return Returns a data frame with the groups for each run (in columns called
#' cv.1, cv.2 or boot.1, boot.2), a presence/absence vector, the geographical
#' coordinates of the occurrence and pseudoabsence points, and the associated
#' environmental variables (either all the layers or the selected ones if
#' \code{select_variables = TRUE}). Function writes on disk (inside subfolder
#' at \code{models_dir} directory) a text file named sdmdata.csv that will be used
#' by \code{\link{do_any}} or \code{\link{do_many}}
#' @examples
#' sp <- names(example_occs)[1]
#' sp_coord <- example_occs[[1]]
#' sp_setup <- setup_sdmdata(species_name = sp,
#'                           occurrences = sp_coord,
#'                           example_vars)
#' head(sp_setup)
#' @references
#'     \insertAllCited{}
#' @seealso \code{\link{create_buffer}}
#' @seealso \code{\link[dismo]{gridSample}}
#' @importFrom utils write.table
#' @importFrom stats cor
#' @export
#'
#'
# tabela de valores
setup_sdmdata <- function(species_name,
                          occurrences,
                          predictors,
                          models_dir = "./models",
                          real_absences = NULL,
                          lon = "lon",
                          lat = "lat",
                          buffer_type = NULL,
                          dist_buf = NULL,
                          env_buffer = FALSE,
                          env_distance = "centroid",
                          dist_min = NULL,
                          buffer_shape = NULL,
                          max_env_dist = 0.5,
                          write_buffer = FALSE,
                          seed = NULL,
                          clean_dupl = FALSE,
                          clean_nas = FALSE,
                          clean_uni = FALSE,
                          geo_filt = FALSE,
                          geo_filt_dist = NULL,
                          select_variables = FALSE,
                          cutoff = 0.8,
                          percent = 0.8,
                          plot_sdmdata = TRUE,
                          n_back = 1000,
                          partition_type = c("bootstrap"),
                          boot_n = 1,
                          boot_proportion = 0.7,
                          cv_n = NULL,
                          cv_partitions = NULL) {
    # replacing characters not welcome in species name
    # characters to avoid in file and dir names
    avoid_chars <- intToUtf8(c(91, 62, 33, 180, 60, 35, 63, 38, 47, 92, 46, 93))
    print_avoid <- intToUtf8(c(62, 33, 180, 60, 35, 63, 38, 47, 92, 46))
    if (grepl(avoid_chars, species_name) == TRUE) {
        species_name <- gsub(avoid_chars, "", species_name)
        warning(cat(paste0('You entered a bad character (any in "',
                           print_avoid,
                           '") in the species name and we removed it for you')))
        }
    if (file.exists(models_dir) == FALSE)
        dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
    if (file.exists(paste(models_dir, species_name, sep = "/")) == FALSE)
        dir.create(paste(models_dir, species_name, sep = "/"))

    #creates a separate folder for sdmdata and metadata
    setup.folder <-
        paste(models_dir, species_name, "present", "data_setup", sep = "/")
    if (file.exists(setup.folder) == FALSE)
        dir.create(setup.folder, recursive = TRUE)

    ## checking latitude and longitude columns
    if (all(c(lon, lat) %in% names(occurrences))) {
        occurrences <- occurrences[, c(lon, lat)]
        names(occurrences) <- c("lon", "lat")
    } else {
        stop("Coordinate column names do not match. Either rename to `lon` and `lat` or specify")
    }
    #creates metadata for this run
    original_n <- nrow(occurrences)
    original_n_back <- n_back
    original_predictors <- paste(names(predictors), collapse = "-")

    metadata_new <- data.frame(
        species_name = as.character(species_name),
        original_predictors = original_predictors,
        select_variables = select_variables,
        original.n = as.integer(original_n),
        original.n.back = as.integer(original_n_back),
        buffer_type = ifelse(is.null(buffer_type), NA, buffer_type),
        dist_buf = ifelse(is.null(dist_buf), NA, dist_buf),
        seed = ifelse(is.null(seed), NA, as.integer(seed)),
        res.x = res(predictors)[1],
        res.y = res(predictors)[2],
        clean_dupl = ifelse(is.null(clean_dupl), NA, clean_dupl),
        clean_nas = ifelse(is.null(clean_nas), NA, clean_nas),
        clean_uni = ifelse(is.null(clean_uni), NA, clean_uni),
        geo_filt = geo_filt,
        geo_filt_dist = ifelse(is.null(geo_filt_dist), NA, as.integer(geo_filt_dist)),
        models_dir = models_dir,
        partition = partition_type,
        boot_proportion = ifelse(is.null(boot_proportion), NA, boot_proportion),
        boot_n = ifelse(is.null(boot_n), NA, as.integer(boot_n)),
        cv_partitions = ifelse(is.null(cv_partitions), NA, as.integer(cv_partitions)),
        cv_n = ifelse(is.null(cv_n), NA, as.integer(cv_n))#,
        #row.names = 1
        )

        #checking metadata----
    if (file.exists(paste(setup.folder, "metadata.csv", sep = "/"))) {
        message("metadata file found, checking metadata \n")
        metadata_old <- read.csv(paste(setup.folder, "metadata.csv", sep = "/"),
                                   as.is = FALSE) #row.names = 1)
        # removes columns that dont exist yet for comparison
        metadata_old <- metadata_old[,
                                     setdiff(names(metadata_old),
                                             c("final.n", "final.n.back", "selected_predictors"))]
        if (all(all.equal(metadata_old, metadata_new) == TRUE)) {
            message("same metadata, no need to run data partition")
            sdmdata <- read.csv(paste(setup.folder, "sdmdata.csv", sep = "/"), as.is = FALSE)
            return(sdmdata)
            }
    }

    ##cleaning occurrences with clean and geo_filt----
    message("running data setup")
    message("cleaning data")
    occurrences <-
        clean(occurrences = occurrences,
              predictors = predictors,
              clean_dupl = clean_dupl,
              clean_nas = clean_nas,
              clean_uni = clean_uni)

    if (geo_filt == TRUE) {
        message("applying a geographical filter")
        occurrences <-
            geo_filt(occurrences = occurrences,
                     min_distance = geo_filt_dist)
    }
    final_n <- nrow(occurrences)


    #background selection:

    #first option: there is a buffer
    if (!is.null(buffer_type)) {
        if (buffer_type %in% c("mean", "maximum", "median", "distance", "user", "environmental_distance")) {
            message("creating buffer")
            pbuffr <- create_buffer(occurrences = occurrences,
                                    models_dir = models_dir,
                                    species_name = species_name,
                                    buffer_type = buffer_type,
                                    predictors = predictors,
                                    dist_buf = dist_buf, #tiene que estar
                                    dist_min = dist_min,
                                    buffer_shape = buffer_shape,
                                    env_distance = env_distance,
                                    max_env_dist = max_env_dist,
                                    write_buffer = write_buffer)

        } else {
        warning("buffer_type not recognized, returning predictors")
            pbuffr <- predictors
            }
        }
    # second option: there is no buffer
    else pbuffr <- predictors # second option, there is no buffer

    # absences
    #first option: user-supplied absences
    if (!is.null(real_absences)) {
        backgr <- real_absences[, c(lon, lat)]
        n_back_mod <- nrow(backgr)
    } else {
        #sampling pseudoabsence points
                #checks if there will be enough cells to sample pseudoabsences from
                vals <- values(pbuffr[[1]])
                available_cells <- sum(!is.na(vals)) - nrow(occurrences)
                # and corrects accordingly
                if (available_cells < n_back) {
                    n_back_mod <- available_cells
                    message(paste(available_cells, "available cells"))
                    message(paste("Using", n_back_mod, "pseudoabsences", "\n"))
                } else {
                    n_back_mod <- n_back
                }
        #Now it does the sampling
                message(paste("sampling pseudoabsence points with", buffer_type, "buffer"))
        set.seed(seed)
                backgr <- dismo::randomPoints(mask = pbuffr,
                                              n = n_back_mod,
                                              p = occurrences,
                                              excludep = TRUE)
    }
    colnames(backgr) <- c("lon", "lat")

    # Seleccionando variables if sel_vars ==TRUE
    if (select_variables == TRUE) {
    message(paste("selecting variables...", "\n"))
        predictors <- select_variables(predictors = predictors,
                                       buffer = pbuffr,
                                       cutoff = cutoff,
                                       percent = percent)
        }

    # edit metadata
    metadata_new$selected_predictors <- paste(names(predictors), collapse = "-")
    metadata_new$final.n <- as.integer(final_n)
    metadata_new$final.n.back <- as.integer(n_back_mod)

    message(paste("saving metadata"), "\n")
    write.table(metadata_new, file = paste(setup.folder, "metadata.csv", sep = "/"),
                sep = ",", col.names = TRUE, row.names = FALSE)

    # cria a tabela de valores
    message("extracting environmental data")
    presvals <- raster::extract(predictors, occurrences)
    # Extraindo dados ambientais dos bckgr
    message("extracting background data")
    backvals <- raster::extract(predictors, backgr)
    if (any(complete.cases(backvals) == FALSE)) {
        backgr   <- backgr[complete.cases(backvals), ]
        backvals <- raster::extract(predictors, backgr)
        warning(paste("Your background data had NA values, ", nrow(backvals),
                      "points were retained"))
        }

    pa <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))
    pres <- cbind(occurrences, presvals)
    back <- cbind(backgr, backvals)
    sdmdata <- cbind(pa, rbind(pres, back))
    # Data partition-----
    message("performing data partition")
    #Crossvalidation, repetated crossvalidation and jacknife
    if (partition_type == "crossvalidation") {
        if (nrow(occurrences) < 11) {
            message("data set has 10 occurrences or less, forcing jacknife")
            #forces jacknife
            cv_partitions <- nrow(occurrences)
            cv_n <- 1
        }
        if (is.null(cv_n)) stop("cv_n must be specified in crossvalidation")
        if (is.null(cv_partitions)) stop("cv_partitions must be specified in crossvalidation")
        if (cv_n == 1) {
            #Crossvalidation
            set.seed(seed)  #reproducibility
            group <- dismo::kfold(occurrences, cv_partitions)
            set.seed(seed)
            bg.grp <- dismo::kfold(backgr, cv_partitions)
            group.all <- c(group, bg.grp)
        }
        if (cv_n > 1) {
            # Repeated CV
            cv.pres <- replicate(n = cv_n,
                                 dismo::kfold(occurrences, cv_partitions))
            dimnames(cv.pres) <- list(NULL, paste0("cv", 1:cv_n))
            cv.back <- replicate(n = cv_n,
                                 dismo::kfold(backgr, cv_partitions))
            dimnames(cv.back) <- list(NULL, paste0("cv", 1:cv_n))
            cv.matrix <- rbind(cv.pres, cv.back)
        }
    }
    # Bootstrap
    if (partition_type == "bootstrap") {
        if (boot_proportion > 1 | boot_proportion <= 0)
            stop("bootstrap training set proportion must be between 0 and 1")
        if (is.null(boot_n))
            stop("boot_n must be specified")
    boot.pres <- replicate(n = boot_n,
                           sample(
                               x = seq_along(1:nrow(occurrences)),
                               size = nrow(occurrences) * boot_proportion,
                               replace = FALSE
                           ))
    boot.back <- replicate(n = boot_n,
                           sample(
                               x = seq_along(1:nrow(backgr)),
                               size = nrow(backgr) * boot_proportion,
                               replace = FALSE
                           ))
    boot_p <- matrix(data = 1,
                     nrow = nrow(occurrences),
                     ncol = boot_n,
                     dimnames = list(NULL, paste0("boot", 1:boot_n)))
    boot_a <- matrix(data = 1,
                     nrow = nrow(backgr),
                     ncol = boot_n,
                     dimnames = list(NULL, paste0("boot", 1:boot_n)))
    for (i in seq_along(1:boot_n)) {
        boot_p[, i][boot.pres[, i]] <- 0
        }
    for (i in seq_along(1:boot_n)) {
        boot_a[, i][boot.back[, i]] <- 0
    }
    boot.matrix <- rbind(boot_p, boot_a)
}

    if (exists("group.all"))   sdmdata <- data.frame(group.all, sdmdata)
    if (exists("cv.matrix"))   sdmdata <- data.frame(cv.matrix, sdmdata)
    if (exists("boot.matrix")) sdmdata <- data.frame(boot.matrix, sdmdata)
    message(paste("saving sdmdata", "\n"))
    write.table(sdmdata, file = paste(setup.folder, "sdmdata.csv", sep = "/"),
                sep = ",", row.names = FALSE, col.names = TRUE)


    if (plot_sdmdata) {
        message(paste("Plotting the dataset...", "\n"))
        png(filename = paste0(setup.folder, "/sdmdata_", species_name, ".png"))
        par(mfrow = c(1, 1), mar = c(5, 4, 3, 0))
        raster::plot(predictors[[1]], legend = FALSE, col = "grey90", colNA = NA)
        points(back, pch = ".", col = "black")
        points(pres, pch = 3, col = "grey50")
        legend("topleft", pch = c("+", "."),
               col = c("grey50", "black"), legend = c("Occs", "Back"))
        dev.off()
    }
    message("DONE!")
    return(sdmdata)
}
