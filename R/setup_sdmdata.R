#' Prepares the dataset to perform ENM.
#'
#' This function takes the occurrence points files and makes the data cleaning,
#' data partitioning and the pseudo-absence point sampling, and saves the
#' metadata and sdmdata files into the hard disk.
#'
#' @inheritParams create_buffer
#' @inheritParams rescale_layer
#' @inheritParams clean
#' @inheritParams geo_filt
#' @param species_name A character string with the species name
#' @param occurrences A data frame with occurrence data
#' @param lon the name of the longitude column. defaults to "lon"
#' @param lat the name of the latitude column. defaults to "lat"
#' @param predictors A RasterStack of predictor variables
#' @param seed For reproducibility purposes
#' @param real_absences User-defined absence points
#' @param geo_filt Logical, delete occurrence that are too close?
#' @param geo_filt_dist The distance of the geographic filter (in kilometers)
#' @param models_dir Folder path to save the output files
#' @param plot_sdmdata Logical, whether png files will be written
#' @param n_back Number of pseudoabsence points
#' @param partition_type Perform bootstrap or k-fold crossvalidation?
#' @param boot_proportion Numerical 0 to 1, proportion of points to be sampled
#' for bootstrap
#' @param boot_n How many bootstrap runs
#' @param cv_partitions Number of partitions in the crossvalidation
#' @param cv_n How many crossvalidation runs
#' @param equalize Logical, whether the number of presences and absences should be
#' equalized in randomForest and brt.
#' @param ... parameters from create_buffer()
#' @return A dataframe called sdmdata with the groups for each run
#' (in columns called cv.1, cv.2 or boot.1, boot.2), a presence/absence vector,
#' the geographical coordinates, of the occurrence and pseudoabsence points, and
#' the associated environmental variables.
#' @author Andrea Sánchez-Tapia
#' @examples
#' predictors <- example_vars
#' species = unique(coordenadas$sp)
#' setup_sdmdata(species_name = species[1], occurrences = coordenadas[-1], predictors = predictors, models_dir = './models_dir',
#' real_absences = NULL, seed = 55, clean_dupl = T, clean_nas = T, clean_uni = T, partition_type = c("crossvalidation"),
#' cv_n = 1, cv_partitions = 3, equalize = T)
#'
#' @seealso \code{\link[dismo]{gridSample}}

#' @export
#'
#'
# tabela de valores
setup_sdmdata <- function(species_name = species_name,
                          occurrences = occurrences,
                          predictors = predictors,
                          models_dir = models_dir,
                          real_absences = NULL,
                          lon = "lon",
                          lat = "lat",
                          buffer_type = NULL,
                          dist_buf = NULL,
                          buffer_shape = NULL,
                          seed = NULL,
                          clean_dupl = T,
                          clean_nas = F,
                          clean_uni = T,
                          geo_filt = F,
                          geo_filt_dist = NULL,
                          select_variables = F,
                          percent_correlation = 0.8,
                          plot_sdmdata = T,
                          n_back = 1000,
                          partition_type = c("bootstrap"),
                          boot_n = 1,
                          boot_proportion = 0.7,
                          cv_n = NULL,
                          cv_partitions = NULL,
                          equalize = NULL,
                          ...) {
    if (file.exists(paste0(models_dir)) == FALSE)
        dir.create(paste0(models_dir), recursive = T, showWarnings = F)
    if (file.exists(paste0(models_dir, "/", species_name)) == FALSE)
        dir.create(paste0(models_dir, "/", species_name))
    partition.folder <-
        paste0(models_dir, "/", species_name, "/present", "/partitions")
    if (file.exists(partition.folder) == FALSE)
        dir.create(partition.folder, recursive = T)

    ## checking latitude and longitude columns
    if (all(c(lon, lat) %in% names(occurrences))) {
    occurrences <- occurrences[, c(lon, lat)]
    names(occurrences) <- c("lon", "lat")
    } else {
        stop("Coordinate column names do not match. Either rename to `lon` and `lat` or specify")
    }
    #writes some original metadata that will be modified later
    original_n <- nrow(occurrences)
    original_n_back <- n_back
    original_predictors <- paste(names(predictors), collapse = '-')

        #checking metadata
    if (file.exists(paste0(partition.folder, "/metadata.txt"))) {
        message("metadata file found, checking metadata \n")
        metadata_old <- read.table(paste0(partition.folder, "/metadata.txt"), as.is = F, row.names = 1)
        metadata_old <- metadata_old[, setdiff(names(metadata_old), c("final.n", "final.n.back", "selected_predictors"))]
        metadata_new <- data.frame(
            species_name = as.character(species_name),
            original_predictors = original_predictors,
            select_variables = select_variables,
            original.n = original_n,
            original.n.back = original_n_back,
            buffer_type = ifelse(is.null(buffer_type), NA, buffer_type),
            dist_buf = ifelse(is.null(dist_buf), NA, dist_buf),
            seed = ifelse(is.null(seed), NA, seed),
            res.x = res(predictors)[1],
            res.y = res(predictors)[2],
            clean_dupl = clean_dupl,
            clean_nas = clean_nas,
            clean_uni = clean_uni,
            geo_filt = geo_filt,
            geo_filt_dist = ifelse(is.null(geo_filt_dist), NA, geo_filt_dist),
            models_dir = models_dir,
            partition = partition_type,
            boot_proportion = ifelse(is.null(boot_proportion), NA, boot_proportion),
            boot_n = ifelse(is.null(boot_n), NA, boot_n),
            cv_partitions = ifelse(is.null(cv_partitions), NA, cv_partitions),
            cv_n = ifelse(is.null(cv_n), NA, cv_n),
            equalize = ifelse(is.null(equalize), NA, equalize),
            row.names = 1
            )

            if (all(all.equal(metadata_old, metadata_new) == T)) {
            message("same metadata, no need to run data partition")
            sdmdata <- read.table(paste0(partition.folder, "/sdmdata.txt"))
            return(sdmdata)
            }
    }
    message("performing data partition")

    occurrences <-
        clean(occurrences,
              lon = lon,
              lat = lat,
              predictors,
              clean_dupl = clean_dupl,
              clean_nas = clean_nas,
              clean_uni = clean_uni)

    if (geo_filt == TRUE) {
        message("applying a geographical filter")
        occurrences <-
            geo_filt(occurrences = occurrences, min_distance = geo_filt_dist)
    }
    final_n <- nrow(occurrences)


    #background selection:
    #first option: user-supplied background
    if (!is.null(real_absences)) {
        backgr <- real_absences[, c(lon, lat)]
    } else {
        #second option: there is a buffer
        if (!is.null(buffer_type)) {
            if (buffer_type %in% c("mean", "max", "median", "distance", "user")) {
                message("creating buffer")
                pbuffr <- create_buffer(occurrences = occurrences,
                                        models_dir = models_dir,
                                        species_name = species_name,
                                        buffer_type = buffer_type,
                                        predictors = predictors,
                                        dist_buf = dist_buf,
                                        dist_min = dist_min,
                                        ...)
            }
            # third option: there is no buffer
        } else {
            buffer_type <- "no"
            pbuffr <- predictors[[1]]
        }
        #before sampling pseudoabsence points
                #checks if there will be enough cells to sample pseudoabsences from
                vals <- values(pbuffr)
                available_cells <- sum(!is.na(vals)) - nrow(occurrences)
                # and corrects accordingly
                if (available_cells < n_back) {
                    n_back_mod <- available_cells
                    message(paste0(available_cells, "available cells"))
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
                                              excludep = T)
    colnames(backgr) <- c("lon", "lat")
    }

    # Seleccionando variables if sel_vars ==T
    if (select_variables == T) {
        predictors <- select_variables(species_name = species_name,
                                       predictors = predictors,
                                       models_dir = models_dir,
                                       buffer = pbuffr,
                                       cutoff = cutoff,
                                       percent = percent_correlation)

        }
    # tabela de valores
    message("extracting environmental data")
    presvals <- raster::extract(predictors, occurrences)
    # Extraindo dados ambientais dos bckgr
    message("extracting background data")
    backvals <- raster::extract(predictors, backgr)
    if (any(complete.cases(backvals) == F)) {
        backvals <- backvals[complete.cases(backvals), ]
        backgr   <- backgr[complete.cases(backvals), ]
        warning(paste("Your background data had NA values, ", nrow(backvals),
                      "points were retained"))
        }

    pa <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))

    pres <- cbind(occurrences, presvals)
    back <- cbind(backgr, backvals)
    coord_env_all <- rbind(pres, back)
    sdmdata <- cbind(pa, coord_env_all)
    # Data partition-----
    #Crossvalidation, repetated crossvalidation and jacknife
    #if (crossvalidation == TRUE) {
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
    #if (bootstrap == TRUE) {
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
    message("saving sdmdata")
    write.table(sdmdata, file = paste0(partition.folder, "/sdmdata.txt"))



    if (plot_sdmdata) {
        message("Plotting the dataset...", '\n')
        png(filename = paste0(partition.folder, "/sdmdata_", species_name, ".png"))
        par(mfrow = c(1, 1), mar = c(5, 4, 3, 0))
        raster::plot(predictors[[1]], legend = F, col = "grey90", colNA = NA)
        points(back, pch = ".", col = "black")
        points(pres, pch = 3, col = "grey50")
        legend("topleft", pch = c("+", "."),
               col = c("grey50", "black"), legend = c("Occs", "Back"))
        dev.off()
    }

    #metadata
    metadata <- data.frame(
        species_name = as.character(species_name),
        select_variables = select_variables,
        selected_predictors = paste(names(predictors), collapse = '-'),
        original.n = original_n,
        final.n = final_n,
        original.n.back = original_n_back,
        final.n.back = n_back_mod,
        buffer_type = ifelse(is.null(buffer_type), NA, buffer_type),
        dist_buf = ifelse(is.null(dist_buf), NA, dist_buf),
        seed = ifelse(is.null(seed), NA, seed),
        res.x = res(predictors)[1],
        res.y = res(predictors)[2],
        clean_dupl = clean_dupl,
        clean_nas = clean_nas,
        clean_uni = clean_uni,
        geo_filt = geo_filt,
        geo_filt_dist = ifelse(is.null(geo_filt_dist), NA, geo_filt_dist),
        models_dir = models_dir,
        partition = partition_type,
        boot_proportion = ifelse(is.null(boot_proportion), NA, boot_proportion),
        boot_n = ifelse(is.null(boot_n), NA, boot_n),
        cv_partitions = ifelse(is.null(cv_partitions), NA, cv_partitions),
        cv_n = ifelse(is.null(cv_n), NA, cv_n),
        equalize = ifelse(is.null(equalize), NA, equalize)
    )
    message("saving metadata")
    write.table(metadata, file = paste0(partition.folder, "/metadata.txt"))


    rm(coord_env_all)
    rm(pres)
    rm(back)
    gc()
    return(sdmdata)
}
