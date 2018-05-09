#' Prepares the dataset to perform ENM
#'
#' @param species_name A character string with the species name
#' @param coordinates A data frame with occurrence data
#' @param real_absences User defined absence points
#' @param lon the name of the longitude column. defaults to "lon"
#' @param lat the name of the latitude column. defaults to "lat"
#' @param buffer Defines if a buffer will be used to sample pseudo-absences
#'        (F, "mean", "median", "max")
#' @param seed For reproducibility purposes
#' @param predictors A RasterStack of predictor variables
#' @param clean_dupl Logical, delete duplicate occurrence points? defaults to
#'  TRUE
#' @param clean_nas Logical, delete occurrence points with no environmental
#' information? Defaults to FALSE and can take a while for large datasets
#' @param geo_filt Logical, delete occurrence that are too close?
#' @param geo_filt_dist The distance of the geographic filter (in kilometers)
#' @param models_dir Folder path to save the output files
#' @param plot_sdmdata Logical, whether png files will be written
#' @param n_back Number of pseudoabsence points
#' @param bootstrap Logical, perform bootstrap
#' @param boot_proportion Numerical 0 to 1, proportion of points to be sampled
#' for bootstrap
#' @param boot_n How many bootstrap runs
#' @param crossvalidation Logical, perform k-fold crossvalidation
#' @param cv_partitions Number of partitions in the crossvalidation
#' @param cv_n How many crossvalidation runs
#' @return A dataframe called sdmdata with the groups for each run
#' (in columns called cv.1, cv.2 or boot.1, boot.2), a presence/absence vector,
#' the geographical coordinates, of the occurrence and pseudoabsence points, and
#' the associated environmental variables.
#' @author Andrea Sánchez-Tapia
#' @seealso \code{\link[dismo]{gridSample}}

#' @export
#'
#'
# tabela de valores
setup_sdmdata <- function(species_name = species_name,
                          coordinates = coordinates,
                          real_absences = NULL,
                          lon = "lon",
                          lat = "lat",
                          buffer = FALSE,
                          seed = 512,
                          predictors = predictors,
                          clean_dupl = T,
                          clean_nas = F,
                          geo_filt = F,
                          geo_filt_dist = NULL,
                          models_dir = models_dir,
                          plot_sdmdata = T,
                          n_back = 1000,
                          bootstrap = F,
                          boot_proportion = 0.8,
                          boot_n = 10,
                          crossvalidation = F,
                          cv_partitions = cv_partitions,
                          cv_n = 10) {
    if (file.exists(paste0(models_dir)) == FALSE)
        dir.create(paste0(models_dir))
    if (file.exists(paste0(models_dir, "/", species_name)) == FALSE)
        dir.create(paste0(models_dir, "/", species_name))
    partition.folder <- paste0(models_dir, "/", species_name, "/present", "/partitions")
    if (file.exists(partition.folder) == FALSE)
        dir.create(partition.folder, recursive = T)

    coordinates <- coordinates[,c(lon,lat)]
    names(coordinates) <- c("lon", "lat")
    # tabela de valores
    presvals <- raster::extract(predictors, coordinates)
    if (clean_dupl == TRUE) {
        coordinates <- coordinates[!base::duplicated(coordinates),]
        presvals <- presvals[!base::duplicated(coordinates),]
    }
    if (clean_nas == TRUE) {
        coordinates <- coordinates[complete.cases(presvals),]
        presvals <- presvals[complete.cases(presvals),]
    }
    if (geo_filt == TRUE) {
        coordinates <- geo_filt(coordinates = coordinates, min.distance = geo_filt_dist)
        presvals <- raster::extract(predictors, coordinates)
    }
    if (!is.null(real_absences)) {
        backgr <- real_absences[,c(lon, lat)]
    } else {
    if (buffer %in% c("mean", "max", "median")) {
        backgr <- create_buffer(coord = coordinates,
                                n_back = n_back,
                                buffer_type = buffer,
                                seed = seed,
                                predictors = predictors)

    } else {
        set.seed(seed + 2)
        backgr <- dismo::randomPoints(mask = predictors,
                                      n = n_back,
                                      p = coordinates,
                                      excludep = T)
    }
    }

    colnames(backgr) <- c("lon", "lat")

    # Extraindo dados ambientais dos bckgr
    backvals <- raster::extract(predictors, backgr)

    pa <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))

    pres <- cbind(coordinates, presvals)
    back <- cbind(backgr, backvals)
    coord_env_all <- rbind(pres, back)
    sdmdata <- cbind(pa,coord_env_all)
    # Data partition-----
    #Crossvalidation, repetated crossvalidation and jacknife
    if (crossvalidation == TRUE) {
        if (nrow(coordinates) < 11) {
            #forces jacknife
            cv_partitions <- nrow(coordinates)
            cv_n <- 1
        }
        if (cv_n == 1) {
            #Crossvalidation
            set.seed(seed)  #reproducibility
            group <- dismo::kfold(coordinates, cv_partitions)
            set.seed(seed + 1)
            bg.grp <- dismo::kfold(backgr, cv_partitions)
            group.all <- c(group, bg.grp)
        }
        if (cv_n > 1) {
            # Repeated CV
            cv.pres <- replicate(n = cv_n,
                                 dismo::kfold(coordinates, cv_partitions))
            dimnames(cv.pres) <- list(NULL, paste0("cv", 1:cv_n))
            cv.back <- replicate(n = cv_n,
                                 dismo::kfold(backgr, cv_partitions))
            dimnames(cv.back) <- list(NULL, paste0("cv", 1:cv_n))
            cv.matrix <- rbind(cv.pres, cv.back)
        }
    }
    # Bootstrap
    if (bootstrap == TRUE) {
    boot.pres <- replicate(n = boot_n,
                           sample(
                               x = seq_along(1:nrow(coordinates)),
                               size = nrow(coordinates) * boot_proportion,
                               replace = FALSE
                           ))
    boot.back <- replicate(n = boot_n,
                           sample(
                               x = seq_along(1:nrow(backgr)),
                               size = nrow(backgr) * boot_proportion,
                               replace = FALSE
                           ))
    boot_p <- matrix(data = 1,
                     nrow = nrow(coordinates),
                     ncol = boot_n,
                     dimnames = list(NULL,paste0("boot",1:boot_n)))
    boot_a <- matrix(data = 1,
                     nrow = nrow(backgr),
                     ncol = boot_n,
                     dimnames = list(NULL,paste0("boot",1:boot_n)))
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
    write.table(sdmdata, file = paste0(partition.folder, "/sdmdata.txt"))



    if (plot_sdmdata) {
        #Creates a .png plot of the initial dataset
        cat(paste("Plotting the dataset...",'\n'))
        png(filename = paste0(partition.folder, "/sdmdata_", species_name,".png"))
        par(mfrow = c(1, 1), mar = c(5, 4, 3, 0))
        raster::plot(predictors[[1]], legend = F, col = "grey90", colNA = NA)
        points(back, pch = ".", col = "black")
        points(pres, pch = 3, col = "grey50")
        legend("topleft", pch = c("+","."),
               col = c("grey50", "black"), legend = c("Occs","Back"))
        dev.off()
    }
    rm(coord_env_all)
    rm(pres)
    rm(back)
    gc()
    return(sdmdata)
}
