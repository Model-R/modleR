#' Creates ensemble models from several algorithms
#'
#' This function reads the output of \code{final_model} for each species and
#' multiple algorithms and builds a simple ensemble model by calculating the
#' mean of the final models in order to obtain one model per species. It also
#' calculates median, standard deviation and range (maximum - minimum)
#'
#' @inheritParams setup_sdmdata
#' @inheritParams final_model
#' @param ensemble_dir Character string, name of the folder to save the output
#' files. A subfolder will be created. Defaults to "\code{ensemble}"
#' @param algorithms Character vector specifying which algorithms will be
#' processed. Note that it can have length > 1, ex. \code{c("bioclim", "rf")}.
#' Defaults to NULL: it no name is given it will process all algorithms present
#' in the final_models folder
#' @param which_ensemble Which method to apply consensus between algorithms will
#' be used? Current options are:
#' \describe{
#'   \item{\code{best}}{Selects models from the best-performing algorithm. A
#'     performance metric must be specified (\code{performance_metric}). Parameter
#'     \code{which_final} indicates which model will be returned}
#'   \item{\code{average}}{Computes the means between models. Parameter
#'     \code{which_final} indicates which model will be returned}
#'   \item{\code{weighted_average}}{Computes a weighted mean between models. A
#'     performance metric must be specified. Parameter \code{which_final} indicates
#'     which model will be returned}
#'   \item{\code{median}}{Computes the median between models. Parameter
#'     \code{which_final} indicates which model will be returned}
#'   \item{\code{frequency}}{Computes the mean between binary models, which is
#'     analogous to calculating a relative consensus}
#'   \item{\code{consensus}}{Computes a binary model with the final consensus area.
#'      A \code{consensus_level} must be specified}
#'   \item{\code{pca}}{Computes a PCA between the models for each algorithm and extract the first axis, that summarizes variation between them}
#' }
#' @param which_final Which \code{final_model} will be used to calculate the
#' average, weighted average or median ensembles? See \code{\link{final_model}}
#' @param performance_metric Which performance metric will be used to define
#' the \code{"best"} algorithm any in \code{c("AUC", "pROC", "TSSmax",
#' "KAPPAmax", "CCR", "F_score", "Jaccard")}
#' @param dismo_threshold Character string indicating threshold (cut-off) to
#' transform raw_mean final models to binary for frequency and consensus methods.
#' The options are from \code{\link[dismo]{threshold}}:
#'  "\code{kappa}", "\code{spec_sens}", "\code{no_omission}",
#'   "\code{prevalence}", "\code{equal_sens_spec}",
#'  "\code{sensitivity}". Default value is "\code{spec_sens}"
#' @param png_ensemble Logical. If \code{TRUE} writes png files of the
#' ensemble models
#' @param write_map Logical. If \code{TRUE} adds a map contour to the png file
#' of the ensemble models
#' @param write_occs Logical. If \code{TRUE} writes the occurrence points on the
#' png file of the ensemble model
#' @param ... Other parameters from \code{\link[raster]{writeRaster}}
#' @param uncertainty Calculates the uncertainty between models, as a range
#' (maximum - minimum)
#' @import raster
#' @importFrom scales alpha
#' @import graphics
#' @importFrom stats sd
#' @importFrom stats prcomp
#' @export
#' @seealso \code{\link{final_model}}
#' @return Retuns a RasterStack with all generated statistics written in the
#' \code{ensemble_dir} subfolder
#' @return Writes on disk raster files with the median, mean and standard
#' deviation and range of the assembled models
#' @return If \code{png_ensemble = TRUE} writes .png figures
#'  in the \code{ensemble_dir} subfolder
#' @examples
#' # run setup_sdmdata
#' sp <- names(example_occs)[1]
#' sp_coord <- example_occs[[1]]
#' sp_setup <- setup_sdmdata(species_name = sp,
#'                           occurrences = sp_coord,
#'                           predictors = example_vars)
#'
#' # run do_many
#' sp_many <- do_many(species_name = sp,
#'                    predictors = example_vars,
#'                    bioclim = TRUE,
#'                    maxnet = TRUE)
#'
#' # run final_model
#' sp_final <- final_model(species_name = sp,
#'                         algorithms = c("bioclim", "maxnet"),
#'                         select_partitions = TRUE,
#'                         select_par = "TSSmax",
#'                         select_par_val = 0,
#'                         which_models = c("raw_mean"),
#'                         consensus_level = 0.5,
#'                         overwrite = TRUE)
#'
#' # run ensemble model
#' sp_ensemble <- ensemble_model(species_name = sp,
#'                               occurrences = sp_coord,
#'                               overwrite = TRUE)
ensemble_model <- function(species_name,
                           occurrences,
                           lon = "lon",
                           lat = "lat",
                           models_dir = "./models",
                           final_dir = "final_models",
                           ensemble_dir = "ensemble",
                           proj_dir = "present",
                           algorithms = NULL,
                           which_ensemble = c("average"),
                           which_final = c("raw_mean"),
                           performance_metric = "TSSmax",
                           dismo_threshold = "spec_sens",
                           consensus_level = 0.5,
                           png_ensemble = TRUE,
                           write_occs = FALSE,
                           write_map = FALSE,
                           scale_models = TRUE,
                           uncertainty = TRUE,
                           ...) {


    ## output folder
    if (file.exists(
        paste0(models_dir, "/", species_name, "/", proj_dir, "/",
               ensemble_dir, "/")) == FALSE) {
        dir.create(paste0(models_dir, "/", species_name, "/", proj_dir, "/",
                          ensemble_dir, "/"))
    }
    print(date())
    cat(paste(species_name, "\n"))
    cat(paste("Reading mean evaluation files for", species_name, "in", proj_dir, "\n"))

    ## reads mean statistics table per algorithm

    stats_summary <- read.csv(paste0(models_dir, "/", species_name, "/present/",
                                     final_dir, "/", species_name,
                                     "_mean_statistics.csv"), header = T,
                              row.names = 1)
    #specify algorithms or not
    if (is.null(algorithms)) {
        algorithms <- unique(stats_summary$algorithm)
    }
    stats_summary <- stats_summary[, stats_summary$algorithm %in% algorithms]
    ensemble_mods <- raster::stack()

    if ("best" %in% which_ensemble) {
        if (is.null(performance_metric))
        stop("A performance metric must be specified to select the 'best' algorithm")
    best <- which.max(stats_summary[,performance_metric])
    message(paste("The best performing algorithm was",
                  stats_summary$algorithm[best], "according to",
                  performance_metric, "values"))
    if (is.null(which_final))
        stop("A which_final model must be specified to return the 'best' algorithm")
    best_mod_files <- list.files(paste0(models_dir, "/",
                                        species_name, "/",
                                        proj_dir, "/",
                                        final_dir),
                                 recursive = TRUE,
                                 full.names = TRUE,
                                 pattern =
                                     paste0(stats_summary$algorithm[best],
                                            "_", which_final, ".tif$"))
    if (length(best_mod_files) == 0)
        stop(paste("No", which_final, "models to ensemble from for", species_name, "\n"))
    best_mod <- raster(best_mod_files)
    names(best_mod) <- "best"
    writeRaster(best_mod,
                filename = paste0(models_dir, "/", species_name,
                                  "/", proj_dir, "/",
                                  ensemble_dir, "/", species_name,
                                  "_", which_final,
                                  "_ensemble", "_best",
                                  ".tif"),
                ...
                )
    ensemble_mods <- raster::addLayer(ensemble_mods, best_mod)

    }
    if ("average" %in% which_ensemble) {
        raw_mean_files <- list.files(paste0(models_dir, "/",
                                            species_name, "/",
                                            proj_dir, "/",
                                            final_dir),
                                     recursive = TRUE,
                                     full.names = TRUE,
                                     pattern =
                                         paste0(
                                                "_raw_mean.tif$"))
        if (length(raw_mean_files) == 0)
            stop(paste("No models to assemble from for", species_name, "\n"))
        raw_mean_models <- raster::stack(raw_mean_files)
        average_ensemble <- mean(raw_mean_models)
        names(average_ensemble) <-  "average"
        writeRaster(average_ensemble,
                    filename = paste0(models_dir, "/", species_name,
                                      "/", proj_dir, "/",
                                      ensemble_dir, "/", species_name,
                                      "_ensemble", "_average",
                                      ".tif"),
                    ...
        )
        ensemble_mods <- raster::addLayer(ensemble_mods, average_ensemble)

    }
    if ("weighted_average" %in% which_ensemble) {
        if (is.null(performance_metric))
            stop("A performance metric must be specified to compute the weighted average")
        w_coefs <- stats_summary[,performance_metric]
        raw_mean_files <- list.files(paste0(models_dir, "/",
                                            species_name, "/",
                                            proj_dir, "/",
                                            final_dir),
                                     recursive = TRUE,
                                     full.names = TRUE,
                                     pattern =
                                         paste0(
                                                "_raw_mean.tif$"))
        if (length(raw_mean_files) == 0)
            stop(paste("No models to ensemble from for", species_name, "\n"))
        raw_mean_models <- raster::stack(raw_mean_files)
        weighted_average <- raster::weighted.mean(raw_mean_models, w_coefs)
        names(weighted_average) <- "weighted_average"
        writeRaster(weighted_average,
                    filename = paste0(models_dir, "/", species_name,
                                      "/", proj_dir, "/",
                                      ensemble_dir, "/", species_name,
                                      "_", performance_metric,
                                      "_ensemble_weighted_average",
                                      ".tif"),
                    ...
        )
        ensemble_mods <- raster::addLayer(ensemble_mods, weighted_average)

    }
    if ("median" %in% which_ensemble) {
        raw_mean_files <- list.files(paste0(models_dir, "/",
                                            species_name, "/",
                                            proj_dir, "/",
                                            final_dir),
                                     recursive = TRUE,
                                     full.names = TRUE,
                                     pattern =
                                         paste0(
                                             "_raw_mean.tif$"))
        if (length(raw_mean_files) == 0)
            stop(paste("No models to ensemble from for", species_name, "\n"))
        raw_mean_models <- raster::stack(raw_mean_files)
        median_ensemble <- raster::calc(raw_mean_models,
                                         fun = function(x) {
                                             stats::median(x, na.rm = TRUE)
                                             }
                                         )
        names(median_ensemble) <- "median"
        writeRaster(median_ensemble,
                    filename = paste0(models_dir, "/", species_name,
                                      "/", proj_dir, "/",
                                      ensemble_dir, "/", species_name,
                                      "_ensemble", "_median",
                                      ".tif"),
                    ...
        )
        ensemble_mods <- raster::addLayer(ensemble_mods, median_ensemble)
    }
    if (any(c("frequency", "consensus") %in% which_ensemble)) {
        #reads raw
        raw_mean_files <- list.files(paste0(models_dir, "/",
                                            species_name, "/",
                                            proj_dir, "/",
                                            final_dir),
                                     recursive = TRUE,
                                     full.names = TRUE,
                                     pattern =
                                         paste0(
                                             "_raw_mean.tif$"))
        if (length(raw_mean_files) == 0)
            stop(paste("No models to ensemble from for", species_name, "\n"))
        raw_mean_models <- raster::stack(raw_mean_files)
        if (is.null(dismo_threshold))
            stop("A dismo_threshold must be specified to create binary models")
        #cuts the models by the mean threshold for each algorithm
        th <- stats_summary[,dismo_threshold]
        bin_mean_models <- raw_mean_models > th
        #calculates the mean
        frequency_ensemble <- mean(bin_mean_models)
        names(frequency_ensemble) <- "frequency"
        if ("frequency" %in% which_ensemble) {
        writeRaster(frequency_ensemble,
                    filename = paste0(models_dir, "/", species_name,
                                      "/", proj_dir, "/",
                                      ensemble_dir, "/", species_name,
                                      "_ensemble", "_frequency",
                                      ".tif"),
                    ...
        )
            ensemble_mods <- raster::addLayer(ensemble_mods, frequency_ensemble)
        }
        if ("consensus" %in% which_ensemble) {
            if (missing("consensus_level"))
                stop("Parameter consensus_level must be specified to calculate consensus average")
            consensus_ensemble <- frequency_ensemble > consensus_level
            names(consensus_ensemble) <- "consensus"
            writeRaster(consensus_ensemble,
                        filename = paste0(models_dir, "/", species_name,
                                          "/", proj_dir, "/",
                                          ensemble_dir, "/", species_name,
                                          "_ensemble_", consensus_level,
                                          "_consensus",
                                          ".tif"),
                        ...
            )
            ensemble_mods <- raster::addLayer(ensemble_mods, consensus_ensemble)

        }
    }
    if ("pca" %in% which_ensemble) {
        raw_mean_files <- list.files(paste0(models_dir, "/",
                                            species_name, "/",
                                            proj_dir, "/",
                                            final_dir),
                                     recursive = TRUE,
                                     full.names = TRUE,
                                     pattern =
                                         paste0(
                                             "_raw_mean.tif$"))
        if (length(raw_mean_files) == 0)
            stop(paste("No models to ensemble from for", species_name, "\n"))
        raw_mean_models <- raster::stack(raw_mean_files)
        vals <- raster::getValues(raw_mean_models)
        vals <- vals[!is.na(rowSums(vals)),]
        vals.st <- scale(vals)
        pca_mod <- prcomp(vals.st)
        summary_pca <- summary(pca_mod)
        #axis_nb <- which(summary_pca$importance["Cumulative Proportion",] >= 0.95)[1]
        expl <- summary_pca$importance["Cumulative Proportion",1]
        first_axis <- predict(raw_mean_models, pca_mod, index = 1)
        first_axis <- rescale_layer(first_axis)
        names(first_axis) <- "pca"
        writeRaster(first_axis,
                    filename = paste0(models_dir, "/", species_name,
                                      "/", proj_dir, "/",
                                      ensemble_dir, "/", species_name,
                                      "_ensemble_pca_", round(expl,3),
                                      ".tif"),
                    ...
        )
        ensemble_mods <- raster::addLayer(ensemble_mods, first_axis)
    }
    if (uncertainty == TRUE) {
        raw_mean_files <- list.files(paste0(models_dir, "/",
                                            species_name, "/",
                                            proj_dir, "/",
                                            final_dir),
                                     recursive = TRUE,
                                     full.names = TRUE,
                                     pattern =
                                         paste0(
                                             "_raw_mean.tif$"))
        raw_mean_models <- raster::stack(raw_mean_files)
        message("Calculating range")
        ensemble_inctz <- raster::calc(raw_mean_models,
                                       fun = function(x) {
                                           max(x) - min(x)
                                           }
                                       )
        names(ensemble_inctz) <- "uncertainty"
        writeRaster(ensemble_inctz,
                    filename = paste0(models_dir, "/", species_name,
                                      "/", proj_dir, "/",
                                      ensemble_dir, "/", species_name,
                                      "_uncertainty.tif"),
                    ...
        )
        ensemble_mods <- raster::addLayer(ensemble_mods, ensemble_inctz)

    }

            #    #scale models to 0-1#ö escalar??
        #    if (scale_models == TRUE) {
        #        mod2 <- rescale_layer(mod2)
        #    }


            coord <- occurrences[, c(lon, lat)]

            if (png_ensemble) {
            message("Writing pngs")
                for (i in 1:dim(ensemble_mods)[3]) {
                png(filename = paste0(models_dir, "/",
                                      species_name, "/",
                                      proj_dir, "/",
                                      ensemble_dir, "/",
                                      species_name, "_",
                                      names(ensemble_mods)[i], ".png"),
                    res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
                raster::plot(ensemble_mods[[i]])
                if (write_map) {
                maps::map("world",
                          add = TRUE,
                          col = "grey")
                }
                if (write_occs) {
                points(coord, pch = 21, cex = 0.6,
                       bg = scales::alpha("cyan", 0.6))
                }
                dev.off()
                }
            }
    # creating and writing ensemble_model metadata
    metadata <- data.frame(
        species_name = as.character(species_name),
        algorithms = paste(algorithms, collapse = "-"),
        which_ensemble = paste(which_ensemble, collapse = "-"),
        which_final = paste(which_final, collapse = "-"),
        performance_metric = as.character(performance_metric),
        dismo_threshold = as.character(dismo_threshold),
        consensus_level = ifelse("consensus" %in% which_ensemble, consensus_level, NA),
        scale_models = ifelse(scale_models, "yes", "no"),
        uncertainty = ifelse(uncertainty, "yes", "no")
    )
    message("writing metadata")
    write.csv(metadata, file = paste0(models_dir, "/", species_name, "/",
                                      proj_dir, "/", ensemble_dir,
                                      "/metadata.csv"))
    print("DONE!")
    print(date())
    return(ensemble_mods)
}
