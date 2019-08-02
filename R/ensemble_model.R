#' Create ensemble models from several algorithms.
#'
#' \code{ensemble} takes the input of \code{final_model} and builds an ensemble
#'  model by calculating the mean of each final model per algorithm.
#'
#' @inheritParams final_model
#' @param species_name A character string with the species name
#' @param occurrences A two-column data frame with the occurrence points
#' @param models_dir Character string. Folder path where the input files
#'                   are located
#' @param final_dir Character string, name of the subfolder where the files
#'                  for the final models are located
#' @param ensemble_dir Character string, name of the folder to save the output
#'                     files. A subfolder will be created.
#' @param proj_dir Character. The name of the subfolder with the projection.
#' Defaults to "present" but can be set according to the other projections (i.e.
#' to execute the function in projected models)
#'
#' @param which_final Which final_model() will be used? Currently it can be:
#' \describe{
#'   \item{\code{weighted_AUC} or \code{weighted_TSS}}{the models weighted
#'   by TSS or AUC}
#'   \item{\code{raw_mean}}{the mean of the selected raw models}
#'   \item{\code{bin_mean_th}}{the binary model created by cutting \code{raw_mean} by the mean of the thresholds that
#'    maximize the selected evaluation metric (e.g. TSS (\code{spec_sens}) or other dismo thresholds)}
#'    \item{\code{cut_mean_th}}{the cut model created by recovering \code{raw_mean} values above the mean threshold that
#'    maximizes the selected evaluation metric (e.g. TSS (\code{spec_sens}) or other dismo thresholds)}
#'   \item{\code{bin_mean}}{the mean of the selected binary models}
#'   \item{\code{bin_consensus}}{the binary consensus from \code{bin_mean}}.
#'   \item{\code{cut_mean}}{the mean of the selected cut models}
#' }
#' @param consensus Logical. Will a consensus between the algorithms be applied?
#' @param write_ensemble Logical. If \code{TRUE} writes png files of the ensemble models.
#' @param lon the name of the longitude column. defaults to "lon"
#' @param lat the name of the latitude column. defaults to "lat"
#' @param ... Other parameters from \code{\link[raster]{writeRaster}}
#'
#' @import raster
#' @importFrom scales alpha
#' @import graphics
#' @importFrom stats sd
#' @export
#' @seealso \link{final_model}
#' @return A rasterStack with the minimun, maximum, median, mean and standard
#' deviation of the assembled models. A set of ensemble models and figures
#' (optional) written in the \code{ensemble_dir} subfolder
ensemble_model <- function(species_name,
                           occurrences,
                           lon = "lon",
                           lat = "lat",
                           models_dir = "./models",
                           final_dir = "final_models",
                           ensemble_dir = "ensemble",
                           proj_dir = "present",
                           which_final = c("raw_mean"),
                           consensus = FALSE,
                           consensus_level = 0.5,
                           write_ensemble = T,
                           scale_models = TRUE, 
                           ...) {


    ## output folder
    if (file.exists(
        paste0(models_dir, "/", species_name, "/", proj_dir, "/", ensemble_dir, "/")) == FALSE) {
        dir.create(paste0(models_dir, "/", species_name, "/", proj_dir, "/", ensemble_dir, "/"))
    }

    ## for each model specified in final_models
    for (whi in which_final) {
        cat(paste(whi, "-", species_name, "\n"))  #lê os arquivos
        tif.files <- list.files(paste0(models_dir, "/", species_name, "/", proj_dir, "/",
                                       final_dir), recursive = T,
                                full.names = T, pattern = paste0(whi, ".*tif$"))

        if (length(tif.files) == 0) {
            cat(paste("No", whi, "models to ensemble from for", species_name, "\n"))
        } else {
            cat(paste(length(tif.files), whi,
                      "models to ensemble from for", species_name, "\n"))
            mod2 <- raster::stack(tif.files)
            #scale models to 0-1
            if (scale_models == T) {
                mod2 <- rescale_layer(mod2)
            }
            message("Calculating mean")
            ensemble.mean <- raster::calc(mod2,
                                          fun = function(x) {
                                              mean(x, na.rm = T)
                                          }
            )
            message("Calculating sd")
            ensemble.sd <- raster::calc(mod2,
                                        fun = function(x) {
                                            sd(x, na.rm = T)
                                        }
            )
            # message("Calculating min")
            # ensemble.min <- raster::calc(mod2,
            #                              fun = function(x) {
            #                                  min(x, na.rm = T)
            #                              }
            # )
            # message("Calculating max")
            # ensemble.max <- raster::calc(mod2,
            #                              fun = function(x) {
            #                                  max(x, na.rm = T)
            #                              }
            # )
            message("Calculating median")
            ensemble.median <- raster::calc(mod2,
                                            fun = function(x) {
                                                stats::median(x, na.rm = T)
                                            }
            )
            message("Calculating range")
            ensemble.inctz <- raster::calc(mod2,
                                           fun = function(x) {
                                               max(x) - min(x)
                                           }
            )
            message("Stack results")
            ensemble.mods <- raster::stack(ensemble.mean,
                                           ensemble.median,
                                           ensemble.sd,
                                           #ensemble.min,
                                           #ensemble.max,
                                           ensemble.inctz
                                           )
            names(ensemble.mods) <- c("mean",
                                      "median",
                                      "sd",
                                      #"min",
                                      #"max",
                                      "range"
                                      )

            coord <- occurrences[, c(lon, lat)]

            if (write_ensemble) {
            message("Writing pngs")
                for (i in 1:dim(ensemble.mods)[3]) {
                png(filename = paste0(models_dir, "/", species_name, "/", proj_dir, "/",
                                      ensemble_dir, "/", species_name, "_", whi,
                                      "_", names(ensemble.mods)[i], ".png"),
                    res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
                raster::plot(ensemble.mods[[i]])
                maps::map("world",
                          c("", "South America"),
                          add = T,
                          col = "grey")
                points(coord, pch = 21, cex = 0.6,
                       bg = scales::alpha("cyan", 0.6))
                dev.off()
                }
            }

            raster::writeRaster(ensemble.mods,
                                filename = paste0(models_dir, "/", species_name,
                                                  "/", proj_dir, "/",
                                                  ensemble_dir, "/", species_name, "_",
                                                  whi),
                                bylayer = T,
                                suffix = "names",
                                format = "GTiff", 
                                ...)

            #### Consensus models
            if (consensus == TRUE) {
                ensemble.consensus <- ensemble.mean >= consensus_level
                raster::writeRaster(ensemble.consensus,
                                    filename = paste0(models_dir, "/", species_name,
                                                      "/", proj_dir, "/",
                                                      ensemble_dir, "/", species_name,
                                                      "_", whi,
                                                      "_ensemble", "_meanconsensus",
                                                      consensus_level * 100,
                                                      ".tif"), 
                                    ...)

                if (write_ensemble) {
                png(filename = paste0(models_dir, "/", species_name, "/", proj_dir, "/",
                                      ensemble_dir, "/",
                                      species_name, "_", whi,
                                      "_ensemble", "_meanconsensus",
                                      consensus_level * 100, ".png"), res = 300,
                    width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
                raster::plot(ensemble.consensus)
                maps::map("world", , add = T, col = "grey")
                points(coord, pch = 19, cex = 0.3,
                       col = scales::alpha("cyan", 0.6))
                dev.off()
                }
            }
        }
    }
    print("DONE!")
    print(date())
    #return(ensemble.mods)
}
