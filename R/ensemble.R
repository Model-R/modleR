#' Create ensemble models from several algorithms.
#'
#' \code{ensemble} takes the input of \code{final_model} and builds an ensemble
#'  model by calculating the mean of each final model per algorithm.
#'
#' @param species.name A character string with the species name
#' @param occs A two-column data frame with the occurrence points
#' @param models.dir Character string. Folder path where the input files
#'                   are located
#' @param final.dir Character string, name of the subfolder where the files
#'                  for the final models are located
#' @param ensemble.dir Character string, name of the folder to save the output
#'                     files. A subfolder will be created.
#' @param which.models Which final_model() will be used? Currently it can be:
#' \describe{
#'   \item{\code{weighted_AUC} or \code{weighted_TSS}}{the models weighted
#'   by TSS or AUC}
#'   \item{\code{final_model_3}}{the binary model created by selecting or not the
#'    partitions, taking their mean and cutting by the mean threshold that
#'    maximizes TSS (or other dismo thresholds)}
#'   \item{\code{final_model_7}}{the mean of the selected binary models}
#'   \item{\code{final_model_8}}{the binary consensus from \code{final_model_7}}
#' }
#' @param consensus Logical. Will a consensus rule be applied?
#' @param consensus.level Threshold for the consensus rule, betwen 0 and 1
#'                        (0.5 means a majority rule).
#' @param write_png Write png? Defaults to TRUE
#' @param write_raw_map Create a mean raw map without margins
#'
#' @import raster
#' @importFrom scales alpha
#' @import graphics
#' @importFrom stats sd
#' @export
#' @seealso \link{final_model}
#' @return A rasterStack with the minimun, maximum, median, mean and standard
#' deviation of the assembled models. A set of ensemble models and figures
#' (optional) written in the \code{ensemble.dir} subfolder
ensemble_model <- function(species.name,
                     occs,
                     models.dir = "./models",
                     final.dir = "final_models",
                     ensemble.dir = "ensemble",
                     which.models = c("final_model_3", "final_model_7"),
                     consensus = FALSE,
                     consensus.level = 0.5,
                     write_png = T,
                     write_raw_map = F) {

    ## output folder
    if (file.exists(
        paste0(models.dir, "/", species.name, "/present/", ensemble.dir, "/")) == FALSE) {
        dir.create(paste0(models.dir, "/", species.name, "/present/", ensemble.dir, "/"))
    }

    ## for each model specified in final_models
    for (whi in which.models) {
        cat(paste(whi, "-", species.name, "\n"))  #lê os arquivos
        tif.files <- list.files(paste0(models.dir, "/", species.name, "/present/",
                                       final.dir),
                                full.names = T, pattern = paste0(whi, ".*tif$"))

        if (length(tif.files) == 0) {
            cat(paste("No models to ensemble from for", species.name, "\n"))
        } else {
            cat(paste(length(tif.files),
                      "models to ensemble from for", species.name, "\n"))
            mod2 <- raster::stack(tif.files)
            #if (length(tif.files) == 1) {
             #   ensemble.mean <- mod2
            #} else {
                ensemble.mean <- raster::overlay(mod2, fun = function(x) {
                    return(mean(x, na.rm = T))
                    }
                    )
                ensemble.sd <- raster::overlay(mod2, fun = function(x) {
                    return(sd(x, na.rm = T))
                    }
                    )
                ensemble.min <- raster::overlay(mod2, fun = function(x) {
                    return(min(x, na.rm = T))
                    }
                    )
                ensemble.max <- raster::overlay(mod2, fun = function(x) {
                    return(max(x, na.rm = T))
                    }
                    )
                ensemble.median <- raster::overlay(mod2, fun = function(x) {
                    return(stats::median(x, na.rm = T))
                    }
                    )
            ensemble.mods <- raster::stack(ensemble.mean, ensemble.median, ensemble.sd,
                                   ensemble.min, ensemble.max)
            names(ensemble.mods) <- c("mean", "median", "sd", "min", "max")

            coord <- occs[occs$sp == species.name, c("lon", "lat")]

            if (write_png) {
                png(filename = paste0(models.dir, "/", species.name, "/present/",
                                      ensemble.dir, "/", species.name, "_", whi,
                                      "_ensemble_mean.png"),
                    res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
                raster::plot(ensemble.mean)
                maps::map("world",
                          c("", "South America"),
                          add = T,
                          col = "grey")
                points(coord, pch = 21, cex = 0.6,
                       bg = scales::alpha("cyan", 0.6))
                dev.off()
            }
            if (write_raw_map) {
                png(filename = paste0(models.dir, "/", species.name, "/present/",
                                  ensemble.dir, "/", species.name, "_", whi,
                                  "_ensemble_without_margins.png"),
                    bg = "transparent",
                res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
                raster::image(ensemble.mean, col = rev(terrain.colors(25)),
                          axes = F, asp = 1)
                dev.off()
                }

            raster::writeRaster(ensemble.mods,
                                filename = paste0(models.dir, "/", species.name,
                                                  "/present/",
                                                  ensemble.dir, "/", species.name, "_",
                                                  whi,
                                                  "_ensemble.tif"),
                                bylayer = T,
                                suffix = "names",
                                overwrite = T)

            #### Consensus models
            if (consensus == TRUE) {
                ensemble.consensus <- ensemble.mean >= consensus.level
                raster::writeRaster(ensemble.consensus,
                                    filename = paste0(models.dir, "/", species.name,
                                                      "/present/",
                                                      ensemble.dir, "/", species.name,
                                                      "_", whi,
                                                      "_ensemble", "_meanconsensus",
                                                      consensus.level * 100,
                                                      ".tif"), overwrite = T)


                if (write_png) {
                png(filename = paste0(models.dir, "/", species.name, "/present/",
                                      ensemble.dir, "/",
                                      species.name, "_", whi,
                                      "_ensemble", "_meanconsensus",
                                      consensus.level * 100, ".png"), res = 300,
                    width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
                raster::plot(ensemble.consensus)
                maps::map("world", c("", "South America"),
                          add = T, col = "grey")
                points(coord, pch = 19, cex = 0.3,
                       col = scales::alpha("cyan", 0.6))
                dev.off()
                }
                if (write_raw_map) {
                png(filename = paste0(models.dir, "/", species.name, "/present/",
                                      ensemble.dir, "/",
                                       species.name, "_", whi, "_ensemble",
                                       consensus.level * 100,
                                      "without_margins.png"),
                    bg = "transparent",
                    res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
                raster::image(ensemble.consensus,
                              col = rev(terrain.colors(25)), axes = F, asp = 1)
                dev.off()
                }
            }
        }
    }
    print(date())
}


