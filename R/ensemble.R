#' Create ensemble models from several algorithms.
#'
#' \code{ensemble} takes the input of \code{final_model} and builds an ensemble
#'  model by calculating the mean of each final model per algorithm.
#'
#' @param species_name A character string with the species name
#' @param occurrences A two-column data frame with the occurrence points
#' @param models_dir Character string. Folder path where the input files
#'                   are located
#' @param final_dir Character string, name of the subfolder where the files
#'                  for the final models are located
#' @param ensemble_dir Character string, name of the folder to save the output
#'                     files. A subfolder will be created.
#' @param which_models Which final_model() will be used? Currently it can be:
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
#' @param consensus_level Threshold for the consensus rule, betwen 0 and 1
#'                        (0.5 means a majority rule).
#' @param write_png Write png? Defaults to TRUE
#' @param scale_models Sets the maximum value of the input models to 1.
#' Defaults to TRUE.
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
#' (optional) written in the \code{ensemble_dir} subfolder
ensemble_model <- function(species_name,
                           occurrences,
                           models_dir = "./models",
                           final_dir = "final_models",
                           ensemble_dir = "ensemble",
                           which_models = c("raw_mean"),
                           consensus = FALSE,
                           consensus_level = 0.5,
                           write_png = T,
                           scale_models = TRUE,
                           write_raw_map = F) {


    ## output folder
    if (file.exists(
        paste0(models_dir, "/", species_name, "/present/", ensemble_dir, "/")) == FALSE) {
        dir.create(paste0(models_dir, "/", species_name, "/present/", ensemble_dir, "/"))
    }

    ## for each model specified in final_models
    for (whi in which_models) {
        cat(paste(whi, "-", species_name, "\n"))  #lê os arquivos
        tif.files <- list.files(paste0(models_dir, "/", species_name, "/present/",
                                       final_dir),
                                full.names = T, pattern = paste0(whi, ".*tif$"))

        if (length(tif.files) == 0) {
            cat(paste("No", whi, "models to ensemble from for", species_name, "\n"))
        } else {
            cat(paste(length(tif.files), whi,
                      "models to ensemble from for", species_name, "\n"))
            mod2 <- raster::stack(tif.files)
            #scale models to 0-1
            if (scale_models == T) {
                mod2 <- mod2/raster::maxValue(mod2)
            }
            ensemble.mean <- raster::overlay(mod2, fun = function(x) {
                return(mean(x, na.rm = T))
                }
                )
            ensemble.sd <- raster::overlay(mod2, fun = function(x) {
                return(sd(x, na.rm = T))
                }
                )
            #ensemble.min <- raster::overlay(mod2, fun = function(x) {
             #   return(min(x, na.rm = T))
              #  }
               # )
            #ensemble.max <- raster::overlay(mod2, fun = function(x) {
             #   return(max(x, na.rm = T))
              #  }
               # )
            ensemble.median <- raster::overlay(mod2, fun = function(x) {
                return(stats::median(x, na.rm = T))
                }
                )
            ensemble.mods <- raster::stack(ensemble.mean,
                                           ensemble.median,
                                           ensemble.sd#,
                                           #ensemble.min,
                                           #ensemble.max
                                           )
            names(ensemble.mods) <- c("mean",
                                      "median",
                                      "sd"#,
                                      #"min", "max"
                                      )

            coord <- occurrences[occurrences$sp == species_name, c("lon", "lat")]

            if (write_png) {
                png(filename = paste0(models_dir, "/", species_name, "/present/",
                                      ensemble_dir, "/", species_name, "_", whi,
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
                png(filename = paste0(models_dir, "/", species_name, "/present/",
                                      ensemble_dir, "/", species_name, "_", whi,
                                      "_ensemble_median.png"),
                    res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
                raster::plot(ensemble.median)
                maps::map("world",
                          c("", "South America"),
                          add = T,
                          col = "grey")
                points(coord, pch = 21, cex = 0.6,
                       bg = scales::alpha("cyan", 0.6))
                dev.off()
                png(filename = paste0(models_dir, "/", species_name, "/present/",
                                      ensemble_dir, "/", species_name, "_", whi,
                                      "_ensemble_sd.png"),
                    res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
                raster::plot(ensemble.sd)
                maps::map("world",
                          c("", "South America"),
                          add = T,
                          col = "grey")
                points(coord, pch = 21, cex = 0.6,
                       bg = scales::alpha("cyan", 0.6))
                dev.off()
            }
            if (write_raw_map) {
                png(filename = paste0(models_dir, "/", species_name, "/present/",
                                  ensemble_dir, "/", species_name, "_", whi,
                                  "_ensemble_without_margins.png"),
                    bg = "transparent",
                res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
                raster::image(ensemble.mean, col = rev(terrain.colors(25)),
                          axes = F, asp = 1)

                dev.off()
                }

            raster::writeRaster(ensemble.mods,
                                filename = paste0(models_dir, "/", species_name,
                                                  "/present/",
                                                  ensemble_dir, "/", species_name, "_",
                                                  whi,
                                                  "_ensemble.tif"),
                                bylayer = T,
                                suffix = "names",
                                overwrite = T)

            #### Consensus models
            if (consensus == TRUE) {
                ensemble.consensus <- ensemble.mean >= consensus_level
                raster::writeRaster(ensemble.consensus,
                                    filename = paste0(models_dir, "/", species_name,
                                                      "/present/",
                                                      ensemble_dir, "/", species_name,
                                                      "_", whi,
                                                      "_ensemble", "_meanconsensus",
                                                      consensus_level * 100,
                                                      ".tif"), overwrite = T)


                if (write_png) {
                png(filename = paste0(models_dir, "/", species_name, "/present/",
                                      ensemble_dir, "/",
                                      species_name, "_", whi,
                                      "_ensemble", "_meanconsensus",
                                      consensus_level * 100, ".png"), res = 300,
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
                png(filename = paste0(models_dir, "/", species_name, "/present/",
                                      ensemble_dir, "/",
                                       species_name, "_", whi, "_ensemble",
                                       consensus_level * 100,
                                      "without_margins.png"),
                    bg = "transparent",
                    res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
                raster::image(ensemble.consensus,
                              col = rev(terrain.colors(25)), axes = F)
                dev.off()
                }
            }
        }
    }
    print(date())
}


