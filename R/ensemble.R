#' Joins ENM from several algorithms, creating a model per species.
#'
#' @param sp A character string with the species name
#' @param models.dir Character string. Folder path where the input files are located
#' @param final.dir Character string, name of the subfolder where the files for the final models are located
#' @param ensemble.dir Character string, name of the folder to save the output files. A subfolder will be created.
#' @param occs A two-column data frame with the occurrence points
#' @param which.models How will the ensemble be built? Currently "Final.bin.mean3" and/or "Final.mean.bin7"
#' @param consensus Logical. Will a consensus rule be applied?
#' @param consensus.level Threshold for the consensus rule, betwen 0 and 1 (0.5 means a majority rule).
#' @param write_png Write png? Defaults to TRUE
#' @import raster
#' @export
ensemble <- function(sp,
                     models.dir = "./models",
                     final.dir = "final_models",
                     ensemble.dir = "ensemble",
                     occs = spp.filt,
                     which.models = c("Final.bin.mean3", "Final.mean.bin7"),
                     consensus = F,
                     consensus.level = 0.5,
                     write_png = T) {

    ## pasta de output
    if (file.exists(paste0(models.dir, "/", sp, "/present/", ensemble.dir, "/")) == FALSE) {
        dir.create(paste0(models.dir, "/", sp, "/present/", ensemble.dir, "/"))
    }

    ## para cada tipo de modelo
    for (whi in which.models) {
        cat(paste(whi, "-", sp, "\n"))  #lê os arquivos
        tif.files <- list.files(paste0(models.dir, "/", sp, "/present/", final.dir),
                                full.names = T, pattern = paste0(whi, ".*tif$"))

        if (length(tif.files) == 0) {
            cat(paste("No models to ensemble from for", sp, "\n"))
        } else {
            cat(paste(length(tif.files), "models to ensemble from for", sp, "\n"))
            mod2 <- raster::stack(tif.files)
            if (length(tif.files) == 1) {
                ensemble.m <- mod2
            } else {
                # plot(mod2)
                #ensemble.m <- raster::mean(mod2)
                ensemble.m <- raster::overlay(mod2, fun = function(x) { return(mean(x, na.rm = T)) })
                ensemble.sd <- raster::overlay(mod2, fun = function(x) { return(sd(x, na.rm = T)) })
            }
            coord <- occs[occs$sp == sp, c("lon", "lat")]

            if (write_png) {
            png(filename = paste0(models.dir, "/", sp, "/present/",
                                  ensemble.dir, "/",
                                  sp, "_", whi, "_ensemble.png"),
                 res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
            par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
            raster::plot(ensemble.m)
            maps::map("world", c("", "South America"), add = T, col = "grey")
            #points(coord, pch = 21, cex = 0.6, bg = scales::alpha("cyan", 0.6))
            dev.off()
            }

            png(filename = paste0(models.dir, "/", sp, "/present/",
                                  ensemble.dir, "/",
                                  sp, "_", whi, "_ensemble_without_margins.png"), bg = "transparent",
                res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
            par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
            raster::image(ensemble.m, col = rev(terrain.colors(25)), axes = F, asp = 1)
            dev.off()

            # o ensemble cru
            raster::writeRaster(ensemble.m,
                                filename = paste0(models.dir, "/", sp,
                                                  "/present/",
                                                  ensemble.dir, "/", sp, "_",
                                                  whi,
                                                  "_ensemble.tif"),
                                overwrite = T)

            #### Consensus models
            if (consensus == TRUE) {
                ensemble.consensus <- ensemble.m >= consensus.level
                raster::writeRaster(ensemble.consensus,
                                    filename = paste0(models.dir, "/", sp,
                                                      "/present/",
                                                      ensemble.dir, "/", sp,
                                                      "_", whi,
                                                      "_ensemble",
                                                      consensus.level * 100,
                                                      ".tif"), overwrite = T)

                if (write_png) {
                png(filename = paste0(models.dir, "/", sp, "/present/",
                                      ensemble.dir, "/",
                                      sp, "_", whi, "_ensemble",
                                      consensus.level * 100, ".png"), res = 300,
                    width = 410 * 300/72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
                raster::plot(ensemble.consensus)
                maps::map("world", c("", "South America"), add = T, col = "grey")
                #points(coord, pch = 19, cex = 0.3, col = scales::alpha("cyan", 0.6))
                dev.off()
                }

                png(filename = paste0(models.dir, "/", sp, "/present/", ensemble.dir, "/",
                                       sp, "_", whi, "_ensemble",
                                       consensus.level * 100, "without_margins.png"), bg = "transparent",
                    res = 300, width = 410 * 300/72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
                raster::image(ensemble.consensus, col = rev(terrain.colors(25)), axes = F, asp = 1)
                dev.off()
            }
        }
    }
    #return(ensemble.m)
    print(date())
}
