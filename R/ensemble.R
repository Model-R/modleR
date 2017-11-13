#' Gera um modelo só por espécie
#'
#' @param sp Um nome de espécie
#' @param models.dir Path do diretório onde estão os modelos
#' @param final.dir Path para subdiretório escrito com finalModel()
#' @param ensemble.dir Path para os arquivos de saída
#' @param occs Pontos de ocorrência da espécie
#' @param which.models Qual tipo de modelo final será juntado (ex. bin.mean3)
#' @param consensus Se será usada uma regra de consenso para cortar o ensemble
#' @param consensus.level Quanto dos modelos será retido (0,5 = maioria)
#' @return NULL
#' @import raster
#' @export
ensemble <- function(sp,
                     models.dir = "./models",
                     final.dir = "presfinal",
                     ensemble.dir = "ensemble",
                     occs = spp.filt,
                     which.models = c("Final.bin.mean3", "Final.mean.bin7"),
                     consensus = F,
                     consensus.level = 0.5) {

    ## pasta de output
    if (file.exists(paste0(models.dir, "/", sp, "/", ensemble.dir, "/")) == FALSE) {
        dir.create(paste0(models.dir, "/", sp, "/", ensemble.dir, "/"))
    }

    ## para cada tipo de modelo
    for (whi in which.models) {
        cat(paste(whi, "-", sp, "\n"))  #lê os arquivos
        tif.files <- list.files(paste0(models.dir, "/", sp, "/", final.dir),
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

            png(filename = paste0(models.dir, "/", sp, "/", ensemble.dir, "/",
                                   sp, "_", whi, "_ensemble.png"),
                 res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
            raster::plot(ensemble.m)
            maps::map("world", c("", "South America"), add = T, col = "grey")
            #points(coord, pch = 21, cex = 0.6, bg = scales::alpha("cyan", 0.6))
            dev.off()
            
            png(filename = paste0(models.dir, "/", sp, "/", ensemble.dir, "/",
                                  sp, "_", whi, "_ensemble_without_margins.png"),
                res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
            par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
            raster::image(ensemble.m, legend = F, axes = FALSE,box = F, col=rev(terrain.colors(25)))
            #points(coord, pch = 21, cex = 0.6, bg = scales::alpha("cyan", 0.6))
            dev.off()

            # o ensemble cru
            raster::writeRaster(ensemble.m,
                                filename = paste0(models.dir, "/", sp, "/",
                                                  ensemble.dir, "/", sp, "_", whi,
                                                  "_ensemble.tif"), overwrite = T)

            #### Consensus models
            if (consensus == TRUE) {
                ensemble.consensus <- ensemble.m >= consensus.level
                raster::writeRaster(ensemble.consensus,
                                    filename = paste0(models.dir, "/", sp, "/",
                                                      ensemble.dir, "/", sp, "_", whi,
                                                      "_ensemble", consensus.level * 100,
                                                      ".tif"), overwrite = T)
                
                png(filename = paste0(models.dir, "/", sp, "/", ensemble.dir, "/",
                                      sp, "_", whi, "_ensemble",
                                      consensus.level * 100, ".png"), res = 300,
                    width = 410 * 300/72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
                raster::plot(ensemble.consensus)
                maps::map("world", c("", "South America"), add = T, col = "grey")
                #points(coord, pch = 19, cex = 0.3, col = scales::alpha("cyan", 0.6))
                dev.off()


                png(filename = paste0(models.dir, "/", sp, "/", ensemble.dir, "/",
                                       sp, "_", whi, "_ensemble",
                                       consensus.level * 100, "without_margins", ".png"), res = 300,
                     width = 410 * 300/72, height = 480 * 300 / 72)
                par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
                raster::image(ensemble.consensus, legend = F, axes = FALSE, box = F, col=rev(terrain.colors(25)))
                raster::plot(ensemble.consensus, legend = F, axes = FALSE, box = F)
                #points(coord, pch = 19, cex = 0.3, col = scales::alpha("cyan", 0.6))
                dev.off()
            }
        }
    }
    #return(ensemble.m)
}
