write_raw_map <- function(x) {
    mod <- raster(x)
    png("ensemble_without_margins.png", bg = "transparent",
        res = 300, width = 410 * 300 / 72, height = 480 * 300 / 72)
    par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
    raster::image(mod, col = rev(terrain.colors(25)),
                  axes = F)
    dev.off()
    }
