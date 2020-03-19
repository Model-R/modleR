rescale_layer <- function(layers) {
    message("Standardizing models from 0 to 1")
    if (missing(layers)) {
        stop("No layers were provided. Please enter a Raster or a RasterStack")
        }
    for (i in 1:dim(layers)[3]) {
        stand <- function(x) {
            (x - min(layers[[i]][], na.rm = TRUE)) /
                (max(layers[[i]][], na.rm = TRUE) - min(layers[[i]][],
                                                     na.rm = TRUE))
      }
        bb <- raster::calc(layers[[i]], stand)
        bb
    if (i == 1) {
      cc <- raster::stack(bb)
      names(cc)[i] <- names(layers)[i]
    }
    else{
      cc <- raster::stack(cc, bb)
      names(cc)[i] <- names(layers)[i]
    }
    if (i == dim(layers)[3]) {
      layers <- cc
      rm(cc)
      return(layers)
    }
  }
}
