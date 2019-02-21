#' Rescaling rasters to values between 0 and 1.
#'
#' This function rescales rasters to values between 0 and 1.
#' @param layers A RasterStack or Rasterlayer of predictor variables to scale
#' @return A RasterStack or Rasterlayer with values between 0 and 1.
#' @author Diogo S. B. Rocha
#' @examples
#' rescale.layer(example_vars)
#' @import raster
#' @export
rescale_layer <- function(layers) {
  cat(paste("Standardizing models from 0 to 1", "\n"))
  if (missing(layers)) {
    stop("No layers were provided. Please enter a Raster layer or a Rasterstack")
  }
  #pb <- txtProgressBar(min = 1, max = dim(layers)[3], style = 3)
  for (i in 1:dim(layers)[3]) {
    # setTxtProgressBar(pb, i)
    stand <-
      function(x) {
        (x - min(layers[[i]][], na.rm = T)) / (max(layers[[i]][], na.rm = T) - min(layers[[i]][], na.rm = T))
      }
    bb <- calc(layers[[i]], stand)
    bb
    if (i == 1) {
      cc <- stack(bb)
      names(cc)[i] = names(layers)[i]
    }
    else{
      cc <- stack(cc, bb)
      names(cc)[i] = names(layers)[i]
    }
    if (i == dim(layers)[3]) {
      layers <- cc
      rm(cc)
      return(layers)
    }
  }
}
