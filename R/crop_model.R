
crop_model <- function(modelo, mask) {
    modelo <- raster::crop(modelo, mask)
    modelo <- raster::mask(modelo, mask)
    return(modelo)
}
