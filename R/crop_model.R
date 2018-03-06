
crop_model <- function(modelo, mascara) {
    modelo <- raster::crop(modelo, mascara)
    modelo <- raster::mask(modelo, mascara)
    return(modelo)
}
