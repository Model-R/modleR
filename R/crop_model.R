#' Crop environmental raster based on speficif buffer
#'
#' @param modelo raster with the ecological niche model 
#' @param mask shapefile object to be used as mask for cropping the model
#' @author Andrea Sánchez-Tapia
#' @return Croppped environmental raster based on specific buffer
#' @seealso \code{\link[raster]{crop}}
#' @seealso \code{\link[dismo]{mask}}
#' @import raster
#' @export

crop_model <- function(modelo, mask) {
    modelo <- raster::crop(modelo, mask)
    modelo <- raster::mask(modelo, mask)
    return(modelo)
}
