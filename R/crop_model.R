#' Crops and masks environmental raster based on specific mask
#'
#' @param modelo RasterStack object of the ecological niche model
#' @param mask A SpatialPolygonsDataFrame to be used to mask the final models
#' @return Croppped environmental raster based on specific buffer
#' @seealso \code{\link[raster]{crop}}
#' @seealso \code{\link[raster]{mask}}
#' @import raster
#' @export

crop_model <- function(modelo,
                       mask) {
    modelo <- raster::crop(modelo, mask)
    modelo <- raster::mask(modelo, mask)
    return(modelo)
}
