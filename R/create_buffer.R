#' Samples pseudoabsences inside a geographic buffer
#'
#' @param occurrences  A data frame with occurrence data. It should contain only two columns:
#' lon and lat, in that order.
#' @param buffer_type Character string indicating whether the buffer should be
#' calculated using the mean, median or maximum distance between occurrence points
#' @param dist_buf Optional, a distance in km for tbe buffer. If set it will override buffer_type
#' @param predictors A RasterStack of predictor variables
#' @return Table of pseudoabsence points sampled within the selected distance
#' @author Felipe Barros
#' @author Andrea Sánchez-Tapia
#' @author Fabrício Vilasboas
#' @return A buffer around the occurrence points
#' @details The sampling is performed by dismo::randomPoints() excluding the presence points (exclupep =TRUE)
#' @references VanDerWal, J., Shoo, L. P., Graham, C., & Williams, S. E. (2009). Selecting pseudo-absence data for presence-only distribution modeling: How far should you stray from what you know? Ecological Modelling, 220(4), 589-594. doi:10.1016/j.ecolmodel.2008.11.010
#' @seealso \code{\link[raster]{buffer}}
#' @seealso \code{\link[dismo]{randomPoints}}
#' @examples
#' library(raster)
#' library(dplyr)
#' species <- sort(unique(coordenadas$sp))
#' occs <- coordenadas %>% filter(sp == species[1]) %>% dplyr::select(lon, lat)
#' buf <- create_buffer(occs, "mean", example_vars)
#' plot(buf)
#'
#' @import raster
#' @importFrom dismo randomPoints
#' @importFrom rgeos gBuffer
#' @export
create_buffer <- function(occurrences,
                          buffer_type,
                          predictors,
                          dist_buf = NULL) {

    sp::coordinates(occurrences) <- ~lon + lat
    raster::crs(occurrences) <- raster::crs(predictors)
    if (buffer_type == "mean")
        dist.buf <- mean(sp::spDists(x = occurrences,
                                     longlat = FALSE,
                                     segments = FALSE))
    if (buffer_type == "max")
        dist.buf <-  max(sp::spDists(x = occurrences,
                                    longlat = FALSE,
                                    segments = FALSE))
    if (buffer_type == "median")
        dist.buf <- stats::median(sp::spDists(x = occurrences,
                                              longlat = FALSE,
                                              segments = FALSE))
    if (buffer_type == "distance")
        dist.buf <- dist_buf

    #creates the buffer - it's a shapefile
    buffer.shape <- rgeos::gBuffer(spgeom = occurrences,
                                   byid = F, width = dist.buf,
                                   )

    #rasterizes to sample the random points
    r_buffer <- raster::rasterize(buffer.shape,
                                  predictors,
                                  field = buffer.shape@plotOrder)
    # masks the buffer to avoid sampling outside the predictors
    r_buffer <- raster::mask(r_buffer, predictors[[1]])

    return(r_buffer)
}
