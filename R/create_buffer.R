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
                          predictors,
                          buffer_type = NULL,
                          dist_buf = NULL,
                          models_dir = "./models",
                          species_name,
                          write_buffer = T) {
    sp::coordinates(occurrences) <- ~lon + lat
    raster::crs(occurrences) <- raster::crs(predictors)
    if (is.null(buffer_type) | !buffer_type %in% c("distance", "mean", "median", "max", "user")) {
      warning("buffer_type NULL or not recognized, returning predictors")
      r_buffer <- predictors
      return(r_buffer)
    }
    if (buffer_type == "user") {
        if (is.null(buffer_shape) | class(buffer_shape) != "SpatialPolygonsDataFrame") {
            stop("if buffer_type == 'user', buffer_shape needs to be specified and to be a shapefile")
        }
        if (class(buffer_shape) == "SpatialPolygonsDataFrame") {
            buffer.shape <- buffer_shape
        }
    }

    if (buffer_type %in% c("distance", "mean", "median", "max")) {
        if (buffer_type %in% c("distance")) {
            if (is.null(dist_buf)) stop("dist_buf must be set when using a distance buffer")
            else dist.buf <- dist_buf
            }
        if (buffer_type %in% c("mean", "median", "max")) {
            dists <- rgeos::gDistance(spgeom1 = occurrences, byid = T)
            if (buffer_type == "mean")
                dist.buf <- mean(dists)
            if (buffer_type == "max")
                dist.buf <-  max(dists)
            if (buffer_type == "median")
                dist.buf <- stats::median(dists)
            }
        #creates the buffer - it's a shapefile
        buffer.shape <- rgeos::gBuffer(spgeom = occurrences, byid = F, width = dist.buf)
        }
    partition.folder <- paste0(models_dir, "/", species_name, "/present", "/partitions")

    #rasterizes to sample the random points
    r_buffer <- raster::crop(predictors, buffer.shape)

    # masks the buffer to avoid sampling outside the predictors
    r_buffer <- raster::mask(r_buffer, buffer.shape)
    if (write_buffer) {
        writeRaster(r_buffer, filename = paste0(partition.folder,"/buffer"), format = "GTiff", overwrite = T)
    }
    return(r_buffer)
}
