#' Samples pseudoabsences inside a geographic buffer
#'
#' @inheritParams setup_sdmdata
#' @param buffer_type Character string indicating whether the buffer should be
#' calculated using the "mean", "median", "maximum" distance between occurrence points, or an absolute "distance". If set to "distance",
#'         "dist_buf" needs to be specified. If set to "user", "buffer_shape" needs to be specified.
#' @param dist_buf Defines the width of the buffer. Needs to be specified if buffer_type = "distance"
#' @param buffer_shape User-defined buffer shapefile. Needs to be specified if buffer_type = "user"
#' @param predictors A RasterStack of predictor variables
#' @param write_buffer Logical. Should the resulting raster file be written? defaults to FALSE
#' @param ... Other parameters from writeRaster
#' @return Table of pseudoabsence points sampled within the selected distance
#' @author Felipe Barros
#' @author Andrea Sánchez-Tapia
#' @author Diogo S.B. Rocha#'
#' @author Sara Mortara
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
                          buffer_type = NULL,
                          predictors,
                          dist_buf = NULL,
                          dist_min = NULL
                          buffer_shape = NULL,
                          models_dir = "./models",
                          species_name,
                          write_buffer = F,
                          ...) {
    sp::coordinates(occurrences) <- ~lon + lat
    raster::crs(occurrences) <- raster::crs(predictors)
    if (is.null(buffer_type) | !buffer_type %in% c("distance", "mean", "median", "max", "user")) {
      warning("buffer_type NULL or not recognized, returning predictors")
      r_buffer <- predictors
      return(r_buffer)
  }
    

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
      
      if (is.numeric(dist_min)) {
      buffer.shape.min <- rgeos::gBuffer(spgeom = occurrences,
                                         byid = F, width = dist_min*0.00833333)
      #rasterizes to sample the random points
      r_buffer <- raster::crop(predictors, buffer.shape)
      # masks the buffer to avoid sampling outside the predictors
      r_buffer <- raster::mask(r_buffer, buffer.shape)
      r_buffer <- raster::mask(r_buffer, buffer.shape.min, inverse = TRUE)
      } else {
      r_buffer <- raster::crop(predictors, buffer.shape)
      # masks the buffer to avoid sampling outside the predictors
      r_buffer <- raster::mask(r_buffer, buffer.shape)
      }
      
      #creates the buffer - it's a shapefile
      buffer.shape <- rgeos::gBuffer(spgeom = occurrences, byid = F, width = dist.buf)
      }

    #rasterizes to sample the random points
    r_buffer <- raster::crop(predictors, buffer.shape)

    # masks the buffer to avoid sampling outside the predictors
    r_buffer <- raster::mask(r_buffer, buffer.shape)
    if (write_buffer) {
        partition.folder <- paste0(models_dir, "/", species_name, "/present", "/partitions")
        writeRaster(r_buffer, filename = paste0(partition.folder,"/buffer"), format = "GTiff", ...)
    }
    return(r_buffer)
}
