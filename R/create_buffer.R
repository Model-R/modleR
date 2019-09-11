#' Samples pseudoabsences inside a geographic buffer.
#'
#' @inheritParams setup_sdmdata
#' @inheritParams euclidean
#' @param buffer_type Character string indicating whether the buffer should be
#' calculated using the "mean", "median", "maximum" distance between occurrence
#' points, or an absolute geographic "distance". If NULL pseudoabsences are randomly generated in the entire area
#' of the RasterStack of predictor variables.
#' filled with predictors. If set to "distance", "dist_buf" needs to
#' be specified. If set to "user", "buffer_shape" needs to be specified.
#' @param dist_buf Defines the width of the buffer. Needs to be specified if buffer_type = "distance".
#' Distance unit is in the same unit of the RasterStack of predictor variables
#' @param env_buffer Logical. Should an euclidean environmental filter be executed? If TRUE, \code{"env_distance"} and \code{"max_env_dist"} need to be specified.
#' @param env_distance Character. Type of environmental distance \code{"centroid", "mindist"}. The default is \code{"centroid"}.
#' @param max_env_dist Numeric. Quantile set as maximum environmental distance. Needs to be specified if env_buffer = TRUE
#' @param dist_min Optional, a distance for the exclusion buffer.
#' Distance unit is in the same unit of the RasterStack of predictor variables
#' @param buffer_shape User-defined buffer shapefile in which pseudoabsences will be generated.
#' Needs to be specified if buffer_type = "user"
#' @param write_buffer Logical. Should the resulting raster file be written? defaults to FALSE
#' @return Table of pseudoabsence points sampled within the selected distance
#' @author Felipe Barros
#' @author Andrea Sánchez-Tapia
#' @author Diogo S.B. Rocha
#' @author Sara Mortara
#' @return A buffer around the occurrence points
#' @details It will return a raster object with the same resolution and extent, and cropped by the predictors stack.
#' @references VanDerWal, J., Shoo, L. P., Graham, C., & Williams, S. E. (2009). Selecting pseudo-absence data for presence-only distribution modeling: How far should you stray from what you know? Ecological Modelling, 220(4), 589-594. doi:10.1016/j.ecolmodel.2008.11.010
#' @seealso \code{\link[raster]{buffer}}
#' @seealso \code{\link[dismo]{randomPoints}}
#' @examples
#' library(raster)
#' sp <- names(coordenadas)[1]
#' occs <- coordenadas[[1]]
#' buf <- create_buffer(species_name=sp, occurrences=occs, predictors=example_vars)
#' plot(buf)
#'
#' @import raster
#' @importFrom dismo randomPoints
#' @importFrom rgeos gBuffer
#' @importFrom stats median
#' @export
create_buffer <- function(species_name,
                          occurrences,
                          lon = "lon",
                          lat = "lat",
                          predictors,
                          buffer_type = "median",
                          dist_buf = NULL,
                          env_buffer = FALSE,
                          env_distance = "centroid",
                          dist_min = NULL,
                          max_env_dist = 0.5,#quantil
                          buffer_shape,
                          models_dir = "./models",
                          write_buffer = F) {
    sp::coordinates(occurrences) <- ~lon + lat
    raster::crs(occurrences) <- raster::crs(predictors)
    if (is.null(buffer_type) |
        !buffer_type %in% c("distance", "mean", "median", "maximum", "user")) {
        message("buffer_type NULL or not recognized, returning predictors")
        r_buffer <- predictors[[1]]
    }
    if (buffer_type %in% c("distance", "mean", "median", "maximum", "user")) {
            if (buffer_type == "user") {
                if (missing(buffer_shape) | !class(buffer_shape) %in% c("SpatialPolygonsDataFrame", "SpatialPolygonsDataFrame")) {
                    stop("if buffer_type == 'user', buffer_shape needs to be specified and to be a shapefile")
                }
                if (class(buffer_shape) %in% c("SpatialPolygonsDataFrame", "SpatialPolygonsDataFrame")) {
                    buffer.shape <- buffer_shape
                }
            }
        if (buffer_type %in% c("distance", "mean", "median", "maximum")) {
            if (buffer_type %in% c("distance")) {
                if (missing(dist_buf)) stop("dist_buf must be set when using a distance buffer")
                else dist.buf <- dist_buf
            }
        if (buffer_type %in% c("mean", "median", "maximum")) {
            dists <- rgeos::gDistance(spgeom1 = occurrences, byid = T)
            if (buffer_type == "mean") {
                dist.buf <- mean(dists)
            }
            if (buffer_type == "maximum") {
                dist.buf <-  max(dists)
            }
            if (buffer_type == "median") {
                dist.buf <- stats::median(dists)
            }
        }
        # creates the buffer - it's a shapefile
        buffer.shape <- rgeos::gBuffer(spgeom = occurrences,
                                       byid = F, width = dist.buf)
        }
    # crops the predictors to that shape to rasterize
    r_buffer <- raster::crop(predictors[[1]], buffer.shape)
    # masks the buffer to avoid sampling outside the predictors
    r_buffer <- raster::mask(r_buffer, buffer.shape)
    }
    if (env_buffer == TRUE) {
        if (missing(env_distance))
            stop(paste("The type of environmental distance ('centroid', 'mindist') must be specified"))
        if (missing(max_env_dist))
            stop(paste("A quantile for maximum environmental distance must be specified"))
        message("Applying environmental filter")

        env.buffer <- euclidean(predictors = predictors,
                                occurrences = occurrences,
                                algo = env_distance)
        q <- quantile(raster::getValues(env.buffer),
                      max_env_dist, names = F, na.rm = T)
        env.buffer[env.buffer <= q] <- NA #era diferente y el sample no funcionaba
        # we create a shapefile so it can be masked like the other types
        env.shape <- raster::rasterToPolygons(env.buffer, dissolve = T)
        #env.shape <- env.shape[env.shape@data$layer == 1,]#no need anymore
        r_buffer <- raster::crop(r_buffer, env.shape)
        r_buffer <- raster::mask(r_buffer, env.shape)
    }
    if (is.numeric(dist_min)) {
        if (exists("dist.buf")) {
            if (dist_min >= dist.buf) {
                warning("dist_min is higher than dist_buf")
            }
        }

        buffer.shape.min <- rgeos::gBuffer(spgeom = occurrences,
                                           byid = F,
                                           width = dist_min)

        r_buffer <- raster::mask(r_buffer, buffer.shape.min, inverse = TRUE)
    }

    if (write_buffer) {
        setup.folder <- paste0(models_dir, "/", species_name, "/present", "/data_setup")
        if (file.exists(setup.folder) == FALSE)
            dir.create(setup.folder, recursive = T)
        writeRaster(r_buffer, filename = paste0(setup.folder, "/buffer"),
                    format = "GTiff", overwrite = T)
    }
    return(r_buffer)
}

