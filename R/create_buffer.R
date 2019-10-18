#' Samples pseudoabsences inside a geographic or environmental buffer
#'
#'This function is used internally by function
#' \code{\link{setup_sdmdata}} to define the area where pseudoabsences will be
#' sampled in different ways. First, it can create a maximum inclusion buffer,
#' within which the pseudoabsences will be sampled, to restrict model evaluation
#'  to accesible areas. This can be performed by either setting a user-defined
#'  shapefile, or by selecting one of several distance measures, including a
#'  numeric fixed value, to draw a buffer around the occurrences. In addition to
#'  this, an euclidean environmental distance can be superimposed to the
#'  previous step. The third step aims to control for overfitting by excluding
#'  areas that are too close to the occurrence points, with parameter "\code{dist_min}".
#'  The function will return the resulting buffer as a RasterStack object with
#'  the same resolution and NA values of the predictors RasterStack.
#'#'
#' @inheritParams setup_sdmdata
#' @inheritParams euclidean
#' @param buffer_type Character string indicating whether the buffer should be
#' calculated using the "\code{mean}", "\code{median}", "\code{maximum}"
#' distance between occurrence points, or an absolute geographic
#' "\code{distance}". If set to "\code{user}", a user-supplied shapefile will be
#' used as a sampling area, and \code{buffer_shape} needs to be specified. If
#' NULL, no distance buffer is applied. If set to "\code{distance}",
#' \code{dist_buf} needs to be specified
#' @param dist_buf Defines the width of the buffer. Needs to be specified if
#' \code{buffer_type = "distance"}. Distance unit is in the same unit of the
#' RasterStack of predictor variables
#' @param env_buffer Logical. Should an euclidean environmental filter be
#' applied? If TRUE, \code{env_distance} and \code{max_env_dist} need to be
#' specified
#' @param env_distance Character. Type of environmental distance, any in
#' "\code{centroid}", "\code{mindist}". Defaults to "\code{centroid}", the
#' distance of each raster pixel to the environmental centroid of the
#' distribution. When set to "\code{mindist}", the minimum distance of each
#' raster pixel to any of the occurrence points is calculated. Needs to be
#' specified if \code{env_buffer = TRUE}. A maximum value needs to be specified
#'  (parameter \code{max_env_dist})
#' @param max_env_dist Numeric. Since large negative values can arise
#'  during the calculation of the euclidean environmental distance, this
#'  parameter sets a maximum value to cut the environmental distance buffer.
#'  Expressed in quantiles, from 0: all values to 1: no values. Defaults to 0.5,
#'  the median value. Needs to be specified if \code{env_buffer = TRUE}
#' @param dist_min Optional, numeric. A distance for the exclusion of areas too
#' close from the occurrence points. Distance unit is in the same unit of the
#' RasterStack of predictor variables
#' @param buffer_shape User-defined buffer shapefile in which pseudoabsences
#' will be generated. Needs to be specified if \code{buffer_type = "user"}
#' @param write_buffer Logical. Should the resulting RasterStack be written?
#' Defaults to FALSE
#' @return Table of pseudoabsence points sampled within the selected distance
#' @return A buffer around the occurrence points
#'#' @references VanDerWal, J., Shoo, L. P., Graham, C., & Williams, S. E. (2009).
#' Selecting pseudo-absence data for presence-only distribution modeling: How
#' far should you stray from what you know? Ecological Modelling, 220(4),
#' 589-594. doi:10.1016/j.ecolmodel.2008.11.010
#' @seealso \code{\link[raster]{buffer}}
#' @seealso \code{\link[dismo]{randomPoints}}
#' @examples
#' library(raster)
#' sp <- names(example_occs)[1]
#' occs <- example_occs[[1]]
#' buf <- create_buffer(species_name = sp,
#'                      occurrences = occs,
#'                      predictors = example_vars)
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
                          max_env_dist = 0.5,
                          buffer_shape,
                          models_dir = "./models",
                          write_buffer = FALSE) {
    sp::coordinates(occurrences) <- ~lon + lat
    raster::crs(occurrences) <- raster::crs(predictors)
    if (is.null(buffer_type) |
        !buffer_type %in% c("distance", "mean", "median", "maximum", "user")) {
        message("buffer_type NULL or not recognized, returning predictors")
        r_buffer <- predictors[[1]]
    }
    if (buffer_type %in% c("distance", "mean", "median", "maximum", "user")) {
            if (buffer_type == "user") {
                if (missing(buffer_shape) | !class(buffer_shape) %in%
                    c("SpatialPolygonsDataFrame", "SpatialPolygonsDataFrame")) {
                    stop("if buffer_type == 'user', buffer_shape needs to be
                         specified and to be a shapefile")
                }
                if (class(buffer_shape) %in%
                    c("SpatialPolygonsDataFrame", "SpatialPolygonsDataFrame")) {
                    buffer.shape <- buffer_shape
                }
            }
        if (buffer_type %in% c("distance", "mean", "median", "maximum")) {
            if (buffer_type %in% c("distance")) {
                if (missing(dist_buf)) stop("dist_buf must be set when using a
                                            distance buffer")
                else dist.buf <- dist_buf
            }
        if (buffer_type %in% c("mean", "median", "maximum")) {
            dists <- rgeos::gDistance(spgeom1 = occurrences, byid = TRUE)
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
                                       byid = FALSE, width = dist.buf)
        }
    # crops the predictors to that shape to rasterize
    r_buffer <- raster::crop(predictors[[1]], buffer.shape)
    # masks the buffer to avoid sampling outside the predictors
    r_buffer <- raster::mask(r_buffer, buffer.shape)
    }
    if (env_buffer == TRUE) {
        if (missing(env_distance))
            stop(paste("The type of environmental distance ('centroid',
                       'mindist') must be specified"))
        if (missing(max_env_dist))
            stop(paste("A quantile for maximum environmental distance must be
                       specified"))
        message("Applying environmental filter")

        env.buffer <- euclidean(predictors = predictors,
                                occurrences = occurrences,
                                env_dist = env_distance)
        q <- quantile(raster::getValues(env.buffer),
                      max_env_dist, names = FALSE, na.rm = TRUE)
        env.buffer[env.buffer <= q] <- NA
        # we create a shapefile so it can be masked like the other types
        env.shape <- raster::rasterToPolygons(env.buffer, dissolve = TRUE)
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
                                           byid = FALSE,
                                           width = dist_min)

        r_buffer <- raster::mask(r_buffer, buffer.shape.min, inverse = TRUE)
    }

    if (write_buffer) {
        setup.folder <- paste0(models_dir, "/", species_name, "/present",
                               "/data_setup")
        if (file.exists(setup.folder) == FALSE)
            dir.create(setup.folder, recursive = TRUE)
        writeRaster(r_buffer, filename = paste0(setup.folder, "/buffer"),
                    format = "GTiff", overwrite = TRUE)
    }
    return(r_buffer)
}
