#' Samples pseudoabsences inside a geographic or environmental buffer
#'
#' This function is used internally by function \code{\link{setup_sdmdata}} to
#' define the area where pseudoabsences will be sampled in different ways.
#' First, it can create an inclusion buffer, within which the pseudoabsences
#' will be sampled, to restrict model evaluation to accesible areas. This can be
#' performed by either setting a user-defined shapefile, or by drawing a buffer
#' around the species occurrences: a geographic distance fixed value, or the
#' mean, median or maximum pairwise distance between occurrences. In addition to
#' this, an euclidean environmental distance filter can be superimposed to the
#' previous step, to control for overfitting by excluding areas that are too
#' close to the occurrence points, either in the environmental space or in the
#' geographic space.
#' The function will return the resulting buffer as a RasterStack object with
#' the same resolution and NA values of the predictors RasterStack.
#'
#' @inheritParams setup_sdmdata
#' @param buffer_type Character string indicating whether the buffer should be
#' calculated using the "\code{mean}", "\code{median}", "\code{maximum}"
#' distance between occurrence points, or an absolute geographic
#' "\code{distance}". If set to "\code{user}", a user-supplied shapefile will be
#' used as a sampling area, and \code{buffer_shape} needs to be specified. If
#' NULL, no distance buffer is applied. If set to "\code{distance}",
#' \code{dist_buf} needs to be specified
#' @param buffer_shape User-defined buffer shapefile in which pseudoabsences
#' will be generated. Needs to be specified if \code{buffer_type = "user"}
#' @param dist_buf Defines the width of the buffer. Needs to be specified if
#' \code{buffer_type = "distance"}. Distance unit is in the same unit of the
#' RasterStack of predictor variables
#' @param env_filter Logical. Should an euclidean environmental filter be
#' applied? If TRUE, \code{env_distance} and
#' \code{min_env_dist} need to be specified. Areas closest than
#' \code{min_env_dist} (expressed in quantiles in the environmental space)will
#' be omitted from the pseudoabsence sampling
#' @param env_distance Character. Type of environmental distance, either
#' "\code{centroid}" or "\code{mindist}". Defaults to "\code{centroid}", the
#' distance of each raster pixel to the environmental centroid of the
#' distribution. When set to "\code{mindist}", the minimum distance of each
#' raster pixel to any of the occurrence points is calculated. Needs to be
#' specified if \code{env_filter = TRUE}. A minimum value needs to be
#' specified (parameter \code{min_env_dist})
#' @param min_env_dist Numeric. Sets a minimum value to exclude the areas
#' closest (in the environmental space) to the occurrences or their centroid,
#' expressed in quantiles, from 0 (the closest) to 1. Defaults to 0.05,
#' excluding areas belonging to the 5% closest environmental values. Note that
#' since this is based on quantiles, and environmental similarity can take large
#' negative values, this is an abitrary value
#' @param min_geog_dist Optional, numeric. A distance for the exclusion of the
#' areas closest to the occurrence points (in the geographical space). Distance
#'  unit is in the same unit of the RasterStack of predictor variables
#' @param write_buffer Logical. Should the resulting buffer RasterLayer be
#' written? Defaults to FALSE
#' @return Table of pseudoabsence points sampled within the selected distance
#' @return A RasterLayer object containing the final buffer around the
#' occurrence points
#' @references VanDerWal, J., Shoo, L. P., Graham, C., & Williams, S. E. (2009).
#' Selecting pseudo-absence data for presence-only distribution modeling: How
#' far should you stray from what you know? Ecological Modelling, 220(4),
#' 589-594. doi:10.1016/j.ecolmodel.2008.11.010
#' @seealso \code{\link[raster]{buffer}}  in \pkg{raster} package
#' @seealso \code{\link[dismo]{randomPoints}}  in \pkg{dismo} package
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
                          buffer_type = "none",
                          buffer_shape,
                          dist_buf = NULL,
                          env_filter = FALSE,
                          env_distance = "centroid",
                          min_env_dist = NULL,
                          min_geog_dist = NULL,
                          models_dir = "./models",
                          write_buffer = FALSE) {
    sp::coordinates(occurrences) <- ~lon + lat
    raster::crs(occurrences) <- raster::crs(predictors)
    if (!buffer_type %in% c("distance", "mean", "median", "maximum", "user")) {
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
        r_buffer <- raster::crop(predictors[[1]],  buffer.shape)
        # masks the buffer to avoid sampling outside the predictors
        r_buffer <- raster::mask(r_buffer, buffer.shape)
    }
    if (env_filter == TRUE) {
        if (missing(env_distance))
            stop(paste("The type of environmental distance ('centroid',
                       'mindist') must be specified"))
        if (missing(min_env_dist))
            stop(paste("A quantile for mininum environmental distance
                       must be specified"))
        message("Applying environmental filter")

        env.filter <- euclidean(predictors = predictors,
                                occurrences = occurrences,
                                env_dist = env_distance)
        if (!missing(min_env_dist)) {
        #min_env_dist
        q_min <- quantile(raster::getValues(env.filter),
                          (1 - min_env_dist), names = FALSE, na.rm = TRUE)
        env.filter[env.filter >= q_min] <- NA
        }
        # we create a shapefile so it can be used to mask like the other types
        env.shape <- raster::rasterToPolygons(env.filter, dissolve = TRUE)
        r_buffer <- raster::crop(r_buffer, env.shape)
        r_buffer <- raster::mask(r_buffer, env.shape)
    }
    if (is.numeric(min_geog_dist)) {
        if (exists("dist.buf")) {
            if (min_geog_dist >= dist.buf) {
                warning("min_geog_dist is higher than dist_buf")
            }
        }

        buffer.shape.min <- rgeos::gBuffer(spgeom = occurrences,
                                           byid = FALSE,
                                           width = min_geog_dist)

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
