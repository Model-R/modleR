#' Samples pseudoabsences inside a geographic buffer
#'
#' @param coord  A data frame with occurrence data. It should contain only two columns:
#' lon and lat, in that order.
#' @param n.back Number of pseudoabsence points
#' @param buffer.type Character string indicating whether the buffer should be
#' calculated using the mean, median or maximum distance between occurrence points
#' @param seed for reproducibility purposes
#' @param predictors A RasterStack of predictor variables
#' @return Table of pseudoabsence points sampled within the selected distance
#' @author Felipe Barros
#' @author Fabrício Vilasboas
#' @author Andrea Sánchez-Tapia
#' @details The sampling is performed by dismo::randomPoints() excluding the presence points (exclupep =TRUE)
#' @references VanDerWal, J., Shoo, L. P., Graham, C., & Williams, S. E. (2009). Selecting pseudo-absence data for presence-only distribution modeling: How far should you stray from what you know? Ecological Modelling, 220(4), 589-594. doi:10.1016/j.ecolmodel.2008.11.010
#' @seealso \code{\link[raster]{buffer}}
#' @seealso \code{\link[dismo]{randomPoints}}
#' @examples
#' create_buffer(coordinates = coordenadas, n.back = 500, buffer.type = "mean",
#'  predictors = variaveis_preditoras)
#'
#' @import raster
#' @importFrom dismo randomPoints
#'
create_buffer <- function(coord,
                          n.back,
                          buffer.type,
                          seed = 512,
                          predictors) {

    sp::coordinates(coord) <- ~lon + lat
    raster::crs(coord) <- raster::crs(predictors)
    if (buffer.type == "mean")
        dist.buf <- mean(sp::spDists(x = coord,
                                     longlat = TRUE,
                                     segments = FALSE))
    if (buffer.type == "max")
        dist.buf <-  max(sp::spDists(x = coord,
                                    longlat = TRUE,
                                    segments = FALSE))
    if (buffer.type == "median")
        dist.buf <- stats::median(sp::spDists(x = coord,
                                              longlat = TRUE,
                                              segments = FALSE))

    #creates the buffer - it's a shapefile
    buffer.shape <- raster::buffer(coord,
                                   width = dist.buf * 1000,
                                   dissolve = TRUE)

    #rasterizes to sample the random points
    r_buffer <- raster::rasterize(buffer.shape,
                                  predictors,
                                  field = buffer.shape@plotOrder)
    # masks the buffer to avoid sampling outside the predictors
    r_buffer <- raster::mask(r_buffer, predictors[[1]])


    # Samples random points
    set.seed(seed + 2)
    backgr <- dismo::randomPoints(mask = r_buffer,
                                  n = n.back,
                                  p = coord,
                                  excludep = T)
    rm(buffer.shape)
    rm(r_buffer)
    gc()
    return(backgr)
}
