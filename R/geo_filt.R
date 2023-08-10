#' Excludes occurrence points within a given geographic distance
#'
#' @inheritParams setup_sdmdata
#' @param min_distance Numeric. Minimum distance between points in the unit of
#' the predictor raster (most likely degrees). This value will be used to create a raster
#' with a \code{min_distance} resolution, from which the redundant occurences
#' will be excluded. 
#' @return Data frame of the occurrence points thinned to have the minimum
#' distance indicated in \code{min_distance} between them
#' @author Diogo S. B. Rocha
#' @references Varela, S., Anderson, R. P., García-Valdés, R., &
#' Fernández-González, F. (2014). Environmental filters reduce the effects of
#' sampling bias and improve predictions of ecological niche models.
#' Ecography, 37(11), 1084-1091. doi:10.1111/j.1600-0587.2013.00441.x
#' @seealso \code{\link[dismo]{gridSample}}
#' @seealso \pkg{spThin}
#'
#' @importFrom dismo gridSample
geo_filt <- function(occurrences,
                     lon = "lon",
                     lat = "lat",
                     min_distance = 1) {
    
    res <- min_distance
    r <- raster::raster(extent(range(occurrences[, lon]),
                               range(occurrences[, lat])) + res)
    res(r) <- res
    pts <- dismo::gridSample(occurrences, r, n = 1)
    message(paste0(dim(pts)[1],
               " Points remaining after the geographic filter of ",
               min_distance, "degrees", "\n"))
    return(pts)
}
