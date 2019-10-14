#' Excludes occurrence points within a geographic distance
#'
#' @inheritParams setup_sdmdata
#' @param min_distance numeric. Minimum distance between points in the unit of the predictor raster.
#' @return Table of occurrence points with minimum distance indicated in min.distance
#' @references Varela, S., Anderson, R. P., García-Valdés, R., &
#' Fernández-González, F. (2014). Environmental filters reduce the effects of
#' sampling bias and improve predictions of ecological niche models.
#' Ecography, 37(11), 1084-1091. doi:10.1111/j.1600-0587.2013.00441.x
#' @seealso \code{\link[dismo]{gridSample}}
#'
#' @importFrom dismo gridSample
#' @export
geo_filt <- function(occurrences,
                     lon = "lon",
                     lat = "lat",
                     min_distance = 10) {
    res <- min_distance
    r <- raster::raster(extent(range(occurrences[, lon]),
                               range(occurrences[, lat])) + res)
    res(r) <- res
    pts <- dismo::gridSample(occurrences, r, n = 1)
    message(paste0(dim(pts)[1],
               " Points remaining after the geographic filter of ",
               min_distance, "km", "\n"))
    return(pts)
}
