#' Excludes occurrence points within a geographic distance.
#'
#' @param occurrences data.frame. Table with the species" data. It should contain only two columns: lon and lat, in that order.
#' @param min_distance numeric. Minimum distance (in Km) between points.
#' @return Table of occurrence points with minimum distance indicated in min.distance
#' @author Diogo S. B. Rocha
#' @references Varela, S., Anderson, R. P., García-Valdés, R., &
#' Fernández-González, F. (2014). Environmental filters reduce the effects of
#' sampling bias and improve predictions of ecological niche models.
#' Ecography, 37(11), 1084-1091. doi:10.1111/j.1600-0587.2013.00441.x
#' @seealso \code{\link[dismo]{gridSample}}
#'
#' @importFrom dismo gridSample
#' @export
geo_filt <- function(occurrences,
                     min_distance = 10) {
    res <- min_distance * 0.008333333
    r <- raster::raster(extent(range(occurrences[, 1]),
                               range(occurrences[, 2])) + res)
    res(r) <- res
    pts <- dismo::gridSample(occurrences, r, n = 1)
    message(paste0(dim(pts)[1],
               " Points remaining after the geographic filter of ",
               min_distance, "km", "\n"))
    return(pts)
}
