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
