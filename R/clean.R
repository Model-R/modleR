clean <- function(occurrences,
                  lon = "lon",
                  lat = "lat",
                  predictors,
                  clean_dupl = FALSE,
                  clean_nas = FALSE,
                  clean_uni = FALSE) {
# define occurrences to be lon and lat
  occurrences <- occurrences[, c(lon, lat)]
    if (exists("predictors")) {
        ori <- nrow(occurrences)
        if (clean_dupl == TRUE) {
            message("cleaning duplicates")
            dupls <- base::duplicated(occurrences)
            if (any(dupls)) {
            occurrences <- occurrences[!dupls, ]
            }
        }
        if (clean_nas == TRUE) {
            message("cleaning occurrences with no environmental data")
            presvals <- raster::extract(predictors, occurrences)
            compl <- complete.cases(presvals)
            if (any(compl)) {
            occurrences <- occurrences[compl, ]
            }
        }
        if (clean_uni == TRUE) {
            if (clean_nas == FALSE) {
                warning("There may be points that are outside the raster")
            }
            mask <- predictors[[1]]
            cell <- raster::cellFromXY(mask, occurrences)
            dup <- duplicated(cell)
            if (any(dup)) {
                occurrences <- occurrences[!dup, ]
                }
        }

        cat(ori - nrow(occurrences), "points removed\n")
        cat(nrow(occurrences), " clean points\n")
        return(occurrences)
    } else
        (cat("Indicate the object with the predictive variables"))
}
