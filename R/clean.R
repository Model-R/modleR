#' Cleans the occurrence records
#'
#' This function performs data cleaning of occurrence records by removing
#' records outside the extension of the raster of predictor variables. It also
#' allows data cleaning of (1) NAs; (2) duplicated records; (3) duplicated
#' records in the same pixel.
#'
#' @inheritParams setup_sdmdata
#' @param clean_dupl Logical. If TRUE, removes points with the same longitude and
#'  latitude
#' @param clean_nas Logical. If TRUE, removes points that are outside the bounds
#' of the raster
#' @param clean_uni Logical. If TRUE, selects only one point per pixel
#'
#' @details Used internally in function \code{\link[modleR]{setup_sdmdata}}.
#'
#' @return A data frame containing longitude and latitude.
#'
#'
#' @seealso \code{\link[raster]{cellFromXY}}, \code{\link[raster]{mask}},
#' \code{\link[raster]{extract}}
#'
#' @examples
#' occs <- coordenadas[[1]]
#' clean(occurrences = occs,
#'       predictors = example_vars,
#'       clean_dupl = TRUE,
#'       clean_nas = TRUE,
#'       clean_uni = TRUE)
#'
#' @import raster
#'
#' @export
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
