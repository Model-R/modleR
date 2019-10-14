#' Automatically selects environmental variables from a stack
#'
#' This function takes a stack of environmental variables and the calibration area.
#' The function calculates Pearson correlations between all environmental variables and
#' returns a new stack with uncorrelated variables based on a cutoff.
#'
#' @inheritParams setup_sdmdata
#' @param buffer Raster specified by the user or output from create_buffer()
#' @param cutoff Cutoff value of correlation between variables to exclude environmental layer.
#' Default is to exclude environmental variables with correlation > 0.8
#' @param percent percentage of the raster values to be sampled to calculate the correlation.
#' Default is to 0.8 but should be useful with high resolution rasters
#' @return A RasterStack of independent environmental variables based on a
#'  specific cutoff
#' @seealso \code{\link[modleR]{create_buffer}}
#' @import raster
#' @importFrom stats cor
#' @export
#' @examples
#'
#' ## selecting data for only sp1
#' coord1sp <- coordenadas[[1]]
#' ## selecting only columns with longitude and latitude
#' occ <- coord1sp[,c(2,3)]
#' ## using coord1sp to create buffer w/ mean distance between points
#' buf <- create_buffer("foo", occ, predictors = example_vars)
#' ## running select_variables w/ output from create_buffer
#' select_variables(predictors = example_vars, buffer = buf)
#'
#' ## selecting variables for the entire area with a different cutoff
#' select_variables(example_vars, cutoff = 0.5)
select_variables <- function(predictors,
                             buffer = NULL,
                             cutoff = 0.8,
                             percent = 0.8) {

    if (!class(predictors) %in% c("RasterBrick","RasterStack")) {
  stop("predictors must be a RasterBrick or RasterStack object")
    }
    if (!is.null(buffer) & class(buffer) %in% c("RasterBrick", "RasterStack")) {
        predictors <- crop_model(predictors, buffer)
    }

  sampled <- dismo::randomPoints(
      mask = predictors,
      n = floor(sum(!is.na(raster::values(predictors[[1]]))) * percent))
  vals <- raster::extract(x = predictors, sampled)
  vals <- vals[complete.cases(vals),]
  exclude_vars <- caret::findCorrelation(cor(vals), cutoff = cutoff)
  if (length(exclude_vars) > 0) {
      excluded <- names(predictors)[exclude_vars]
      retained <- setdiff(names(predictors), excluded)
      final_vars <- raster::subset(predictors, retained, drop = FALSE)
      message(paste(paste(excluded, collapse = ","),
                    "excluded with cutoff =", cutoff))
      } else {
          final_vars <- predictors
          message(paste("No variables were excluded with cutoff =", cutoff))
          }
      return(final_vars)
}
