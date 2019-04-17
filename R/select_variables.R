#' Helper function to select environmental variables from a stack.
#'
#' This function takes a stack of environmental variables and the calibration area. The function calculates Pearson correlations between all environmental variables and returns a new stack with uncorrelated variables based on a cutoff.
#'
#' @param species_name A character string with the species name
#' @param predictors Stack of environmental variables
#' @param buffer Raster specified by the user or output from create_buffer()
#' @param cutoff Cutoff value of correlation between variables to exclude environmental layer. Default is to exclude environmental variables with correlation > 0.8.
#' @param percent percentage of the raster values to be sampled to calculate the correlation. Defaults to 0.8 but should be useful with high resolution rasters
#' @param ... parameters from create_buffer()
#' @return A raster stack of independent environmental variables based on a
#'  specific cutoff
#' @author Andrea Sánchez-Tapia and Sara Mortara
#' @seealso \code{\link[ModelR]{create_buffer}}
#' @import raster
#' @importFrom stats cor
#' @export
#' @examples
#'
#' ## selecting data for only sp1
#' coord1sp <- coordenadas[coordenadas$sp == unique(coordenadas$sp)[1],]
#' ## selecting only columns with longitude and latitude
#' occ <- coord1sp[,c(2,3)]
#' # using coord1sp to create buffer w/ mean distance between points
#' buf <- create_buffer("foo", occ, example_vars)
#' # running select_variables w/ output from create_buffer
#' select_variables("foo", predictors = example_vars, buffer = buf)
select_variables <- function(species_name,
                             predictors,
                             buffer = NULL,
                             cutoff = 0.8,
                             percent = 0.8,
                             ...) {

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
  exclude.vars <- caret::findCorrelation(cor(vals), cutoff = cutoff)
  if (length(exclude.vars) > 0) {
      excluded <- names(predictors)[exclude.vars]
      retained <- setdiff(names(predictors), excluded)
      final_vars <- raster::subset(predictors, retained, drop = F)
      message(paste(paste(excluded, collapse = ","),
                    "excluded with cutoff =", cutoff))
      } else {
          final_vars <- predictors
          message(paste("No variables were excluded with cutoff =", cutoff))
          }
      return(final_vars)
}
