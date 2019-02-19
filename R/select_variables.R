#' Helper function to select environmental variables from a stack
#'
#' This function takes a stack of environmental variables and the calibration area. The function calculates Pearson correlations between all environmental variables and returns a new stack with uncorrelated variables based on a cutoff.
#'
#' @param species_name A character string with the species name
#' @param predictors Stack of environmental variables
#' @param models_dir Character name of output directory where the
#'  buffer created with create_buffer() is
#' @param buffer Shapefile specified by user or use NULL to use the output from create_buffer()
#' @param cutoff Cutoff value of correlation between variables to exclude environmental layer. Default is to exclude environmental variables with correlation > 0.8.
#' @param ... parameters from create_buffer()
#' @return A raster stack of independent environmental variables based on a
#'  specific cutoff
#' @author Andrea Sánchez-Tapia and Sara Mortara
#' @seealso \code{\link[ModelR]{create_buffer}}
#' @import raster
#' @export
#' @examples
#' # using shapefile from create_buffer()
#' ## selecting data for only sp1
#' coord1sp <- coordenadas[coordenadas$sp == unique(coordenadas$sp)[1],]
#' ## selecting only columns with longitude and latitude
#' occ <- coord1sp[,c(2,3)]
#' ## using coord1sp to create buffer w/ mean distance between points
#' buf <- create_buffer(occ, buffer_type="mean", example_vars)
#' # running select_variables w/ output from create_buffer
#' select_variables(predictors=example_vars)
select_variables <- function(species_name = species_name,
                             predictors = predictors,
                             models_dir = models_dir,
                             buffer = NULL,
                             cutoff = 0.8,
                             ...) {
    partition.folder <- paste0(models_dir, "/", species_name, "/present", "/partitions")
    if (is.null(buffer)) {
        # checks if there is a buffer.shp file in partition folder
        if (list.files(path = partition.folder, pattern = "*.shp", full.names = TRUE) != paste0(partition.folder, "/", "buffer.shp")) {
stop("could not find buffer.shp in output directory, specify your own buffer")
            }
        }
    else buffer <- rgdal::readOGR(paste0(partition.folder, "/", "buffer.shp"))
if (class(predictors) != "RasterBrick") {
  stop("predictors must be a RasterBrick object")
    }
    else crop.buffer <- crop(predictors, buffer) # cutting raster by shapefile
# correlation matrix for all layers
    cor.layer <- raster::layerStats(crop.buffer, "pearson", na.rm = TRUE)$`pearson correlation coefficient`
# selecting only values > cutoff
  cor.layer[upper.tri(cor.layer, diag = TRUE)] <- 0
# if no correlation value > cutoff, returns predictors
  if (sum(abs(cor.layer) > cutoff) == 0) {
    final_vars <- predictors
     selected_vars <- list(final_vars,
                           message("all variables retained")) }
# excluding variables w/ cor > cutoff
  else exclude.vars <- which(abs(cor.layer) > cutoff, arr.ind = TRUE)
  names.vars <-
      names(predictors)[!names(predictors) %in% row.names(exclude.vars)]
    final_vars <- predictors[[names.vars]]
    selected_vars <- list(final_vars,
                          message(paste0("excluded variable(s) on layers ",
                                         row.names(exclude.vars))))
    return(selected_vars)
}
