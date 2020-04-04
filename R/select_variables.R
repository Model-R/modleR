select_variables <- function(predictors,
                             buffer = NULL,
                             cutoff = 0.8,
                             sample_proportion = 0.8) {

    if (!class(predictors) %in% c("RasterBrick","RasterStack")) {
  stop("predictors must be a RasterBrick or RasterStack object")
    }
    if (!is.null(buffer) & class(buffer) %in% c("RasterBrick", "RasterStack")) {
        predictors <- crop_model(predictors, buffer)
    }

  sampled <- dismo::randomPoints(
      mask = predictors,
      n = floor(sum(!is.na(raster::values(predictors[[1]]))) * sample_proportion))
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
