#' @rdname fit

#' @examples
#' # run setup_sdmdata
#' sp <- names(coordenadas)[1]
#' sp_coord <- coordenadas[[1]]
#' sp_setup <- setup_sdmdata(species_name = sp, occurrences = sp_coord, example_vars)
#'
#' # run do_many
#' sp_many <- do_many(species_name = sp,
#'                      predictors = example_vars,
#'                      bioclim = TRUE,
#'                      maxent = TRUE)
#' @export
#'
do_many <- function(species_name,
                    bioclim = FALSE,
                    domain = FALSE,
                    glm = FALSE,
                    mahal = FALSE,
                    maxent = FALSE,
                    maxnet = FALSE,
                    rf = FALSE,
                    svmk = FALSE,
                    svme = FALSE,
#                    mindist = FALSE,
#                    centroid = FALSE,
                    brt = FALSE,
                    ...) {

  if (bioclim == TRUE) {
    do_any(
      species_name,
      algo = "bioclim",
      ...)
  }
  if (domain == TRUE) {
    do_any(
      species_name,
      algo = "domain",
      ...)
  }
  if (glm == TRUE) {
    do_any(
      species_name,
      algo = "glm",
      ...)
  }
  if (mahal == TRUE) {
    do_any(
      species_name,
      algo = "mahal",
      ...)
  }
  if (maxent == TRUE) {
    do_any(
      species_name,
      algo = "maxent",
      ...)
  }
  if (maxnet == TRUE) {
    do_any(
      species_name,
      algo = "maxnet",
      ...)
  }
  if (rf == TRUE) {
    do_any(
      species_name,
      algo = "rf",
      ...)
  }
  if (svmk == TRUE) {
    do_any(
      species_name,
      algo = "svmk",
      ...)
  }
  if (svme == TRUE) {
    do_any(
      species_name,
      algo = "svme",
      ...)
  }
  # if (mindist == TRUE) {
  #   do_any(
  #     species_name,
  #     algo = "mindist",
  #     ...)
  # }
  # if (centroid == TRUE) {
  #   do_any(
  #     species_name,
  #     algo = "centroid",
  #     ...)
  # }
  if (brt == TRUE) {
    do_any(
      species_name,
      algo = "brt",
      ...)
  }
}
