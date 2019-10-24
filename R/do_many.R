#' @rdname model_fit
#' @export
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
                    brt = FALSE,
                    ...) {

  if (bioclim == TRUE) {
    do_any(
      species_name,
      algorithm = "bioclim",
      ...)
  }
  if (domain == TRUE) {
    do_any(
      species_name,
      algorithm = "domain",
      ...)
  }
  if (glm == TRUE) {
    do_any(
      species_name,
      algorithm = "glm",
      ...)
  }
  if (mahal == TRUE) {
    do_any(
      species_name,
      algorithm = "mahal",
      ...)
  }
  if (maxent == TRUE) {
    do_any(
      species_name,
      algorithm = "maxent",
      ...)
  }
  if (maxnet == TRUE) {
    do_any(
      species_name,
      algorithm = "maxnet",
      ...)
  }
  if (rf == TRUE) {
    do_any(
      species_name,
      algorithm = "rf",
      ...)
  }
  if (svmk == TRUE) {
    do_any(
      species_name,
      algorithm = "svmk",
      ...)
  }
  if (svme == TRUE) {
    do_any(
      species_name,
      algorithm = "svme",
      ...)
  }
  if (brt == TRUE) {
    do_any(
      species_name,
      algorithm = "brt",
      ...)
  }
}
