#' Fits, predicts and evaluates ecological niche models for various algorithms
#' 
#' Runs internally \code{\link{do_any}}. Can be used to run multiple algorithms at time. See \strong{Details} in \code{\link{do_any}} for a description of how each algorithm is implemented.
#'
#' @inheritParams do_any
#' @param bioclim Execute bioclim algorithm from the \pkg{dismo} implementation with \code{\link[dismo]{bioclim}} function.
#' @param brt Execute Boosted Regression Trees with \code{\link[dismo]{gbm.step}} from \pkg{dismo}.
#' @param domain Execute domain from the \pkg{dismo} implementation with \code{\link[dismo]{domain}} function.
#' @param glm Execute GLM as suggested by the \pkg{dismo} documentation with \code{\link{glm}} and \code{\link{step}}.
#' @param mahal Execute Mahalanobis distance from the \pkg{dismo} implementation with \code{\link[dismo]{mahal}}.
#' @param maxent Execute Maxent algorithm from the \pkg{dismo} implementation with \code{\link[dismo]{maxent}} function.
#' @param maxnet Execute Maxent algorithm from the \pkg{maxnet} implementation with \code{\link[maxnet]{maxnet}} function.
#' @param rf Execute Random forests algorithm from \pkg{randomForest} package with function\code{\link[randomForest]{tuneRF}} as suggested by the \pkg{dismo} documentation.
#' @param svme Execute Support Vector Machines (SVM) algorithm from \pkg{e1071} package with \code{\link[e1071]{best.tune}} function.
#' @param svmk Execute Support Vector Machines (SVM) algorithm from \pkg{kernlab} package with \code{\link[kernlab]{ksvm}} function.
#' @param ... Any parameter from \link{do_any}
#' @return A set of ecological niche models for each partition and algorithm,
#'         written in the \code{models_dir} subfolder.
#' @author Andrea Sánchez-Tapia
#' @examples 
#' # run setup_sdmdata 
#' sp <- names(coordenadas)[1]
#' sp_coord <- coordenadas[[1]]
#' sp_setup <- setup_sdmdata(species_name=sp, occurrences=sp_coord, example_vars)
#' 
#' # run do_many
#' sp_many <- do_many(species_name=sp,
#'                      predictors=example_vars,
#'                      bioclim=TRUE, 
#'                      maxent=TRUE)
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