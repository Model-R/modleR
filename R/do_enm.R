#' Fits ecological niche models for various algorithms.
#'
#' @inheritParams do_any
#' @param bioclim Execute bioclim from the dismo implementation
#' @param domain Execute domain from the dismo implementation
#' @param mahal Execute mahalanobis distance from the dismo implementation
#' @param maxent Execute maxent from the dismo implementation
#' @param glm Execute GLM as suggested by the dismo documentation
#' @param rf Execute random forests from randomForest() as suggested
#' by the dismo documentation
#' @param svm.k Execute svm from kernlab package
#' @param svm.e Execute svm from e1071 package
#' @param mindist Execute minimum euclidean distance
#' @param centroid Execute euclidean distance to the environmental centroid
#' @return A set of ecological niche models for each partition and algorithm,
#'         written in the \code{models_dir} subfolder
#' @author Andrea Sánchez-Tapia
#' @export
#'
do_enm <- function(species_name,
                   coordinates,
                   predictors,
                   models_dir = "./models",
                   bioclim = FALSE,
                   domain = FALSE,
                   glm = FALSE,
                   mahal = FALSE,
                   maxent = FALSE,
                   rf = FALSE,
                   svm.k = FALSE,
                   svm.e = FALSE,
                   mindist = FALSE,
                   centroid = FALSE,
                   ...) {

    if (bioclim == T) {
        do_any(
            species_name,
            algo = "bioclim",
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (domain == T) {
        do_any(
            species_name,
            algo = "domain",
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (glm == T) {
        do_any(
            species_name,
            algo = "glm",
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (mahal == T) {
        do_any(
            species_name,
            algo = "mahal",
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (maxent == T) {
        do_any(
            species_name,
            algo = "maxent",
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (rf == T) {
        do_any(
            species_name,
            algo = "rf",
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (svm.k == T) {
        do_any(
            species_name,
            algo = "svm.k",
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (svm.e == T) {
        do_any(
            species_name,
            algo = "svm.e",
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (mindist == T) {
        do_any(
            species_name,
            algo = "mindist",
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (centroid == T) {
        do_any(
            species_name,
            algo = "centroid",
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
}
