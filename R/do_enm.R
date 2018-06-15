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
#' @param ... Any parameter from \link{setup_sdmdata}
#' @return A set of ecological niche models for each partition and algorithm,
#'         written in the \code{models_dir} subfolder
#' @author Andrea Sánchez-Tapia
#' @export
#'
do_enm <- function(species_name,
                   occurrences,
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
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (domain == T) {
        do_any(
            species_name,
            algo = "domain",
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (glm == T) {
        do_any(
            species_name,
            algo = "glm",
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (mahal == T) {
        do_any(
            species_name,
            algo = "mahal",
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (maxent == T) {
        do_any(
            species_name,
            algo = "maxent",
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (rf == T) {
        do_any(
            species_name,
            algo = "rf",
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (svm.k == T) {
        do_any(
            species_name,
            algo = "svm.k",
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (svm.e == T) {
        do_any(
            species_name,
            algo = "svm.e",
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (mindist == T) {
        do_any(
            species_name,
            algo = "mindist",
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
    if (centroid == T) {
        do_any(
            species_name,
            algo = "centroid",
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            ...)
    }
}
