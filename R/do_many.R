#' Fits ecological niche models for various algorithms.
#'
#' @inheritParams do_any
#' @param bioclim Execute bioclim from the dismo implementation
#' @param domain Execute domain from the dismo implementation
#' @param mahal Execute mahalanobis distance from the dismo implementation
#' @param maxent Execute maxent from the dismo implementation
#' @param maxnet Execute maxent from the maxnet implementation
#' @param glm Execute GLM as suggested by the dismo documentation
#' @param rf Execute random forests from randomForest() as suggested by the dismo documentation
#' @param svmk Execute svm from kernlab package
#' @param svme Execute svm from e1071 package
#' @param mindist Execute minimum euclidean distance
#' @param centroid Execute euclidean distance to the environmental centroid
#' @param brt Execute boosted regression trees
#' @param ... Any parameter from \link{setup_sdmdata}
#' @return A set of ecological niche models for each partition and algorithm,
#'         written in the \code{models_dir} subfolder
#' @author Andrea Sánchez-Tapia
#' @export
#'
do_many <- function(species_name,
                    occurrences,
                    predictors,
                    models_dir = "./models",
                    bioclim = FALSE,
                    domain = FALSE,
                    glm = FALSE,
                    mahal = FALSE,
                    maxent = FALSE,
                    maxnet = FALSE,
                    rf = FALSE,
                    svmk = FALSE,
                    svme = FALSE,
                    mindist = FALSE,
                    centroid = FALSE,
                    brt = FALSE,
                    ...) {

    if (bioclim == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "bioclim",
            ...)
    }
    if (domain == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "domain",
            ...)
    }
    if (glm == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "glm",
            ...)
    }
    if (mahal == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "mahal",
            ...)
    }
    if (maxent == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "maxent",
            ...)
    }
    if (maxnet == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "maxnet",
            ...)
    }
    if (rf == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "rf",
            ...)
    }
    if (svmk == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "svmk",
            ...)
    }
    if (svme == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "svme",
            ...)
    }
    if (mindist == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "mindist",
            ...)
    }
    if (centroid == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "centroid",
            ...)
    }
    if (brt == T) {
        do_any(
            species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            algo = "brt",
            ...)
    }
}
