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
#' @param ... Any parameter from \link{do_any}
#' @return A set of ecological niche models for each partition and algorithm,
#'         written in the \code{models_dir} subfolder
#' @author Andrea Sánchez-Tapia
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
                    mindist = FALSE,
                    centroid = FALSE,
                    brt = FALSE,
                    ...) {

    if (bioclim == T) {
        do_any(
            species_name,
            algo = "bioclim",
            ...)
    }
    if (domain == T) {
        do_any(
            species_name,
            algo = "domain",
            ...)
    }
    if (glm == T) {
        do_any(
            species_name,
            algo = "glm",
            ...)
    }
    if (mahal == T) {
        do_any(
            species_name,
            algo = "mahal",
            ...)
    }
    if (maxent == T) {
        do_any(
            species_name,
            algo = "maxent",
            ...)
    }
    if (maxnet == T) {
        do_any(
            species_name,
            algo = "maxnet",
            ...)
    }
    if (rf == T) {
        do_any(
            species_name,
            algo = "rf",
            ...)
    }
    if (svmk == T) {
        do_any(
            species_name,
            algo = "svmk",
            ...)
    }
    if (svme == T) {
        do_any(
            species_name,
            algo = "svme",
            ...)
    }
    if (mindist == T) {
        do_any(
            species_name,
            algo = "mindist",
            ...)
    }
    if (centroid == T) {
        do_any(
            species_name,
            algo = "centroid",
            ...)
    }
    if (brt == T) {
        do_any(
            species_name,
            algo = "brt",
            ...)
    }
}
