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
#' @return A set of ecological niche models for each partition and algorithm,
#'         written in the \code{models.dir} subfolder
#' @author Andrea Sánchez-Tapia
#' @export
#'
do_enm <- function(species_name,
                   coordinates,
                   partitions,
                   buffer = FALSE,
                   seed = 512,
                   predictors,
                   models.dir = "./models",
                   project.model = FALSE,
                   projections = NULL,
                   mask = NULL,
                   write_png = FALSE,
                   n.back,
                   bioclim = FALSE,
                   domain = FALSE,
                   glm = FALSE,
                   mahal = FALSE,
                   maxent = FALSE,
                   rf = FALSE,
                   svm.k = FALSE,
                   svm.e = FALSE) {

    if (bioclim == T) {
        do_any(
            species_name,
            algo = "bioclim",
            coordinates = coordinates,
            partitions = partitions,
            buffer = buffer,
            seed = seed,
            predictors = predictors,
            models.dir = models.dir,
            project.model = project.model,
            projections = projections,
            mask = mask,
            n.back = n.back,
            write_png = write_png)
    }
    if (domain == T) {
        do_any(
            species_name,
            algo = "domain",
            coordinates = coordinates,
            partitions = partitions,
            buffer = buffer,
            seed = seed,
            predictors = predictors,
            models.dir = models.dir,
            project.model = project.model,
            projections = projections,
            mask = mask,
            n.back = n.back,
            write_png = write_png)
    }
    if (glm == T) {
        do_any(
            species_name,
            algo = "glm",
            coordinates = coordinates,
            partitions = partitions,
            buffer = buffer,
            seed = seed,
            predictors = predictors,
            models.dir = models.dir,
            project.model = project.model,
            projections = projections,
            mask = mask,
            n.back = n.back,
            write_png = write_png)
    }
    if (mahal == T) {
        do_any(
            species_name,
            algo = "mahal",
            coordinates = coordinates,
            partitions = partitions,
            buffer = buffer,
            seed = seed,
            predictors = predictors,
            models.dir = models.dir,
            project.model = project.model,
            projections = projections,
            mask = mask,
            n.back = n.back,
            write_png = write_png)
    }
    if (maxent == T) {
        do_any(
            species_name,
            algo = "maxent",
            coordinates = coordinates,
            partitions = partitions,
            buffer = buffer,
            seed = seed,
            predictors = predictors,
            models.dir = models.dir,
            project.model = project.model,
            projections = projections,
            mask = mask,
            n.back = n.back,
            write_png = write_png)
    }
    if (rf == T) {
        do_any(
            species_name,
            algo = "rf",
            coordinates = coordinates,
            partitions = partitions,
            buffer = buffer,
            seed = seed,
            predictors = predictors,
            models.dir = models.dir,
            project.model = project.model,
            projections = projections,
            mask = mask,
            n.back = n.back,
            write_png = write_png)
    }
    if (svm.k == T) {
        do_any(
            species_name,
            algo = "svm.k",
            coordinates = coordinates,
            partitions = partitions,
            buffer = buffer,
            seed = seed,
            predictors = predictors,
            models.dir = models.dir,
            project.model = project.model,
            projections = projections,
            mask = mask,
            n.back = n.back,
            write_png = write_png)
    }
    if (svm.e == T) {
        do_any(
            species_name,
            algo = "svm.e",
            coordinates = coordinates,
            partitions = partitions,
            buffer = buffer,
            seed = seed,
            predictors = predictors,
            models.dir = models.dir,
            project.model = project.model,
            projections = projections,
            mask = mask,
            n.back = n.back,
            write_png = write_png)
    }
}
