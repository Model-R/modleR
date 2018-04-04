#' Fits ecological niche models for various algorithms.
#'
#' @inheritParams do_bioclim
#' @param bioclim Execute bioclim from the dismo implementation
#' @param domain Execute domain from the dismo implementation
#' @param mahal Execute mahalanobis distance from the dismo implementation
#' @param maxent Execute maxent from the dismo implementation
#' @param glm Execute GLM as suggested by the dismo documentation
#' @param rf Execute random forests from randomForest() as suggested by the dismo documentation
#' @param svm Execute svm from kernlab package
#' @param svm2 Execute svm from e1071 package
#' @return A set of ecological niche models for each partition and algorithm,
#'         written in the \code{models.dir} subfolder
#' @export
#'
do_enm <- function(sp,
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
                   svm = FALSE,
                   svm2 = FALSE) {
    #ocorrencias <- coordinates[coordinates$sp == sp, c("lon", "lat")]

    if (bioclim == T) {
        do_bioclim(
            sp,
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
        do_domain(
            sp,
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
        do_GLM(
            sp,
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
        do_mahal(
            sp,
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
        do_maxent(
            sp,
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
        do_rf(
            sp,
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
    if (svm == T) {
        do_SVM(
            sp,
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
    if (svm2 == T) {
        do_SVM2(
            sp,
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
