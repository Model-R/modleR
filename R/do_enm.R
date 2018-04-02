#' Fits ecological niche models for various algorithms.
#'
#' @inheritParams do_bioclim
#' @param bioclim Execute bioclim from its dismo implementation
#' @param domain Execute domain from its dismo implementation
#' @param mahal Execute mahalanobis distance from its dismo implementation
#' @param maxent Execute maxent from its dismo implementation
#' @param glm Execute glm from the dismo implementation
#' @param rf Execute random forests from its dismo implementation
#' @param svm Execute svm from kernlab package
#' @param svm2 Execute svm from e1071 package
#' @return A set of ecological niche models for each partiion and algorithm,
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
                   mask,
                   write_png = FALSE,
                   n.back,
                   bioclim = TRUE,
                   domain = TRUE,
                   glm = TRUE,
                   mahal = TRUE,
                   maxent = TRUE,
                   rf = TRUE,
                   svm = TRUE,
                   svm2 = TRUE) {
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
        do_randomForest(
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
