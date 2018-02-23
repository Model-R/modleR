#' Fits ecological niche models for various algorithms.
#'
#' @param sp A character string with the species name
#' @param coordinates A three-column data frame with the occurrence points of all species (sp, lat, lon)
#' @param buffer Defines if a buffer will be used to sample pseudo-absences (F, "mean", "median", "max")
#' @param seed For reproducibility purposes
#' @param predictors A RasterStack of predictor variables
#' @param models.dir Folder path to save the output files
#' @param mask A SpatialPolygonsDataFrame to be used as a mask to cut the final models
#' @param write_png Logical, whether png files will be written, defaults to F
#' @param n.back Number of pseudoabsence points
#' @return A data frame with the evaluation statistics (TSS, AUC, and their respective thresholds)
#' @export
#'
do_enm <- function(especie,
                      coordinates,
                      partitions,
                      buffer,
                      seed,
                      predictors,
                      models.dir,
                      project.model,
                      projections,
                      mask,
                      n.back,
                      write_png,
                      maxent = T,
                      domain = T,
                      bioclim = T,
                      rf = T,
                      glm = T,
                      svm = T,
                      svm2 = T,
                      mahal = T
                      ) {
    ocorrencias <-
        coordinates[coordinates$sp == especie, c("lon", "lat")]
    if (maxent == T) {
        do_maxent(
            especie,
            coordinates = ocorrencias,
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
            especie,
            coordinates = ocorrencias,
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
    if (bioclim == T) {
        do_bioclim(
            especie,
            coordinates = ocorrencias,
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
            especie,
            coordinates = ocorrencias,
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
            especie,
            coordinates = ocorrencias,
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
            especie,
            coordinates = ocorrencias,
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
            especie,
            coordinates = ocorrencias,
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
            especie,
            coordinates = ocorrencias,
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
