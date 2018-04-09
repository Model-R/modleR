#' Fits ecological niche models using any dismo algorithm.
#'
#' @param species_name A character string with the species name
#' @param algo The algorithm to be fitted \code{c("bioclim", "maxent","domain",
#'                                        "mahal", "glm", "svm.k","svm.e","rf")}
#' @param coordinates A two-column data frame with the occurrence points
#' @param partitions The number of partitions for a cross validation
#' @param buffer Defines if a buffer will be used to sample pseudo-absences
#'        (F, "mean", "median", "max")
#' @param seed For reproducibility purposes
#' @param predictors A RasterStack of predictor variables
#' @param models.dir Folder path to save the output files
#' @param project.model Logical, whether to perform a projection
#' @param projections The RasterStack of projeciton variables
#' @param mask A SpatialPolygonsDataFrame to be used to mask the final models
#' @param write_png Logical, whether png files will be written
#' @param n.back Number of pseudoabsence points
#' @return A data frame with the evaluation statistics (TSS, AUC, etc.)
#' @import grDevices
#' @importFrom utils write.table
#' @export
do_any <- function(species_name,
                   algo = c("bioclim", "maxent", "mahal", "domain"),
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
                   bootstrap = T,
                   boot_proportion = 0.8,
                   n_boot = 10,
                   crossvalidation = F,
                   n_cv = 10) {

    cat(paste(algo, "\n"))

    #sdmdatasetup
    partition.folder <- paste0(models.dir, "/", species_name, "/present", "/partitions")
    if (file.exists(paste0(partition.folder, "/sdmdata.txt"))) {
        sdmdata <- read.table(paste0(partition.folder, "/sdmdata.txt"))
    } else {
        sdmdata <- setup_sdmdata(
            species_name = species_name,
            coordinates = coordinates,
            partitions = partitions,
            buffer = buffer,
            seed = seed,
            predictors = predictors,
            models.dir = models.dir,
            plot_sdmdata = plot_sdmdata,
            n.back = n.back,
            bootstrap = bootstrap,
            boot_proportion = boot_proportion,
            n_boot = n_boot,
            crossvalidation = crossvalidation,
            n_cv = n_cv
        )
    }

    ##### Hace los modelos
    runs <- which(names(sdmdata) == "pa") - 1

    #para cada columna de la matriz de diseño
    for (i in seq_along(1:runs)) {
        group.all <- sdmdata[,i]
        group <- group.all[sdmdata$pa == 1]
        bg.grp <- group.all[sdmdata$pa == 0]
        backgr <- sdmdata[sdmdata$pa == 0, c("lon", "lat")]
        #para cada grupo
        for (g in unique(group)) {
            cat(paste(species_name,"run number", i, "partition number", g, "\n"))
            pres_train <- coordinates[group != g, ]
            if (nrow(coordinates) == 1)
                pres_train <- coordinates[group == g, ]
            pres_test <- coordinates[group == g, ]
            backg_test <- backgr[bg.grp == g, ]
            sdmdata_train <- sdmdata[group.all != g,]#presences and absences
            envtrain <-  sdmdata_train[,grep("layer", names(sdmdata_train))]
            sdmdata_test  <- sdmdata[group.all == g,]#presences and absences
            env_pres_test <- sdmdata_test[sdmdata_test$pa == 1,
                                          grep("layer", names(sdmdata_test))]
            env_backg_test <- sdmdata_test[sdmdata_test$pa == 0,
                                           grep("layer", names(sdmdata_test))]

            if (algo == "bioclim") mod <- dismo::bioclim(predictors, pres_train)
            if (algo == "maxent")  mod <- dismo::maxent(predictors, pres_train)
            if (algo == "mahal")   mod <- dismo::mahal(predictors, pres_train)
            if (algo == "domain")  mod <- dismo::domain(predictors, pres_train)
            if (algo == "rf") {
                mod <- randomForest::randomForest(sdmdata_train$pa ~ ., data = envtrain)
            }
            if (algo == "glm") {
                ## GLM solo usa a partir de sdmdata BIENs
                null.model <- glm(sdmdata_train$pa ~ 1, data = envtrain,
                                  family = "binomial")
                full.model <- glm(sdmdata_train$pa ~ ., data = envtrain,
                                  family = "binomial")
                mod <- step(object = null.model, scope = formula(full.model),
                            direction = "both", trace = F)
                }
            if (algo == "svm.k") {
            mod <- kernlab::ksvm(sdmdata_train$pa ~ ., data = envtrain)
            }
            if (algo == "svm.e") {
                mod <-
                    e1071::best.tune("svm", envtrain, sdmdata_train$pa, data = envtrain)
}

            eval_mod <- dismo::evaluate(pres_test, backg_test, mod, predictors)
           #if (algo %in% c("rf", "glm", "svm", "svm2")) {
            #eval_mod2 <- dismo::evaluate(env_pres_test, env_backg_test, mod)
            #}

            th_mod   <- eval_mod@t[which.max(eval_mod@TPR + eval_mod@TNR)]
            th_table <- dismo::threshold(eval_mod)
            mod_TSS  <- max(eval_mod@TPR + eval_mod@TNR) - 1
            mod_cont <- dismo::predict(predictors, mod, progress = "text")
            mod_bin  <- mod_cont > th_mod
            mod_cut  <- mod_cont * mod_bin
            th_table$AUC <- eval_mod@auc
            th_table$TSS <- mod_TSS
            th_table$algoritmo <- algo
            th_table$run <- i
            th_table$partition <- g
            row.names(th_table) <- paste(species_name, i, g, algo)

            write.table(
                th_table,
                file = paste0(partition.folder,
                              "/evaluate",
                              species_name,
                              "_",
                              i,
                              "_",
                              g,
                              "_",
                              algo,
                              ".txt")
            )

            if (class(mask) == "SpatialPolygonsDataFrame") {
                mod_cont <- crop_model(mod_cont, mask)
                mod_bin <- crop_model(mod_bin, mask)
                mod_cut <- crop_model(mod_cut, mask)
            }
            raster::writeRaster(
                x = mod_cont,
                filename = paste0(
                    partition.folder,
                    "/",
                    algo,
                    "_cont_",
                    species_name,
                    "_",
                    i,
                    "_",
                    g,
                    ".tif"
                ),
                overwrite = T
            )
            raster::writeRaster(
                x = mod_bin,
                filename = paste0(
                    partition.folder,
                    "/",
                    algo,
                    "_bin_",
                    species_name,
                    "_",
                    i,
                    "_",
                    g,
                    ".tif"
                ),
                overwrite = T
            )
            raster::writeRaster(
                x = mod_cut,
                filename = paste0(
                    partition.folder,
                    "/",
                    algo,
                    "_cut_",
                    species_name,
                    "_",
                    i,
                    "_",
                    g,
                    ".tif"
                ),
                overwrite = T
            )

            if (write_png == T) {
                png(paste0(
                    partition.folder,
                    "/",
                    algo,
                    "_cont",
                    species_name,
                    "_",
                    i,
                    "_",
                    g,
                    ".png"
                ))
                raster::plot(mod_cont,
                             main = paste(
                                 algo,
                                 "raw",
                                 "\n",
                                 "AUC =",
                                 round(eval_mod@auc, 2),
                                 "-",
                                 "TSS =",
                                 round(mod_TSS, 2)
                             ))
                dev.off()
                png(paste0(
                    partition.folder,
                    "/",
                    algo,
                    "_bin",
                    species_name,
                    "_",
                    i,
                    "_",
                    g,
                    ".png"
                ))
                raster::plot(mod_bin,
                             main = paste(
                                 algo,
                                 "bin",
                                 "\n",
                                 "AUC =",
                                 round(eval_mod@auc, 2),
                                 "-",
                                 "TSS =",
                                 round(mod_TSS, 2)
                             ))
                dev.off()
                png(paste0(
                    partition.folder,
                    "/",
                    algo,
                    "cut",
                    species_name,
                    "_",
                    i,
                    "_",
                    g,
                    ".png"
                ))
                raster::plot(mod_cut,
                             main = paste(
                                 algo,
                                 "cut",
                                 "\n",
                                 "AUC =",
                                 round(eval_mod@auc, 2),
                                 "-",
                                 "TSS =",
                                 round(mod_TSS, 2)
                             ))
                dev.off()
            }

            if (project.model == T) {
                for (proj in projections) {
                    projection.folder <- paste0(models.dir, "/", species_name, "/", proj)
                    if (file.exists(projection.folder) == FALSE)
                        dir.create(paste0(projection.folder), recursive = T)
                    data <- list.files(paste0("./env/", proj), pattern = proj)
                    data2 <- stack(data)
                    mod_proj <- predict(data2, mod, progress = "text")
                    mod_proj_bin <- mod_proj > th_mod
                    mod_proj_cut <- mod_proj_bin * mod_proj
                    # Normaliza o modelo cut
                    mod_proj_cut <- mod_proj_cut / maxValue(mod_proj_cut)
                    if (class(mask) == "SpatialPolygonsDataFrame") {
                        mod_proj     <- crop_model(mod_proj, mask)
                        mod_proj_bin <- crop_model(mod_proj_bin, mask)
                        mod_proj_cut <- crop_model(mod_proj_cut, mask)
                    }
                    writeRaster(
                        x = mod_proj,
                        filename = paste0(
                            projection.folder,
                            "/",
                            algo,
                            "_cont_",
                            species_name,
                            "_",
                            i,
                            "_",
                            g,
                            ".tif"
                        ),
                        overwrite = T
                    )
                    writeRaster(
                        x = mod_proj_bin,
                        filename = paste0(
                            projection.folder,
                            "/",
                            algo,
                            "_bin_",
                            species_name,
                            "_",
                            i,
                            "_",
                            g,
                            ".tif"
                        ),
                        overwrite = T
                    )
                    writeRaster(
                        x = mod_proj_cut,
                        filename = paste0(
                            projection.folder,
                            "/",
                            algo,
                            "_cut_",
                            species_name,
                            "_",
                            i,
                            "_",
                            g,
                            ".tif"
                        ),
                        overwrite = T
                    )
                    rm(data2)
                }
            }
        }
    }

    return(th_table)
}
