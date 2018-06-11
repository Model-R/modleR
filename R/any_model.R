#' Fits ecological niche models using any dismo algorithm.
#'
#' @inheritParams setup_sdmdata
#' @param algo The algorithm to be fitted \code{c("bioclim", "maxent","domain",
#'                                        "mahal", "glm", "svm.k","svm.e","rf")}
#' @param project_model Logical, whether to perform a projection
#' @param projections The RasterStack of projeciton variables
#' @param mask A SpatialPolygonsDataFrame to be used to mask the final models
#' @param write_png Logical, whether png files will be written
#' @return A data frame with the evaluation statistics (TSS, AUC, etc.)
#' @author Andrea Sánchez-Tapia
#' @seealso \code{\link[dismo]{bioclim}}
#' @seealso \code{\link[dismo]{maxent}}
#' @seealso \code{\link[dismo]{domain}}
#' @seealso \code{\link[dismo]{mahal}}
#' @import raster
#' @import grDevices
#' @importFrom utils write.table
#' @importFrom stats complete.cases formula glm step dist
#' @export
do_any <- function(species_name,
                   coordinates,
                   predictors,
                   models_dir = "./models",
                   algo = c("bioclim"), #um só
                   project_model = FALSE,
                   projections = NULL,
                   mask = NULL,
                   write_png = FALSE,
                   ...) {
    message(paste(algo, "\n"))

    #sdmdatasetup
    partition.folder <- paste0(models_dir, "/", species_name, "/present", "/partitions")

        sdmdata <- setup_sdmdata(
            species_name = species_name,
            coordinates = coordinates,
            predictors = predictors,
            models_dir = models_dir,
            ...)

    ##### Hace los modelos
    runs <- which(names(sdmdata) == "pa") - 1

    #para cada columna de la matriz de diseño
    for (i in seq_along(1:runs)) {
        group.all <- sdmdata[,i]
        group <- group.all[sdmdata$pa == 1]
        bg.grp <- group.all[sdmdata$pa == 0]
        backgr <- sdmdata[sdmdata$pa == 0, c("lon", "lat")]
        #para cada grupo
        for (g in setdiff(unique(group), 0)) {
            #excluding the zero allows for bootstrap. only 1 partition will run
            cat(paste(species_name, algo, "run number", i, "partition number",
                      g, "\n"))
            pres_train <- coordinates[group != g, ]
            if (nrow(coordinates) == 1)
                pres_train <- coordinates[group == g, ]
            pres_test <- coordinates[group == g, ]
            backg_test <- backgr[bg.grp == g, ]
            sdmdata_train <- sdmdata[group.all != g,]#presences and absences
            envtrain <-  sdmdata_train[,grep("layer", names(sdmdata_train))]
            sdmdata_test  <- sdmdata[group.all == g,]#presences and absences

            if (algo == "bioclim") mod <- dismo::bioclim(predictors, pres_train)
            if (algo == "maxent")  mod <- dismo::maxent(predictors, pres_train)
            if (algo == "mahal")   mod <- dismo::mahal(predictors, pres_train)
            if (algo == "domain")  mod <- dismo::domain(predictors, pres_train)
            if (algo == "rf") {
                mod <- randomForest::randomForest(sdmdata_train$pa ~ .,
                                                  data = envtrain)
            }
            if (algo == "glm") {
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
                mod <- e1071::best.tune("svm", envtrain, sdmdata_train$pa,
                                        data = envtrain)
                }
            if (algo %in% c("centroid", "mindist")) {
                ec_cont <- predictors[[1]]
                ec_cont[!is.na(ec_cont)] <- 1
                predictors.st <- raster::scale(predictors)
                pres.vals <- raster::extract(predictors.st, pres_train)#extrae valores
                pred.vals <- raster::values(predictors.st)
                dist.vals <- vector(length = nrow(pred.vals))

                #esto es centroid
                if (algo == "centroid") {
                    cat(paste("Euclidean environmental distance to centroid",'\n'))
                    #calcula la media ambiental de los puntos de train
                    centroid.val <- apply(pres.vals, 2, mean, na.rm = TRUE)
                    #calcula la distancia a ese centroide
                    dist.vals <- apply(pred.vals, 1, FUN = function(x) {
                        dist(rbind(centroid.val, x))
                        }
                        )
                raster::values(ec_cont) <- -dist.vals
                    }

                #esto es mindist
                #no está haciendo bien el cut
                if (algo == "mindist") {
                    pred.vals.cc <- pred.vals[complete.cases(pred.vals), ]
                    min.vals <- vector(length = nrow(pred.vals.cc))
                    for (m in 1:nrow(pred.vals.cc)) {
                        mindata <- rbind(pred.vals.cc[m, ], pres.vals)
                        min.vals[m] <- min(as.matrix(dist(mindata))[1,][-1], na.rm = T)
                    }
                ec_cont[complete.cases(pred.vals)] <- -min.vals

                }
                #rescales the output so the output is between
                if (raster::minValue(ec_cont) < 0) {
                    ec_cont <-
                        (ec_cont - raster::minValue(ec_cont)) /
                        raster::maxValue(ec_cont - raster::minValue(ec_cont))
                }
                # evaluates the model to extract LPT
                p <- raster::extract(ec_cont, y = pres_test)
                a <- raster::extract(ec_cont, y = backg_test)
                eec <- dismo::evaluate(p = p, a = a)
                if (!nrow(coordinates) %in% c(1, 2)) {
                    #sólo corta por LTP si hay más de dos puntos...
                    LPTec <- dismo::threshold(eec, 'no_omission')
                    ec_cont[ec_cont < LPTec] <- LPTec
                    ec_cont <- (ec_cont - raster::minValue(ec_cont)) /
                        raster::maxValue(ec_cont - raster::minValue(ec_cont))
                    #ec_cont[ec_cont < LPTec] <- 0 #éste debe ser el problema
                }
                #evaluates again because the values changed
                p <- raster::extract(ec_cont, y = pres_test)
                a <- raster::extract(ec_cont, y = backg_test)
                eec <- dismo::evaluate(p = p, a = a)
        }

            if (algo %in% c("mindist", "centroid")) {
                eval_mod <- eec
                mod_cont <- ec_cont
            } else {
                eval_mod <- dismo::evaluate(pres_test, backg_test, mod, predictors)
                mod_cont <- dismo::predict(predictors, mod, progress = "text")
            }

            th_mod   <- eval_mod@t[which.max(eval_mod@TPR + eval_mod@TNR)]
            th_table <- dismo::threshold(eval_mod)
            mod_TSS  <- max(eval_mod@TPR + eval_mod@TNR) - 1

            mod_bin  <- mod_cont > th_mod
            mod_cut  <- mod_cont * mod_bin
            th_table$AUC <- eval_mod@auc
            th_table$TSS <- mod_TSS
            th_table$algoritmo <- algo
            th_table$run <- i
            th_table$partition <- g
            row.names(th_table) <- paste(species_name, i, g, algo)

            write.table(th_table, file = paste0(partition.folder, "/evaluate_",
                                                species_name, "_", i, "_", g,
                                                "_", algo, ".txt"))

            if (class(mask) == "SpatialPolygonsDataFrame") {
                mod_cont <- crop_model(mod_cont, mask)
                mod_bin <- crop_model(mod_bin, mask)
                mod_cut <- crop_model(mod_cut, mask)
            }
            raster::writeRaster(x = mod_cont,
                                filename = paste0(partition.folder, "/", algo,
                                                  "_cont_", species_name, "_",
                                                  i, "_", g, ".tif"),
                                overwrite = T)
            raster::writeRaster(x = mod_bin,
                                filename = paste0(partition.folder,  "/", algo,
                                                  "_bin_", species_name, "_",
                                                  i, "_", g, ".tif"),
                                overwrite = T)
            raster::writeRaster(x = mod_cut,
                                filename = paste0(partition.folder, "/", algo,
                                                  "_cut_", species_name, "_",
                                                  i, "_", g, ".tif"),
                overwrite = T)

            if (write_png == T) {
                png(paste0(partition.folder, "/", algo, "_cont_", species_name,
                           "_", i, "_", g, ".png"))
                raster::plot(mod_cont,
                             main = paste(algo, "raw", "\n", "AUC =",
                                          round(eval_mod@auc, 2), "-", "TSS =",
                                          round(mod_TSS, 2)))
                dev.off()
                png(paste0(partition.folder, "/", algo, "_bin_", species_name,
                           "_", i, "_", g, ".png"))
                raster::plot(mod_bin,
                             main = paste(algo, "bin", "\n", "AUC =",
                                          round(eval_mod@auc, 2), "-", "TSS =",
                                          round(mod_TSS, 2)))
                dev.off()
                png(paste0(partition.folder, "/", algo, "_cut_", species_name,
                           "_", i, "_", g, ".png"))
                raster::plot(mod_cut,
                             main = paste(algo, "cut", "\n", "AUC =",
                                          round(eval_mod@auc, 2), "-", "TSS =",
                                          round(mod_TSS, 2)))
                dev.off()
            }

            if (project_model == T) {
                for (proj in projections) {
                    projection.folder <- paste0(models_dir, "/", species_name,
                                                "/", proj)
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
                    writeRaster(x = mod_proj,
                                filename = paste0(projection.folder, "/", algo,
                                                  "_cont_", species_name, "_",
                                                  i, "_", g, ".tif"),
                                overwrite = T)
                    writeRaster(x = mod_proj_bin,
                                filename = paste0(projection.folder, "/", algo,
                                                  "_bin_", species_name, "_",
                                                  i, "_", g, ".tif"),
                                overwrite = T)
                    writeRaster(x = mod_proj_cut,
                                filename = paste0(projection.folder, "/", algo,
                                                  "_cut_", species_name, "_", i,
                                                  "_", g, ".tif"),
                                overwrite = T)
                    rm(data2)
                }
            }
        }
    }

    return(th_table)
}
