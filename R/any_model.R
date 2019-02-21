#' Fits ecological niche models using several algorithms.
#'
#' @inheritParams setup_sdmdata
#' @param algo The algorithm to be fitted \code{c("bioclim", "maxent","domain",
#'                                        "mahal", "glm", "svm.k", "svm.e",
#'                                         "rf", "brt", "mindist", "centroid")}
#' @param project_model Logical, whether to perform a projection
#' @param proj_data_folder the path to projections -containing one or more
#'  folders with the projection datasets, ex. "./env/proj/proj1"
#' @param mask A SpatialPolygonsDataFrame to be used to mask the final models
#' @param write_bin_cut Logical, whether binary and cut model files(.tif, .png) should be written
#' @param write_png Logical, whether png files will be written
#' @param conf_mat Logical, whether confusion tables should be written in the HD
#' @param ... Any parameter from \link{setup_sdmdata}
#' @return A data frame with the evaluation statistics (TSS, AUC, etc.)
#' @author Andrea Sánchez-Tapia
#' @seealso \code{\link[dismo]{bioclim}}
#' @seealso \code{\link[dismo]{maxent}}
#' @seealso \code{\link[dismo]{domain}}
#' @seealso \code{\link[dismo]{mahal}}
#' @import raster
#' @import grDevices
#' @importFrom utils write.table
#' @importFrom maxnet maxnet
#' @importFrom stats complete.cases formula glm step dist
#' @export
do_any <- function(species_name,
                   occurrences,
                   predictors,
                   models_dir = "./models",
                   algo = c("bioclim"), #um só
                   project_model = FALSE,
                   proj_data_folder = "./data/proj",
                   #proj_cut = NULL,
                   mask = NULL,
                   write_png = FALSE,
                   write_bin_cut = TRUE,
                   buffer_type = NULL,
                   dist_buf = NULL,
                   conf_mat = TRUE,
                   equalize = TRUE,
                   ...) {
    message(paste(algo, "\n"))

    #sdmdatasetup
    partition.folder <- paste0(models_dir, "/", species_name, "/present", "/partitions")

        sdmdata <- setup_sdmdata(
            species_name = species_name,
            occurrences = occurrences,
            predictors = predictors,
            models_dir = models_dir,
            buffer_type = buffer_type,
            dist_buf = dist_buf,
            equalize = equalize,
            ...)

    ##### Hace los modelos
    runs <- which(names(sdmdata) == "pa") - 1

    #para cada columna de la matriz de diseño
    for (i in seq_along(1:runs)) {
        group.all <- sdmdata[, i]
        group <- group.all[sdmdata$pa == 1]
        bg.grp <- group.all[sdmdata$pa == 0]
        backgr <- sdmdata[sdmdata$pa == 0, c("lon", "lat")]
        #para cada grupo
        for (g in setdiff(unique(group), 0)) {
            #excluding the zero allows for bootstrap. only 1 partition will run
            message(paste(species_name, algo, "run number", i, "partition number",
                      g, "\n"))
            pres_train <- occurrences[group != g, ]
            if (nrow(occurrences) == 1)
                pres_train <- occurrences[group == g, ]
            pres_test <- occurrences[group == g, ]
            backg_test <- backgr[bg.grp == g, ]
            sdmdata_train <- sdmdata[group.all != g,]#presences and absences
            envtrain <-  sdmdata_train[, names(predictors)] #ö ajeitar isto com grep.

            message("fitting models...")
            if (algo == "bioclim") mod <- dismo::bioclim(predictors, pres_train)
            if (algo == "maxent")  mod <- maxnet::maxnet(sdmdata_train$pa, envtrain)
            if (algo == "mahal")   mod <- dismo::mahal(predictors, pres_train)
            if (algo == "domain")  mod <- dismo::domain(predictors, pres_train)
            if (algo == "rf") {
              if (equalize == T) {
                #balanceando as ausencias
                abs <- nrow(sdmdata_train[sdmdata_train$pa == 0,])
                pres <- nrow(sdmdata_train[sdmdata_train$pa == 1,])
                prop <- pres:abs
                aus.eq <- sample(prop[-1], pres)
                envtrain.eq <- envtrain[c(1:pres, aus.eq),]
                sdmdata_train.eq <- sdmdata_train[c(1:pres, aus.eq),]
              } else {
                  envtrain.eq <- envtrain
                  sdmdata_train.eq <- sdmdata_train
              }
                #mod <- randomForest::randomForest(sdmdata_train.eq$pa ~ .,
                 #                               data = envtrain.eq,
                  #                              importance = T)
                mod <- randomForest::tuneRF(envtrain.eq,
                                            sdmdata_train.eq$pa,
                                            trace = F,
                                            plot = F,
                                            doBest = T,
                                            importance = F)
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
            if (algo == "brt") {
              if (equalize == T) {
              #balanceando as ausencias
              aus <- dim(sdmdata_train[sdmdata_train$pa == 0,])[1]
              pres <- dim(sdmdata_train[sdmdata_train$pa == 1,])[1]
              prop <- pres:aus
              aus.eq <- sample(prop[-1], pres)
              envtrain.eq <- envtrain[c(1:pres, aus.eq),]
              sdmdata_train.eq <- sdmdata_train[c(1:pres, aus.eq),]
              } else {
                  envtrain.eq <- envtrain
                  sdmdata_train.eq <- sdmdata_train
              }
                ind <- names(predictors)
                mod <- dismo::gbm.step(data = sdmdata_train.eq,
                                       gbm.x = ind,
                                       gbm.y = "pa",
                                       family = "bernoulli",
                                       tree.complexity = 5,
                                       learning.rate = 0.005,
                                       bag.fraction = 0.5,
                                       plot.main = FALSE)
                n.trees <- mod$n.trees
                }
            if (algo %in% c("centroid", "mindist")) {
                ec_cont <- predictors[[1]]
                ec_cont[!is.na(ec_cont)] <- 1
                predictors.st <- raster::scale(predictors)
                pres.vals <- raster::extract(predictors.st, pres_train)
                pred.vals <- raster::values(predictors.st)
                dist.vals <- vector(length = nrow(pred.vals))

                #esto es centroid
                if (algo == "centroid") {
                    cat(paste("Euclidean environmental distance to centroid","\n"))
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
                if (!nrow(occurrences) %in% c(1, 2)) {
                    #sólo corta por LTP si hay más de dos puntos...
                    LPTec <- dismo::threshold(eec, "no_omission")
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
                th_mod <- eval_mod@t[which.max(eval_mod@TPR + eval_mod@TNR)]
            } else if (algo == "brt") {
                eval_mod <- dismo::evaluate(pres_test, backg_test, mod,
                                            predictors, n.trees = n.trees)
                th_mod   <- eval_mod@t[which.max(eval_mod@TPR + eval_mod@TNR)]
                conf <- dismo::evaluate(pres_test, backg_test, mod, predictors,
                                        n.trees = n.trees, tr = th_mod)

                mod_cont <- dismo::predict(predictors, mod, n.trees = n.trees)
            } else if (algo %in% c("bioclim",
                                   "domain",
                                   "mahal")) {
                eval_mod <- dismo::evaluate(pres_test, backg_test, mod, predictors)
                th_mod   <- eval_mod@t[which.max(eval_mod@TPR + eval_mod@TNR)]
                conf <- dismo::evaluate(pres_test, backg_test, mod, predictors, tr = th_mod)
                mod_cont <- raster::predict(mod, predictors)
            } else if (algo %in% c("svm.k", "svm.e", "rf")) {
                eval_mod <- dismo::evaluate(pres_test, backg_test, mod, predictors)
                th_mod   <- eval_mod@t[which.max(eval_mod@TPR + eval_mod@TNR)]
                conf <- dismo::evaluate(pres_test, backg_test, mod, predictors, tr = th_mod)
                mod_cont <- raster::predict(predictors, mod)
            } else if (algo %in% "glm") {
                eval_mod <- dismo::evaluate(pres_test, backg_test, mod, predictors, type = "response")
                th_mod   <- eval_mod@t[which.max(eval_mod@TPR + eval_mod@TNR)]
                conf <- dismo::evaluate(pres_test, backg_test, mod, predictors,
                                        tr = th_mod)
                mod_cont <- raster::predict(predictors, mod, type = "response")
            }else if (algo %in% c( "maxent")) {
              eval_mod <- dismo::evaluate(pres_test, backg_test, mod, predictors, type = "logistic")
              th_mod   <- eval_mod@t[which.max(eval_mod@TPR + eval_mod@TNR)]
              conf <- dismo::evaluate(pres_test, backg_test, mod, predictors, tr = th_mod)
              mod_cont <- raster::predict(predictors, mod, type = "logistic")
              }


            message("evaluating the models...")
            th_table <- dismo::threshold(eval_mod)
            mod_TSS  <- max(eval_mod@TPR + eval_mod@TNR) - 1


            th_table$AUC <- eval_mod@auc
            th_table$TSS <- mod_TSS
            th_table$algoritmo <- algo
            th_table$run <- i
            th_table$partition <- g
            row.names(th_table) <- paste(species_name, i, g, algo)

            #confusion matrix
                #conf <- dismo::evaluate(pres_test, backg_test, mod, predictors, tr = th_mod)
            if (conf_mat == TRUE) {
              if(!algo %in% c("mindist", "centroid")){
                conf_res <- data.frame(presence_record = conf@confusion[,c("tp", "fp")],
                                       absence_record = conf@confusion[,c("fn", "tn")])
                rownames(conf_res) <- c("presence_predicted", "absence_predicted")
                write.csv(conf_res, file = paste0(partition.folder,
                                                  "/confusion_matrices_",
                                                  species_name, "_", i, "_", g,
                                                  "_", algo, ".csv"))
                }
              }

            th_table$presence <- eval_mod@np
            th_table$absence <- eval_mod@na
            th_table$correlation <- eval_mod@cor
            th_table$pvaluecor <- eval_mod@pcor


            if (conf_mat == TRUE) {
              if (!algo %in% c("mindist", "centroid")) {
                th_table$prevalence.value <- conf@prevalence
                th_table$PPP <- conf@PPP
                th_table$NPP <- conf@NPP
                th_table$sensitivity.value <- conf@TPR / (conf@TPR + conf@FPR)
                th_table$specificity.value <- conf@TNR / (conf@FNR + conf@TNR)
                th_table$comission <- conf@FNR / (conf@FNR + conf@TNR)
                th_table$omission <- conf@FPR / (conf@TPR + conf@FPR)
                th_table$accuracy <- (conf@TPR + conf@TNR) / (conf@TPR + conf@TNR + conf@FNR + conf@FPR)
                th_table$KAPPA.value <- conf@kappa
              }
            }

            #writing evaluation tables

            message("writing evaluation tables...")
            write.csv(th_table, file = paste0(partition.folder, "/evaluate_",
                                                species_name, "_", i, "_", g,
                                                "_", algo, ".csv"))


            if (class(mask) == "SpatialPolygonsDataFrame") {
                mod_cont <- crop_model(mod_cont, mask)
            }
            message("writing raster files...")
            raster::writeRaster(x = mod_cont,
                                filename = paste0(partition.folder, "/", algo,
                                                  "_cont_", species_name, "_",
                                                  i, "_", g, ".tif"),
                                overwrite = T)
            if (write_bin_cut == T) {
                message("writing binary and cut raster files...")
                mod_bin  <- mod_cont > th_mod
                mod_cut  <- mod_cont * mod_bin
                if (class(mask) == "SpatialPolygonsDataFrame") {
                    mod_bin <- crop_model(mod_bin, mask)
                    mod_cut <- crop_model(mod_cut, mask)
                }
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
            }


            if (write_png == T) {
                message("writing png files...")
                png(paste0(partition.folder, "/", algo, "_cont_", species_name,
                           "_", i, "_", g, ".png"))
                raster::plot(mod_cont,
                             main = paste(algo, "raw", "\n", "AUC =",
                                          round(eval_mod@auc, 2), "-", "TSS =",
                                          round(mod_TSS, 2)))
                dev.off()

                if (write_bin_cut == T){
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

            }

            if (project_model == T) {
                    pfiles <- list.dirs(proj_data_folder, recursive = F)
                    for (proje in pfiles) {
                        v <- strsplit(proje, "/")
                        name_proj <- v[[1]][length(v[[1]])]
                    projection.folder <- paste0(models_dir, "/", species_name,
                                                "/", name_proj,"/partitions")
                    if (file.exists(projection.folder) == FALSE)
                        dir.create(paste0(projection.folder), recursive = T, showWarnings = FALSE)
                    pred_proj <- raster::stack(list.files(proje, full.names = T))
                    pred_proj <- raster::subset(pred_proj, names(predictors))
                    #names(pred_proj) <- names(predictors)
                    message(name_proj)

                    message("projecting models")
                    if (algo == "brt") {
                        mod_proj_cont <- dismo::predict(pred_proj, mod, n.trees = n.trees, ...)
                    } else if (algo %in% c("bioclim",
                                           "domain",
                                           "maxent",
                                           "mahal")) {
                    mod_proj_cont <- dismo::predict(pred_proj, mod, ...)
                        } else if (algo %in% c("glm", "svm.k", "svm.e", "rf")) {
                            mod_proj_cont <- raster::predict(pred_proj, mod)
                        }

                    #mod_proj_bin <- ifelse(is.null(proj_cut), mod_proj_cont > th_mod, mod_proj_cont > proj_cut)
                    mod_proj_bin <- mod_proj_cont > th_mod
                    mod_proj_cut <- mod_proj_bin * mod_proj_cont
                    # Normaliza o modelo cut
                    #mod_proj_cut <- mod_proj_cut / maxValue(mod_proj_cut)
                    if (class(mask) == "SpatialPolygonsDataFrame") {
                        mod_proj     <- crop_model(mod_proj_cont, mask)
                        mod_proj_bin <- crop_model(mod_proj_bin, mask)
                        mod_proj_cut <- crop_model(mod_proj_cut, mask)
                    }
                    message("writing projected models raster")
                    raster::writeRaster(x = mod_proj_cont,
                                        filename = paste0(projection.folder,
                                                          "/", algo, "_cont_",
                                                          species_name, "_",
                                                          i, "_", g, ".tif"),
                                        overwrite = T)
                    if(write_bin_cut == T){
                      raster::writeRaster(x = mod_proj_bin,
                                          filename = paste0(projection.folder,
                                                            "/", algo, "_bin_",
                                                            species_name, "_",
                                                            i, "_", g, ".tif"),
                                          overwrite = T)
                      raster::writeRaster(x = mod_proj_cut,
                                          filename = paste0(projection.folder,
                                                            "/", algo, "_cut_",
                                                            species_name, "_",
                                                            i, "_", g, ".tif"),
                                          overwrite = T)
                    }

                    if (write_png == T) {
                        message("writing projected models .png")
                        png(paste0(projection.folder, "/", algo, "_cont_",
                                   species_name, "_", i, "_", g, ".png"))
                        raster::plot(mod_proj_cont,
                                     main = paste(algo, "proj_raw", "\n",
                                                  "AUC =",
                                                  round(eval_mod@auc, 2), "-",
                                                  "TSS =", round(mod_TSS, 2)))
                        dev.off()

                        if(write_bin_cut == T){
                          png(paste0(projection.folder, "/", algo, "_bin_",
                                     species_name, "_", i, "_", g, ".png"))
                          raster::plot(mod_proj_bin,
                                       main = paste(algo, "proj_bin", "\n",
                                                    "AUC =",
                                                    round(eval_mod@auc, 2), "-",
                                                    "TSS =", round(mod_TSS, 2)))
                          dev.off()
                          png(paste0(projection.folder, "/", algo, "_cut_",
                                     species_name, "_", i, "_", g, ".png"))
                          raster::plot(mod_proj_cut,
                                       main = paste(algo, "proj_cut", "\n",
                                                    "AUC =",
                                                    round(eval_mod@auc, 2), "-",
                                                    "TSS =", round(mod_TSS, 2)))
                          dev.off()
                        }

                    }
                    rm(pred_proj)

                }
            }
        }
    }

    return(th_table)
}
