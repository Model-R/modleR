#' Model fitting and predicting of ecological niche models using several algorithms
#'
#' This function reads the output from \code{\link{setup_sdmdata}} and
#' computes ecological niche models for a species based on algorithm specified by the user. See details for
#' a description of all algorithms supported in this package.
#'
#' @inheritParams setup_sdmdata
#' @inheritParams crop_model
#' @param algo The algorithm to be fitted \code{c("bioclim", "maxent", "domain",
#'                                        "mahal", "glm", "svmk", "svme",
#'                                         "rf", "brt", "mindist", "centroid")}.
#' @param project_model Logical, whether to perform a projection.
#' @param proj_data_folder Path to projections -containing one or more
#'  folders with the projection datasets, ex. "./env/proj/proj1".
#' @param write_png Logical, whether png files will be written.
#' @param write_bin_cut Logical, whether binary and cut model files(.tif, .png) should be written.
#' @param threshold Character string indicating threshold (cut-off) to transform model predictions
#' to a binary score.
#' as in \code{\link[dismo]{threshold}}: "kappa", "spec_sens", "no_omission", "prevalence",
#' "equal_sens_spec", "sensitivity". Default value is "spec_sens".
#' @param conf_mat Logical, whether confusion tables should be written in the HD.
#' @param equalize Logical, whether the number of presences and absences should be
#' equalized in randomForest and brt.
#' @param proc_threshold Numeric, value from 0 to 100 that will be used as (E) for
#' partialROC calculations in \code{\link[kuenm]{kuenm_proc}} default = 5.
#' @param ... Other arguments from \code{\link[kuenm]{kuenm_proc}}
#' @return A data frame with the evaluation statistics (TSS, AUC, etc).
#' @details See below for a description on the implementation of the algorithms supported in this package.
#' \describe{
#' \item{Bioclim}{Specified by \code{algo="bioclim"} uses \code{\link[dismo]{bioclim}} function in dismo
#' package \insertCite{hijmans_dismo_2017}{modleR}. Bioclim is the climate-envelope-model implemented by Henry Nix
#' \insertCite{nix_biogeographic_1986}{modleR}, the first species
#' distribution modelling package. It is based on climate interpolation methods and despite its limitations
#' it is still used in ecological niche modeling, specially for exploration and teaching purposes
#' \insertCite{@see also @booth_bioclim_2014}{modleR}.
#' In this package it is implemented by the function \code{\link[dismo]{bioclim}}, evaluated and predicted
#' using \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} also from dismo package.
#' }
#' \item{Maximum Entropy (Maxent)}{Specified either by \code{algo="maxent"} or \code{algo="maxnet"}
#' corresponding to implementation by dismo \insertCite{hijmans_dismo_2017}{modleR} and maxnet
#' \insertCite{maxnet}{modleR} packages respectively. Maxent is a machine learning method for modeling
#' species distributions based on incomplete data allowing ENM with presence-only data
#' \insertCite{phillips_maximum_2006}{modleR}. If \code{algo="maxent"} model is fitted by the function
#' \code{\link[dismo]{maxent}}, evaluated and predicted using  \code{\link[dismo]{evaluate}} and
#' \code{\link[dismo]{predict}} also in dismo package. If \code{algo="maxnet"} model is fitted by the
#' function \code{\link[maxnet]{maxnet}} from maxnet package, evaluated using \code{\link[dismo]{evaluate}}
#' from dismo package with argument \code{type="logistic"} and predicted using \code{\link[raster]{predict}}
#' function from raster package.
#' }
#' \item{Mahalanobis}{Specified by \code{algo="mahal"} uses \code{\link[dismo]{mahal}} function from dismo
#' package. Corresponds to a distribution model based on Mahalanobis distance, a measure of the distance
#' between a point P and a distribution D \insertCite{mahalanobis_generalized_1936}{modleR}. In this package
#' it is implemented by the function \code{\link[dismo]{mahal}}, evaluated and predicted
#' using \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} also from dismo package.
#' }
#' \item{Domain}{Specified by \code{algo="domain"} uses \code{\link[dismo]{domain}} function from dismo
#' package. Computes  point-to-point similarity based on Gower distance between environmental variables
#' \insertCite{carpenter_domain_1993}{modleR}. \insertCite{hijmans_dismo_2017}{modleR} state that
#' one should use it with caution because it does not perform well compared to other algorithms
#' \insertCite{elith_novel_2006,hijmans_ability_2006}{modleR}.
#' We add that it is a slow algorithm.
#' In this package it is implemented by the function \code{\link[dismo]{domain}}, evaluated and predicted
#' using \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} also from dismo package.
#' }
#' \item{Support Vector Machines (SVM)}{
#' }
#' \item{GLM}{
#' }
#' \item{Random Forest}{
#' }
#' \item{Euclidean environmental distance}{
#' to do or not to do that is the question
#' }
#' \item{Boosted Regression Trees (BRT)}{
#' }
#' }
#' @references
#' \insertAllCited{}
#' @author Andrea Sánchez-Tapia & Sara Mortara
#' @seealso \code{\link[dismo]{bioclim}}
#' @seealso \code{\link[dismo]{maxent}}
#' @seealso \code{\link[maxnet]{maxnet}}
#' @seealso \code{\link[dismo]{domain}}
#' @seealso \code{\link[dismo]{mahal}}
#' @seealso \code{\link[dismo]{evaluate}}
#' @seealso \code{\link[dismo]{predict}}
#' @seealso \code{\link[raster]{predict}}
#' @import raster
#' @import grDevices
#' @importFrom utils write.csv
#' @importFrom maxnet maxnet
#' @importFrom stats complete.cases formula glm step dist
#' @importFrom Rdpack reprompt
#' @export
do_any <- function(species_name,
                   predictors,
                   models_dir = "./models",
                   algo = c("bioclim"),
                   project_model = FALSE,
                   proj_data_folder = "./data/proj",
                   mask = NULL,
                   write_png = FALSE,
                   write_bin_cut = FALSE,
                   threshold = "spec_sens",
                   conf_mat = TRUE,
                   equalize = TRUE,
                   proc_threshold = 0.5,
                   #probs,
                   ...) {
  # replacing characters not welcome in species name
  # characters to avoid in file and dir names
  avoid_chars <- intToUtf8(c(91, 62, 33, 180, 60, 35, 63, 38, 47, 92, 46, 93))
  print_avoid <- intToUtf8(c(62, 33, 180, 60, 35, 63, 38, 47, 92, 46))
  if (grepl(avoid_chars, species_name) == TRUE) {
    species_name <- gsub(avoid_chars, "", species_name)
    warning(cat(paste0('You entered a bad character (any in "',
                        print_avoid,
                        '") in the species name and we removed it for you')))
  }
    partition.folder <-
        paste(models_dir, species_name, "present", "partitions", sep = "/")
    if (file.exists(partition.folder) == FALSE)
        dir.create(partition.folder, recursive = T)
    setup.folder <-
        paste(models_dir, species_name, "present", "data_setup", sep = "/")

    # reads sdmdata from HD
    if (file.exists(paste(setup.folder, "sdmdata.txt", sep = "/"))) {
        sdmdata <- read.table(paste(setup.folder, "sdmdata.txt", sep = "/"))
    } else {
        stop("sdmdata.txt file not found, run setup_sdmdata() or check your folder settings")
        }

    message(paste(algo, "\n"))

    retained_predictors <-
        names(sdmdata)[(which(names(sdmdata) == "lat") + 1):ncol(sdmdata)]

    if (length(setdiff(names(predictors), retained_predictors)) > 0) {
        message(paste("Remember a variable selection was performed", "\n",
                      "retained variables:", paste(retained_predictors, collapse = "-"), "\n"))
    }
    predictors <- raster::subset(predictors, retained_predictors)


    ##### Hace los modelos
    runs <- which(names(sdmdata) == "pa") - 1

    #para cada coluna da matriz de desenho experimental
    for (i in seq_along(1:runs)) {
        group.all <- sdmdata[, i]
        group  <- group.all[sdmdata$pa == 1]
        bg.grp <- group.all[sdmdata$pa == 0]
        occurrences <- sdmdata[sdmdata$pa == 1, c("lon", "lat")]#isto recria ocorrencias
        backgr      <- sdmdata[sdmdata$pa == 0, c("lon", "lat")]
        #para cada grupo
        for (g in setdiff(unique(group), 0)) {
            #excluding the zero allows for bootstrap. only 1 partition will run
            message(paste(species_name, algo, "run number", i, "part. nb.",
                          g, "\n"))
            pres_train <- occurrences[group != g, ]
            if (nrow(occurrences) == 1) #only distance algorithms can be run
                pres_train <- occurrences[group == g, ]
            pres_test  <- occurrences[group == g, ]
            backg_test <- backgr[bg.grp == g, ]
            sdmdata_train <- sdmdata[group.all != g, ]#presences and absences
            envtrain <-  sdmdata_train[, names(predictors)] #presences and absences

            message("fitting models")
            if (algo == "bioclim") mod <- dismo::bioclim(predictors, pres_train)
            if (algo == "mahal")   mod <- dismo::mahal(predictors, pres_train)
            if (algo == "domain")  mod <- dismo::domain(predictors, pres_train)
            if (algo == "maxent")  mod <- dismo::maxent(envtrain, sdmdata_train$pa)
            if (algo == "maxnet")  mod <- maxnet::maxnet(sdmdata_train$pa, envtrain)
            if (algo == "glm") {
                null.model <- glm(sdmdata_train$pa ~ 1, data = envtrain,
                                  family = "binomial")
                full.model <- glm(sdmdata_train$pa ~ ., data = envtrain,
                                  family = "binomial")
                mod <- step(object = null.model, scope = formula(full.model),
                            direction = "both", trace = F)
            }
            if (algo == "svmk") {

                mod <- kernlab::ksvm(sdmdata_train$pa ~ ., data = envtrain)

            }
            if (algo == "svme") {
                sv <- 1
                while(!exists("mod")) {

                    mod <- e1071::best.tune("svm", envtrain, sdmdata_train$pa,
                                            data = envtrain)
                    sv <- sv + 1
                    message(paste("Trying svme", sv," times"))
                    if (sv == 10 & !exists("mod")) {
                        break
                        message("svme algorithm did not find a solution in 10 runs")
                    }
                }
            }
            if (algo == "rf" | algo == "brt") {
                if (equalize == T) {
                    #balanceando as ausencias
                    pres_train_n <- nrow(sdmdata_train[sdmdata_train$pa == 1, ])
                    abs_train_n  <- nrow(sdmdata_train[sdmdata_train$pa == 0, ])
                    prop <- pres_train_n:abs_train_n
                    aus.eq <- sample(prop[-1], pres_train_n)
                    envtrain.eq <- envtrain[c(1:pres_train_n, aus.eq), ]
                    sdmdata_train.eq <- sdmdata_train[c(1:pres_train_n, aus.eq), ]
                } else {
                    envtrain.eq <- envtrain
                    sdmdata_train.eq <- sdmdata_train
                }
                if (algo == "rf") {
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
                if (algo == "brt") {
                    mod <- dismo::gbm.step(data = sdmdata_train.eq,
                                           gbm.x = names(predictors),
                                           gbm.y = "pa",
                                           family = "bernoulli",
                                           tree.complexity = 1,
                                           learning.rate = 0.01,
                                           bag.fraction = 0.75,
                                           plot.main = FALSE,
                                           n.minobsinnode = 5)
                    n.trees <- mod$n.trees
                }
            }
            if (algo %in% c("centroid", "mindist")) {
                #esto es centroid
                if (algo == "centroid") {
                    cat(paste("Euclidean environmental distance to centroid",'\n'))
                    #calcula la media ambiental de los puntos de train
                    mod <- euclidean(predictors = predictors,
                                     occurrences = pres_train,
                                     #probs = probs,
                                     algo = "centroid",
                                     )
                    }
                if (algo == "mindist") {
                    cat(paste("Minimum Euclidean environmental distance",'\n'))
                    #calcula la media ambiental de los puntos de train
                    mod <- euclidean(predictors = predictors,
                                     occurrences = pres_train,
                                     #probs = probs,
                                     algo = "mindist"
                                     )
                    }
}
            message("projecting the models")
            if (exists("mod")) {
                if (algo == "brt") {
                    eval_mod <- dismo::evaluate(pres_test, backg_test, mod,
                                                predictors, n.trees = n.trees)
                    mod_cont <- dismo::predict(predictors, mod, n.trees = n.trees)
                }
                if (algo == "glm") {
                    eval_mod <- dismo::evaluate(pres_test, backg_test, mod,
                                                predictors, type = "response")
                    mod_cont <- raster::predict(predictors, mod, type = "response")
                }
                if (algo %in% c("bioclim",
                                "domain",
                                "maxent",
                                "mahal")) {
                    eval_mod <- dismo::evaluate(pres_test, backg_test, mod, predictors)
                    mod_cont <- dismo::predict(mod, predictors)
                }
                if (algo %in% c("svmk", "svme", "rf")) {
                    eval_mod <- dismo::evaluate(pres_test, backg_test, mod, predictors)
                    mod_cont <- raster::predict(predictors, mod)
                }
                if (algo == "maxnet") {
                    eval_mod <- dismo::evaluate(pres_test, backg_test, mod,
                                                predictors, type = "logistic")
                    mod_cont <- raster::predict(predictors, mod, type = "logistic")
                }
                if (algo %in% c("centroid", "mindist")) {
                mod_cont <- mod #não dá para projetar mas mod já é um raster continuo - a ideia é separar e que mod seja um objecto DistMod como aqueles gerados pelo dismo e que aqui dê para executar predict e evaluate.

                #este trecho de corte pelo LPT é velho
                #if (!nrow(occurrences) %in% c(1, 2)) {
                    #soh pode cortar se tiver mais de dos pontos
                 # p <- raster::extract(mod_cont, y = pres_test)
                  #a <- raster::extract(mod_cont, y = backg_test)
                  #eval_mod <- dismo::evaluate(p = p, a = a)
                  #cut pelo plt - opcional
                  #LPTec <- dismo::threshold(eval_mod, 'no_omission')
                  #[ö] este corte pelo LPT original poderia ser por uma porcentagem da distância ou por um valor - no caso nãõ seria necessário rodar o eval_mod e o lpt mas pegar o intervalo de distâncias e calcular a porcentagem desejada
                  #mod_cont[mod_cont < LPTec] <- 0
                  #mod_cont <-
                   #   (mod_cont - raster::minValue(mod_cont)) /
                    #  raster::maxValue(mod_cont - raster::minValue(mod_cont))
                    #ö is this the same as rescale?
                    #}
                p <- raster::extract(mod_cont, y = pres_test)
                a <- raster::extract(mod_cont, y = backg_test)
                eval_mod <- dismo::evaluate(p = p, a = a)#por enquanto só pode assim
                }

            message("evaluating the models")
            th_table <- dismo::threshold(eval_mod) #sensitivity 0.9
            #names(th_table) <- paste0(names(th_table), "_th")
            mod_TSS  <- max(eval_mod@TPR + eval_mod@TNR) - 1
            #PROC kuenm
            proc <- kuenm::kuenm_proc(occ.test = pres_test,
                                      model = mod_cont,
                                      threshold = proc_threshold,
                                      ...)

            #threshold-independent values
            th_table$AUC <- eval_mod@auc
            th_table$AUCratio <- eval_mod@auc/0.5
            th_table$pROC <- proc$pROC_summary[1]
            th_table$pval_pROC <- proc$pROC_summary[2]
            th_table$TSS <- mod_TSS
            th_table$algoritmo <- algo
            th_table$run <- i
            th_table$partition <- g
            th_table$presencenb <- eval_mod@np
            th_table$absencenb <- eval_mod@na
            th_table$correlation <- eval_mod@cor
            th_table$pvaluecor <- eval_mod@pcor
            row.names(th_table) <- species_name

            # threshold dependent values
            #which threshold? any value from function threshold() in dismo
            #th_mod <- eval_mod@t[which.max(eval_mod@TPR + eval_mod@TNR)]#tss?
            th_mod <- th_table[, threshold]
            th_table$threshold <- as.character(threshold)
            #confusion matrix
            if (algo == "brt") {
                conf <- dismo::evaluate(pres_test, backg_test, mod, predictors,
                                            n.trees = n.trees, tr = th_mod)
            } else if (algo %in% c("centroid", "mindist")) {
                conf <- dismo::evaluate(p = p, a = a, predictors,
                                        tr = th_mod)
            } else {
                conf <- dismo::evaluate(pres_test, backg_test, mod, predictors,
                                        tr = th_mod)
            }
            th_table$prevalence.value <- conf@prevalence
            th_table$PPP <- conf@PPP
            th_table$NPP <- conf@NPP
            th_table$sensitivity.value <- conf@TPR / (conf@TPR + conf@FPR)
            th_table$specificity.value <- conf@TNR / (conf@FNR + conf@TNR)
            th_table$comission <- conf@FNR / (conf@FNR + conf@TNR)
            th_table$omission <- conf@FPR / (conf@TPR + conf@FPR)
            th_table$accuracy <- (conf@TPR + conf@TNR) / (conf@TPR + conf@TNR + conf@FNR + conf@FPR)
            th_table$KAPPA.value <- conf@kappa

            #confusion matrix
            if (conf_mat == TRUE) {
                conf_res <-
                    data.frame(presence_record = conf@confusion[, c("tp", "fp")],
                               absence_record = conf@confusion[, c("fn", "tn")])
                rownames(conf_res) <- c("presence_predicted", "absence_predicted")
                write.csv(conf_res, file = paste0(partition.folder,
                                                  "/confusion_matrices_",
                                                  species_name, "_", i, "_", g,
                                                  "_", algo, ".csv"))
            }


            #writing evaluation tables

            message("writing evaluation tables")
            write.csv(th_table, file = paste0(partition.folder, "/evaluate_",
                                              species_name, "_", i, "_", g,
                                              "_", algo, ".csv"))

                # apply mask (optional)
                if (class(mask) %in% c("SpatialPolygonsDataFrame",
                                       "SpatialPolygons")) {
                    mod_cont <- crop_model(mod_cont, mask)
                }
                message("writing raster files")
                raster::writeRaster(x = mod_cont,
                                    filename = paste0(partition.folder, "/", algo,
                                                      "_cont_", species_name, "_",
                                                      i, "_", g, ".tif"),
                                    overwrite = T)
                if (write_bin_cut == T) {
                    message("writing binary and cut raster files")
                    mod_bin  <- mod_cont > th_mod
                    mod_cut  <- mod_cont * mod_bin
                    if (class(mask) == "SpatialPolygonsDataFrame") {
                        mod_bin <- crop_model(mod_bin, mask)
                        mod_cut <- crop_model(mod_cut, mask)
                    }
                    raster::writeRaster(x = mod_bin,
                                        filename = paste0(partition.folder, "/", algo,
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
                    message("writing png files")
                    png(paste0(partition.folder, "/", algo, "_cont_", species_name,
                               "_", i, "_", g, ".png"))
                    raster::plot(mod_cont,
                                 main = paste(algo, "raw", "\n", "AUC =",
                                              round(eval_mod@auc, 2), "-", "TSS =",
                                              round(mod_TSS, 2)))
                    dev.off()

                    if (write_bin_cut == T) {
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
                                                    "/", name_proj, "/partitions")
                        if (file.exists(projection.folder) == FALSE)
                            dir.create(paste0(projection.folder), recursive = T, showWarnings = FALSE)
                        pred_proj <- raster::stack(list.files(proje, full.names = T))
                        pred_proj <- raster::subset(pred_proj, names(predictors))
                        message(name_proj)

                        message("projecting models")
                        if (algo == "brt") {
                            mod_proj_cont <- dismo::predict(pred_proj, mod, n.trees = n.trees)
                        }
                        if (algo == "glm") {
                            mod_proj_cont <- raster::predict(pred_proj, mod, type = "response")
                        }
                        if (algo %in% c("bioclim",
                                        "domain",
                                        "maxent",
                                        "mahal")) {
                            mod_proj_cont <- dismo::predict(pred_proj, mod)
                        }
                        if (algo %in% c("svmk",
                                        "svme",
                                        "rf")) {
                            mod_proj_cont <- raster::predict(pred_proj, mod)
                        }

                        if (write_bin_cut == T) {
                            mod_proj_bin <- mod_proj_cont > th_mod
                            mod_proj_cut <- mod_proj_bin * mod_proj_cont
                            # Normaliza o modelo cut
                            #mod_proj_cut <- mod_proj_cut / maxValue(mod_proj_cut)
                        }
                        if (class(mask) == "SpatialPolygonsDataFrame") {
                            mod_proj_cont <- crop_model(mod_proj_cont, mask)
                            mod_proj_bin  <- crop_model(mod_proj_bin, mask)
                            mod_proj_cut  <- crop_model(mod_proj_cut, mask)
                        }
                        message("writing projected models raster")
                        raster::writeRaster(x = mod_proj_cont,
                                            filename = paste0(projection.folder,
                                                              "/", algo, "_cont_",
                                                              species_name, "_",
                                                              i, "_", g, ".tif"),
                                            overwrite = T)
                        if (write_bin_cut == T) {
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

                            if (write_bin_cut == T) {
                                png(paste0(projection.folder, "/", algo, "_bin_",
                                           species_name, "_", i, "_", g, ".png"))
                                raster::plot(mod_proj_bin,
                                             main = paste(algo, "proj_bin", "\n",
                                                          "AUC =",
                                                          round(eval_mod@auc, 2),
                                                          "-", "TSS =",
                                                          round(mod_TSS, 2)))
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
            } else message(paste(species_name, algo, "run number", i, "part. nb.",
                                 g, "could not be fit \n"))
        }

    }
    return(th_table)
    message("DONE!")
    print(date())
}
