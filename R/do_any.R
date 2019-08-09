#' Model fitting and predicting of ecological niche models using several algorithms
#' 
#' This function reads the output from \code{\link{setup_sdmdata}} and 
#' computes ecological niche models for a species based on an algorithm specified by the user. It fits the model, calculate the predicted values and basic statistics for model evaluation. Besides from more commonly adopted metrics such as AUC and TSS, this package also calculates partial ROC (pROC) \insertCite{@for detais on model evaluation see @phillips_maximum_2006, @peterson_ecological_2011}{modleR}. Performs one algoritm at time, for runs with multiple algorithms see \code{\link{do_many}}. Given that there is "no silver bullets in correlative ecological niche modeling" \insertCite{qiao_no_2015}{modleR} the choice of which algorithm to run is on the user. See details for a description of all algorithms supported in this package.
#' 
#' @inheritParams setup_sdmdata
#' @inheritParams crop_model
#' @param algo The algorithm to be fitted \code{c("bioclim", "brt", "domain",
#'                                                 "glm", "maxent", "mahal", 
#'                                                 "svme", "svmk", "rf")}.
#' @param project_model Logical, whether to perform a projection.
#' @param proj_data_folder Path to directory with projections containing one or more folders with the projection datasets (e.g. "./env/proj/proj1"). Projection diretctory should only contain raster files corresponding to the environmental variables. If more than one projection, each projection should be at one directory (e.g. "./env/proj/proj1" and "./env/proj/proj2") and equivalent raster files at diferent subdirectories must have the same names (e.g. "./env/proj/proj1/layer1.asc" and "./env/proj/proj2/layer1.asc"). 
#' @param write_png Logical, whether png files will be written.
#' @param write_bin_cut Logical, whether binary and cut model files(.tif, .png) should be written.
#' @param threshold Character string indicating threshold (cut-off) to transform model predictions
#' to a binary score.
#' as in \code{\link[dismo]{threshold}}: "kappa", "spec_sens", "no_omission", "prevalence",
#' "equal_sens_spec", "sensitivity". Default value is "spec_sens".
#' @param conf_mat Logical, whether confusion tables should be written in the HD.
#' @param equalize Logical, whether the number of presences and absences should be
#' equalized in randomForest and brt.
#' @param proc_threshold Numeric, value from 0 to 100 that will be used as (E) for partialROC calculations in \code{\link[kuenm]{kuenm_proc}}. Default is \code{proc_threshold = 5}.
#' @param ... Other arguments from \code{\link[kuenm]{kuenm_proc}}.
#' @return Writes on disk model for each partition, a .csv file with evaluation statistics (TSS, AUC, etc).
#' @examples
#' # run setup_sdmdata first from one species in coordenadas data 
#' sp <- names(coordenadas)[1]
#' sp_coord <- coordenadas[[1]]
#' sp_setup <- setup_sdmdata(species_name=sp, occurrences=sp_coord, example_vars)
#' 
#' # run bioclim algorith for one species
#' do_any(species_name=sp,
#'        predictors=example_vars,
#'        algo = "bioclim")
#'         
#' @details See bellow for a description on the implementation of the algorithms supported in this package.
#' \describe{
#' \item{Bioclim}{
#' Specified by \code{algo="bioclim"} uses \code{\link[dismo]{bioclim}} function in \pkg{dismo} package \insertCite{hijmans_dismo_2017}{modleR}. Bioclim is the  climate-envelope-model implemented by Henry Nix \insertCite{nix_biogeographic_1986}{modleR}, the first species distribuition modelling package. It is based on climate interpolation methods and despite its limiations it is still used in ecological niche modeling, specially for exploration and teaching purposes \insertCite{@see also @booth_bioclim_2014}{modleR}. In this package it is implemented by the function \code{\link[dismo]{bioclim}}, evaluated and predicted using \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} also from \pkg{dismo} package.  
#' }
#' \item{Boosted Refression Trees (BRT)}{
#' Specified by \code{algo="brt"} uses \code{\link[dismo]{gbm.step}} function from \pkg{dismo} package. Runs the cross-validation procedure of \insertCite{hastie_elements_2001;textual}{modleR} \insertCite{@see also @elith_working_2009}{modleR}. It consists in a regression modeling technique combined with the boosting method, a method for combining many simple models. It is implemented by the function \code{\link[dismo]{gbm.step}} as a regression with the response variable set to bernoulli distribution. evaluated and predicted using \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} from \pkg{dismo} package.
#' }
#' \item{Domain}{
#' Specified by \code{algo="domain"} uses \code{\link[dismo]{domain}} function from \pkg{dismo} package. Computes point-to-point similarity based on Gower distance between environmental variables \insertCite{carpenter_domain_1993}{modleR}. \insertCite{hijmans_dismo_2017}{modleR} state that one should use it with caution because it does not perform well compared to other algorithms \insertCite{elith_novel_2006,hijmans_ability_2006}{modleR}. We add that it is a slow algorithm. In this package it is implemented by the function \code{\link[dismo]{domain}}, evaluated and predicted using \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} also from \pkg{dismo} package.
#' }
#' \item{Euclidean algorithms}{
#' To do or not to do.
#' }
#' \item{Generalized Linear Model (GLM)}{
#' Specified by \code{algo="glm"} runs a GLM with modeling presence and absences as a response variable following a binomial error distribution. It runs runs a step-wise model selection based on AIC both backward and forward considering all possible combinations of predictor variables in the rasterStack. In this package it is implemented using functions \code{glm} and \code{step} to fit a model and choose a model by AIC in a stepwise procedure. Model is evaluated and predicted using \code{\link[dismo]{evaluate}} function from \pkg{dismo} and \code{\link[raster]{predict}} function from \pkg{raster} package both with argument \code{type="response"} to return values in the scale of the response variable. 
#' }
#' \item{Mahalanobis}{
#' Specified by \code{algo="mahal"} uses \code{\link[dismo]{mahal}} function from \pkg{dismo} package. Corresponds to a distribution model based on Mahalanobis distance, a measure of the distance between a point P and a distribution D \insertCite{mahalanobis_generalized_1936}{modleR}. In this package it is implemented by the function \code{\link[dismo]{mahal}}, evaluated and predicted using \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} also from \pkg{dismo} package. 
#' }
#' \item{Maximum Entropy (Maxent)}{
#' Specified either by \code{algo="maxent"} or \code{algo="maxnet"} corresponding to implementation by \pkg{dismo} \insertCite{hijmans_dismo_2017}{modleR} and \pkg{maxnet} \insertCite{phillips_maxnet_2017}{modleR} packages respectivelly. Maxent is a machine learning method for modeling species distributions based in incomplete data allowing ENM with presence-only data \insertCite{phillips_maximum_2006}{modleR}. If \code{algo="maxent"} model is fitted by the function \code{\link[dismo]{maxent}}, evaluated and predicted using  \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} also in \pkg{dismo} package. If \code{algo="maxnet"} model is fitted by the function \code{\link[maxnet]{maxnet}} from \pkg{maxnet} package, evaluated using \code{\link[dismo]{evaluate}} from \pkg{dismo} package with argument \code{type="logistic"} and predicted using \code{\link[raster]{predict}} function from \pkg{raster} package. 
#' }
#' \item{Random Forest}{
#' Specified by \code{algo="rf"} uses \code{\link[randomForest]{tuneRF}} function from \pkg{ramdomForest} package \insertCite{liaw_classification_2002}{modleR}. Corresponds to machine learning regression based on decision trees. In this package uses \code{\link[randomForest]{tuneRF}} function with the optimal number of variables available for splitting at each tree node (i.e. mtry) found as set by parameter \code{doBest=TRUE}. Random Forest model is evaluated with \code{\link[dismo]{evaluate}} function from \pkg{dismo} and predicted with \code{\link[raster]{predict}} function from \pkg{raster} package.
#' } 
#' \item{Support Vector Machines (SVM)}{
#' Specified either by \code{algo="svme"} or \code{algo="svmk"} corresponding to implementation on \pkg{e1071} \insertCite{meyer_e1071_2017}{modleR} and \pkg{kernlab} \insertCite{karatzoglou_kernlab_2004}{modleR} packages respectivelly. SVM are supervised learning models that use learning algorithms for classification and regression analysis. In \pkg{e1071} package SVM is implemented through function \code{\link[e1071]{best.tune}} with method set to \code{"svm"} which uses RBF-kernel (radial basis function kernel) for classification. In \pkg{kernlab} package SVM is implemented through function \code{\link[kernlab]{ksvm}} also with RBF-kernel method (in this case the default method \code{"kbfdot"}). We expect  both implementations to differ only in performance. Both \code{svme} and \code{svmk} are evaluated with \code{\link[dismo]{evaluate}} function from dismo and predicted with \code{\link[raster]{predict}} function from \pkg{raster} package.
#' }
#' }
#' @references
#'     \insertAllCited{}
#' @author Andrea Sánchez-Tapia
#' @seealso \code{\link[dismo]{bioclim}}
#' @seealso \code{\link[dismo]{domain}}
#' @seealso \code{\link{do_any}}
#' @seealso \code{\link[dismo]{evaluate}}
#' @seealso \code{\link[dismo]{maxent}}
#' @seealso \code{\link[maxnet]{maxnet}}
#' @seealso \code{\link[dismo]{mahal}}
#' @seealso \code{\link[dismo]{predict}} in \pkg{dismo} package
#' @seealso \code{\link[raster]{predict}} in \pkg{raster} package
#' @import raster
#' @import grDevices
#' @importFrom utils write.csv
#' @importFrom maxnet maxnet
#' @importFrom stats complete.cases formula glm step dist
#' @importFrom Rdpack reprompt
#' @importFrom kuenm kuenm_proc
#' @export
do_any <- function(species_name,
                   predictors,
                   models_dir = "./models",
                   algo = c("bioclim"), #um so
                   project_model = FALSE,
                   proj_data_folder = "./data/proj",
                   mask = NULL,
                   write_png = FALSE,
                   write_bin_cut = FALSE,
                   threshold = "spec_sens",
                   conf_mat = TRUE,
                   equalize = TRUE,
                   proc_threshold = 0.5,
                   ...) {
  # replacing characters not welcome in species name
  # characters to avoid in file and dir names
  avoid_chars <- intToUtf8(c(91, 62, 33, 180, 60, 35, 63, 38, 47, 92, 46, 93))
  print_avoid <- intToUtf8(c(62, 33, 180, 60, 35, 63, 38, 47, 92, 46))
  if(grepl(avoid_chars, species_name)==TRUE){
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

    #para cada columna de la matriz de diseno
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

            message("fitting models...")
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

            message("projecting the models...")
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
            message("evaluating the models...")
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
            #which threshold? any value from function threhold() in dismo
            #th_mod <- eval_mod@t[which.max(eval_mod@TPR + eval_mod@TNR)]#tss?
            th_mod <- th_table[, threshold]
            th_table$threshold <- as.character(threshold)
            #confusion matrix
            if (algo == "brt") {
                conf <- dismo::evaluate(pres_test, backg_test, mod, predictors,
                                            n.trees = n.trees, tr = th_mod)
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

            message("writing evaluation tables...")
            write.csv(th_table, file = paste0(partition.folder, "/evaluate_",
                                              species_name, "_", i, "_", g,
                                              "_", algo, ".csv"))

                # apply mask (optional)
                if (class(mask) %in% c("SpatialPolygonsDataFrame",
                                       "SpatialPolygons")) {
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
                    message("writing png files...")
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
