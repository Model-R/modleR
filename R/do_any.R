#' @rdname model_fit
#' @export
do_any <- function(species_name,
                   predictors,
                   models_dir = "./models",
                   algorithm = c("bioclim"),
                   project_model = FALSE,
                   proj_data_folder = "./data/proj",
                   mask = NULL,
                   write_rda = FALSE,
                   png_partitions = FALSE,
                   write_bin_cut = FALSE,
                   dismo_threshold = "spec_sens",
                   equalize = TRUE,
                   sensitivity = 0.9,
                   proc_threshold = 0.5,
                   ...) {
  # replacing characters not welcome in species name
  # characters to avoid in file and dir names
  avoid_chars <- intToUtf8(c(91, 62, 33, 180, 60, 35, 63, 38, 47, 92, 46, 93))
  print_avoid <- intToUtf8(c(62, 33, 180, 60, 35, 63, 38, 47, 92, 46))
  if (grepl(avoid_chars, species_name) == TRUE) {
    species_name <- gsub(avoid_chars, "", species_name)
    warning(cat(paste0('You entered a bad character (any in "', print_avoid, '")
                       in the species name and we removed it for you')))
  }
  partition_folder <-
    paste(models_dir, species_name, "present", "partitions", sep = "/")
    if (file.exists(partition_folder) == FALSE)
      dir.create(partition_folder, recursive = TRUE)
  setup_folder <-
    paste(models_dir, species_name, "present", "data_setup", sep = "/")

    #writes session info
    write_session_info(partition_folder)

    # reads sdmdata from HD
    if (file.exists(paste(setup_folder, "sdmdata.csv", sep = "/"))) {
        sdmdata <- read.csv(paste(setup_folder, "sdmdata.csv", sep = "/"))
    } else {
        stop("sdmdata.csv file not found, run setup_sdmdata() or check your
             folder settings")
        }

    message(paste(algorithm))

    retained_predictors <-
        names(sdmdata)[(which(names(sdmdata) == "lat") + 1):ncol(sdmdata)]

    if (length(setdiff(names(predictors), retained_predictors)) > 0) {
        message(paste("Remember a variable selection was performed", "\n",
                      "retained variables:", paste(retained_predictors,
                                                   collapse = "-"), "\n"))
    }
    predictors <- raster::subset(predictors, retained_predictors)


    ##### Hace los modelos
    runs <- which(names(sdmdata) == "pa") - 1

    #para cada coluna da matriz de desenho experimental
    for (i in seq_along(1:runs)) {
        group.all <- sdmdata[, i]
        group  <- group.all[sdmdata$pa == 1]
        bg.grp <- group.all[sdmdata$pa == 0]
        occurrences <- sdmdata[sdmdata$pa == 1, c("lon", "lat")]#recria occs
        backgr      <- sdmdata[sdmdata$pa == 0, c("lon", "lat")]
        #para cada grupo
        for (g in setdiff(unique(group), 0)) {
            #excluding the zero allows for bootstrap. only 1 partition will run
            message(paste(species_name, algorithm, "run number", i, "part. nb.",
                          g))
            pres_train <- occurrences[group != g, ]
            if (nrow(occurrences) == 1) #only distance algorithms can be run
                pres_train <- occurrences[group == g, ]
            pres_test  <- occurrences[group == g, ]
            backg_test <- backgr[bg.grp == g, ]
            sdmdata_train <- sdmdata[group.all != g, ]
            envtrain <-  sdmdata_train[, names(predictors)]

            message("fitting models")
            if (algorithm == "bioclim") mod <- dismo::bioclim(predictors, pres_train)
            if (algorithm == "mahal")   mod <- dismo::mahal(predictors, pres_train)
            if (algorithm == "domain")  mod <- dismo::domain(predictors, pres_train)
            if (algorithm == "maxent")
                mod <- dismo::maxent(envtrain, sdmdata_train$pa)
            if (algorithm == "maxnet")
                mod <- maxnet::maxnet(sdmdata_train$pa, envtrain)
            if (algorithm == "glm") {
                null.model <- glm(sdmdata_train$pa ~ 1, data = envtrain,
                                  family = "binomial")
                full.model <- glm(sdmdata_train$pa ~ ., data = envtrain,
                                  family = "binomial")
                mod <- step(object = null.model, scope = formula(full.model),
                            direction = "both", trace = FALSE)
            }
            if (algorithm == "svmk") {
                mod <- kernlab::ksvm(sdmdata_train$pa ~ ., data = envtrain)
            }
            if (algorithm == "svme") {
                sv <- 1
                while (!exists("mod")) {
                    mod <- e1071::best.tune("svm", envtrain,
                                            sdmdata_train$pa,
                                            data = envtrain)
                    sv <- sv + 1
                    message(paste("Trying svme", sv, "times"))
                    if (sv == 10 & !exists("mod")) {
                        break
                        message("svme algorithm did not converge to a solution in 10 runs")
                    }
                }
            }
            if (algorithm == "rf" | algorithm == "brt") {
                if (equalize == TRUE) {
                    #balanceando as ausencias
                    pres_train_n <- nrow(sdmdata_train[sdmdata_train$pa == 1, ])
                    abs_train_n  <- nrow(sdmdata_train[sdmdata_train$pa == 0, ])
                    prop <- pres_train_n:abs_train_n
                    aus.eq <- sample(prop[-1], pres_train_n)
                    envtrain.eq <- envtrain[c(1:pres_train_n, aus.eq), ]
                    sdmdata_train.eq <- sdmdata_train[c(1:pres_train_n, aus.eq),]
                } else {
                    envtrain.eq <- envtrain
                    sdmdata_train.eq <- sdmdata_train
                }
                if (algorithm == "rf") {
                    mod <- randomForest::tuneRF(envtrain.eq,
                                                sdmdata_train.eq$pa,
                                                trace = FALSE,
                                                plot = FALSE,
                                                doBest = TRUE,
                                                importance = FALSE)
                }
                if (algorithm == "brt") {
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

            message("projecting the models")
            if (exists("mod")) {
                if (algorithm == "brt") {
                    eval_mod <- dismo::evaluate(pres_test, backg_test, mod,
                                                predictors, n.trees = n.trees)
                    mod_cont <- dismo::predict(predictors, mod, n.trees = n.trees)
                }
                if (algorithm == "glm") {
                    eval_mod <- dismo::evaluate(pres_test, backg_test, mod,
                                                predictors, type = "response")
                    mod_cont <- raster::predict(predictors, mod, type = "response")
                }
                if (algorithm %in% c("bioclim",
                                "domain",
                                "maxent",
                                "mahal")) {
                    eval_mod <- dismo::evaluate(pres_test, backg_test, mod, predictors)
                    mod_cont <- dismo::predict(mod, predictors)
                }
                if (algorithm %in% c("svmk", "svme", "rf")) {
                    eval_mod <- dismo::evaluate(pres_test, backg_test, mod, predictors)
                    mod_cont <- raster::predict(predictors, mod)
                }
                if (algorithm == "maxnet") {
                    eval_mod <- dismo::evaluate(pres_test, backg_test, mod,
                                                predictors, type = "logistic")
                    mod_cont <- raster::predict(predictors, mod, type = "logistic")
                }
              message("evaluating the models")

              #evaluate as a complete data.frame
              eval_df <- data.frame(threshold = eval_mod@t,
                                    eval_mod@confusion,
                                    prevalence = eval_mod@prevalence,
                                    ODP = eval_mod@ODP,
                                    CCR = eval_mod@CCR,
                                    TPR = eval_mod@TPR,
                                    TNR = eval_mod@TNR,
                                    FPR = eval_mod@FPR,
                                    FNR = eval_mod@FNR,
                                    PPP = eval_mod@PPP,
                                    NPP = eval_mod@NPP,
                                    MCR = eval_mod@MCR,
                                    OR = eval_mod@OR,
                                    kappa = eval_mod@kappa,
                                    TSS = (eval_mod@TPR + eval_mod@TNR) - 1,
                                    FScore = 1/((1/eval_mod@TPR + 1/eval_mod@PPP)/2))
              eval_df$Jaccard  <- eval_df$tp / (eval_df$fn + eval_df$tp + eval_df$fp)
              #eval_df$Sorensen <- 2 * eval_df$tp / (eval_df$fn + 2 * eval_df$tp + eval_df$fp)#same as Fscore

            th_table <- dismo::threshold(eval_mod, sensitivity = sensitivity)
            #PROC kuenm
            proc <- kuenm::kuenm_proc(occ.test = pres_test,
                                      model = mod_cont,
                                      threshold = proc_threshold,
                                      ...)

            #threshold-independent values
            th_table$species_name <- species_name
            th_table$algorithm <- algorithm
            th_table$run <- i
            th_table$partition <- g
            th_table$presencenb  <- eval_mod@np
            th_table$absencenb   <- eval_mod@na
            th_table$correlation <- eval_mod@cor
            th_table$pvaluecor   <- eval_mod@pcor
            th_table$AUC         <- eval_mod@auc
            th_table$AUC_pval    <- ifelse(length(eval_mod@pauc) == 0, NA, eval_mod@pauc)
            th_table$AUCratio    <- eval_mod@auc / 0.5
            th_table$pROC        <- proc$pROC_summary[1]
            th_table$pROC_pval   <- proc$pROC_summary[2]
            th_table$TSSmax      <- max(eval_df$TSS)
            th_table$KAPPAmax    <- max(eval_df$kappa)

            # threshold dependent values
            #which threshold? any value from function threshold() in dismo
            th_table$dismo_threshold <- as.character(dismo_threshold)
            th_mod <- th_table[, dismo_threshold]

            th_table$prevalence.value <- eval_mod@prevalence[which(eval_mod@t == th_mod)]#a prevalencia desse threshold
            th_table$PPP <- eval_mod@PPP[which(eval_mod@t == th_mod)] #precision
            th_table$NPP <- eval_mod@NPP[which(eval_mod@t == th_mod)]
            th_table$TPR <- eval_mod@TPR[which(eval_mod@t == th_mod)] #sensitivity
            th_table$TNR <- eval_mod@TNR[which(eval_mod@t == th_mod)]#specificity
            th_table$FPR <- eval_mod@FPR[which(eval_mod@t == th_mod)]#comission
            th_table$FNR <- eval_mod@FNR[which(eval_mod@t == th_mod)]#omission
            th_table$CCR <- eval_mod@CCR[which(eval_mod@t == th_mod)] #accuracy
            th_table$Kappa <- eval_mod@kappa[which(eval_mod@t == th_mod)] #kappa
            th_table$F_score <- eval_df$FScore[which(eval_mod@t == th_mod)]
            th_table$Jaccard <- eval_df$Jaccard[which(eval_mod@t == th_mod)]
            #for tests: export th_table as a global object
            #th_table <<- th_table



            #writing evaluation tables
            message("writing evaluation tables")
            write.csv(eval_df, file = paste0(partition_folder, "/eval_mod_",
                                              species_name, "_", i, "_", g,
                                              "_", algorithm, ".csv"))
            write.csv(th_table, file = paste0(partition_folder, "/evaluate_",
                                              species_name, "_", i, "_", g,
                                              "_", algorithm, ".csv"))

                # apply mask (optional)
                if (class(mask) %in% c("SpatialPolygonsDataFrame",
                                       "SpatialPolygons")) {
                    mod_cont <- crop_model(mod_cont, mask)
                }
                message("writing raster files")
                raster::writeRaster(x = mod_cont,
                                    filename = paste0(partition_folder, "/", algorithm,
                                                      "_cont_", species_name, "_",
                                                      i, "_", g, ".tif"),
                                    overwrite = TRUE)
                if (write_bin_cut == TRUE) {
                    message("writing binary and cut raster files")
                    mod_bin  <- mod_cont > th_mod
                    mod_cut  <- mod_cont * mod_bin
                    if (class(mask) == "SpatialPolygonsDataFrame") {
                        mod_bin <- crop_model(mod_bin, mask)
                        mod_cut <- crop_model(mod_cut, mask)
                    }
                    raster::writeRaster(x = mod_bin,
                                        filename = paste0(partition_folder, "/", algorithm,
                                                          "_bin_", species_name, "_",
                                                          i, "_", g, ".tif"),
                                        overwrite = TRUE)
                    raster::writeRaster(x = mod_cut,
                                        filename = paste0(partition_folder, "/", algorithm,
                                                          "_cut_", species_name, "_",
                                                          i, "_", g, ".tif"),
                                        overwrite = TRUE)
                }

                if (write_rda == TRUE) {
                    message("writing .rda objects")
                    save(mod, file = paste0(partition_folder, "/", algorithm,
                                            "_model_", species_name, "_", i,
                                            "_", g, ".rda"))
                }

                if (png_partitions == TRUE) {
                    message("writing png files")
                    png(paste0(partition_folder, "/", algorithm, "_cont_", species_name,
                               "_", i, "_", g, ".png"))
                    raster::plot(mod_cont,
                                 main = paste(algorithm, "raw", "\n", "AUC =",
                                              round(eval_mod@auc, 2), "-", "TSS =",
                                              round(th_table$TSSmax, 2)))
                    dev.off()

                    if (write_bin_cut == TRUE) {
                        png(paste0(partition_folder, "/", algorithm, "_bin_", species_name,
                                   "_", i, "_", g, ".png"))
                        raster::plot(mod_bin,
                                     main = paste(algorithm, "bin", "\n", "AUC =",
                                                  round(eval_mod@auc, 2), "-", "TSS =",
                                                  round(th_table$TSSmax, 2)))
                        dev.off()
                        png(paste0(partition_folder, "/", algorithm, "_cut_", species_name,
                                   "_", i, "_", g, ".png"))
                        raster::plot(mod_cut,
                                     main = paste(algorithm, "cut", "\n", "AUC =",
                                                  round(eval_mod@auc, 2), "-", "TSS =",
                                                  round(th_table$TSSmax, 2)))
                        dev.off()
                    }

                }

                if (project_model == TRUE) {
                    pfiles <- list.dirs(proj_data_folder, recursive = FALSE)
                    for (proje in pfiles) {
                        v <- strsplit(proje, "/")
                        name_proj <- v[[1]][length(v[[1]])]
                        projection_folder <- paste0(models_dir, "/", species_name,
                                                    "/", name_proj, "/partitions")
                        if (file.exists(projection_folder) == FALSE)
                            dir.create(paste0(projection_folder),
                                       recursive = TRUE, showWarnings = FALSE)
                        pred_proj <- raster::stack(list.files(proje,
                                                              full.names = TRUE))
                        pred_proj <- raster::subset(pred_proj, names(predictors))
                        message(name_proj)

                        message("projecting models")
                        if (algorithm == "brt") {
                            mod_proj_cont <- dismo::predict(pred_proj,
                                                            mod,
                                                            n.trees = n.trees)
                        }
                        if (algorithm == "glm") {
                            mod_proj_cont <- raster::predict(pred_proj, mod,
                                                             type = "response")
                        }
                        if (algorithm %in% c("bioclim",
                                        "domain",
                                        "maxent",
                                        "mahal")) {
                            mod_proj_cont <- dismo::predict(pred_proj, mod)
                        }
                        if (algorithm %in% c("svmk",
                                        "svme",
                                        "rf")) {
                            mod_proj_cont <- raster::predict(pred_proj, mod)
                        }

                        if (write_bin_cut == TRUE) {
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
                                            filename = paste0(projection_folder,
                                                              "/", algorithm, "_cont_",
                                                              species_name, "_",
                                                              i, "_", g, ".tif"),
                                            overwrite = TRUE)

                        if (write_bin_cut == TRUE) {
                            raster::writeRaster(x = mod_proj_bin,
                                                filename = paste0(projection_folder,
                                                                  "/", algorithm, "_bin_",
                                                                  species_name, "_",
                                                                  i, "_", g, ".tif"),
                                                overwrite = TRUE)
                            raster::writeRaster(x = mod_proj_cut,
                                                filename = paste0(projection_folder,
                                                                  "/", algorithm, "_cut_",
                                                                  species_name, "_",
                                                                  i, "_", g, ".tif"),
                                                overwrite = TRUE)
                        }


                      # creating and writing do_any metadata
            metadata <- data.frame(
              species_name = as.character(species_name),
              algorithm = algorithm,
              project_model_dir = ifelse(project_model, proj_data_folder, NA),
              projections = ifelse(project_model, paste(pfiles, collapse = "-"), NA),
              mask = ifelse(is.null(mask), NA, mask),
              dismo_threshold = dismo_threshold,
              equalized_occ = ifelse(equalize, "yes", "no"),
              n_equalized_occ = ifelse(algorithm %in% c("rf", "brt")
                                       & equalize == TRUE,
                                       length(aus.eq), NA),
              proc_threshold = proc_threshold
            )
            message("writing metadata")
            write.csv(metadata, file = paste0(partition_folder, "/metadata_",
                                              algorithm, ".csv"))

                        if (png_partitions == TRUE) {
                            message("writing projected models .png")
                            png(paste0(projection_folder, "/", algorithm, "_cont_",
                                       species_name, "_", i, "_", g, ".png"))
                            raster::plot(mod_proj_cont,
                                         main = paste(algorithm, "proj_raw", "\n",
                                                      "AUC =",
                                                      round(eval_mod@auc, 2), "-",
                                                      "TSS =", round(th_table$TSSmax, 2)))
                            dev.off()

                            if (write_bin_cut == TRUE) {
                                png(paste0(projection_folder, "/", algorithm, "_bin_",
                                           species_name, "_", i, "_", g, ".png"))
                                raster::plot(mod_proj_bin,
                                             main = paste(algorithm, "proj_bin", "\n",
                                                          "AUC =",
                                                          round(eval_mod@auc, 2),
                                                          "-", "TSS =",
                                                          round(th_table$TSSmax, 2)))
                                dev.off()
                                png(paste0(projection_folder, "/", algorithm, "_cut_",
                                           species_name, "_", i, "_", g, ".png"))
                                raster::plot(mod_proj_cut,
                                             main = paste(algorithm, "proj_cut", "\n",
                                                          "AUC =",
                                                          round(eval_mod@auc, 2), "-",
                                                          "TSS =", round(th_table$TSSmax, 2)))
                                dev.off()
                            }

                        }
                        rm(pred_proj)
                    }
                }
            } else message(paste(species_name, algorithm, "run number", i, "part. nb.",
                                 g, "could not be fit"))
        }

    }
    return(th_table)
    message("DONE!")
    print(date())
}
