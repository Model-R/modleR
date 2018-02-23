#' Fits ecological niche models using SVM from package e1071.
#'
#' @inheritParams do_bioclim
#' @return A data frame with the evaluation statistics (TSS, AUC, etc.)
#' @export
do_SVM2 <- function(sp,
                    coordinates,
                    partitions,
                    buffer,
                    seed,
                    predictors,
                    models.dir,
                    project.model,
                    projections,
                    mask,
                    write_png,
                    n.back) {
  cat(paste("SVM2", "\n"))

  if (file.exists(paste0(models.dir)) == FALSE)
       dir.create(paste0(models.dir))
    if (file.exists(paste0(models.dir, "/", sp)) == FALSE)
     dir.create(paste0(models.dir, "/", sp))
    partition.folder <- paste0(models.dir, "/", sp, "/present", "/partitions")
    if (file.exists(partition.folder) == FALSE)
        dir.create(partition.folder, recursive = T)

  # tabela de valores
  presvals <- raster::extract(predictors, coordinates)

  if (buffer %in% c("mean", "max", "median")) {
    backgr <- create_buffer(coord = coordinates,
                           n.back = n.back,
                           buffer.type = buffer,
                           seed = seed,
                           predictors = predictors)
  } else {
    set.seed(seed + 2)
    backgr <- dismo::randomPoints(mask = predictors,
                                  n = n.back,
                                  p = coordinates,
                                  excludep = T)
  }

  colnames(backgr) <- c("lon", "lat")

  # Extraindo dados ambientais dos bckgr
  backvals <- raster::extract(predictors, backgr)
  pa <- c(rep(1, nrow(presvals)), rep(0, nrow(backvals)))

  # Data partition
  if (nrow(coordinates) < 11)
    partitions <- nrow(coordinates)
  set.seed(seed)  #reproducibility
  group <- dismo::kfold(coordinates, partitions)
  set.seed(seed + 1)
  bg.grp <- dismo::kfold(backgr, partitions)
  group.all <- c(group, bg.grp)

  pres <- cbind(coordinates, presvals)
  back <- cbind(backgr, backvals)
  rbind_1 <- rbind(pres, back)
  sdmdata <- data.frame(cbind(group.all, pa, rbind_1))
  rm(rbind_1)
  rm(pres)
  rm(back)
  gc()
  write.table(sdmdata, file = paste0(partition.folder, "/sdmdata.txt"))


  ##### Hace los modelos
  for (i in unique(group)) {
    cat(paste(sp, "partition number", i, "\n"))
    pres_train <- coordinates[group != i, ]
    if (nrow(coordinates) == 1)
      pres_train <- coordinates[group == i, ]
    pres_test <- coordinates[group == i, ]

    backg_train <- backgr[bg.grp != i, ]  #not used?
    backg_test <- backgr[bg.grp == i, ]  #new

    sdmdata_train <- subset(sdmdata, group != i)  #new
    sdmdata_test <- subset(sdmdata, group == i)  #new

    envtrain <- subset(sdmdata_train, select = c(-group, -lon, -lat))  #new
    envtest <- subset(sdmdata_test, select = c(-group, -lon, -lat))
    envtest_pre <- subset(sdmdata_test, pa == 1, select = c(-group, -lon, -lat,
      -pa))  #new
    envtest_back <- subset(sdmdata_test, pa == 0, select = c(-group, -lon, -lat,
      -pa))  #new


    svm2 <- e1071::best.tune("svm", envtrain, sdmdata_train$pa, data = envtrain)  ##svm deve ser com a variável resposta binária ou contínua, eu acho que binária
    esvm2 <- dismo::evaluate(envtest_pre, envtest_back, svm2)
    # esvm <- evaluate(pres_test, backg_test, model = svm, x = predictors)
    svm2_TSS <- max(esvm2@TPR + esvm2@TNR) - 1
    thresholdsvm2 <- esvm2@t[which.max(esvm2@TPR + esvm2@TNR)]
    thsvm2 <- dismo::threshold(esvm2)
    thsvm2$AUC <- esvm2@auc
    thsvm2$TSS <- svm2_TSS
    thsvm2$algoritmo <- "svm2"
    thsvm2$partition <- i
    row.names(thsvm2) <- paste(sp, i, "svm2")
    svm2_cont <- dismo::predict(predictors, svm2, progress = "text")
    svm2_bin <- svm2_cont > thresholdsvm2
    svm2_cut <- svm2_bin * svm2_cont

    # TRANSFORMA 0 A 1
    svm2_cont <- svm2_cont / raster::maxValue(svm2_cont)
    svm2_cut <- svm2_cut / raster::maxValue(svm2_cut)

    write.table(thsvm2, file = paste0(partition.folder, "/evaluate",
      sp, "_", i, "_svm2.txt"))

    if (class(mask) == "SpatialPolygonsDataFrame") {
      svm2_cont <- crop_model(svm2_cont, mask)
      svm2_bin <- crop_model(svm2_bin, mask)
      svm2_cut <- crop_model(svm2_cut, mask)
    }

    raster::writeRaster(x = svm2_cont,
                        filename = paste0(partition.folder, "/svm2_cont_",
                                          sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = svm2_bin,
                        filename = paste0(partition.folder, "/svm2_bin_",
                                          sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = svm2_cut,
                        filename = paste0(partition.folder, "/svm2_cut_",
                                          sp, "_", i, ".tif"), overwrite = T)

        if (write_png == T) {
            png(filename = paste0(partition.folder,
                                  "/svm2", sp, "_", i, "%03d.png"))
            raster::plot(svm2_cont,
                         main = paste("SVM2 raw", "\n",
                                      "AUC =", round(esvm2@auc, 2), "-",
                                      "TSS =", round(svm2_TSS, 2)))
            raster::plot(svm2_bin,
                         main = paste("SVM2 P/A", "\n",
                                      "AUC =", round(esvm2@auc, 2), "-",
                                      "TSS =", round(svm2_TSS, 2)))
            raster::plot(svm2_cut,
                         main = paste("SVM2 cut", "\n",
                                      "AUC =", round(esvm2@auc, 2), "-",
                                      "TSS =", round(svm2_TSS, 2)))
            dev.off()
            }

    if (project.model == T) {
      for (proj in projections) {
      projection.folder <- paste0(models.dir, "/", sp, "/", proj)
            if (file.exists(projection.folder) == FALSE)
                dir.create(paste0(projection.folder), recursive = T)

        data <- list.files(paste0("./env/", proj), pattern = proj)
        data2 <- stack(data)
        svm2_proj <- predict(data2, svm2, progress = "text")
        svm2_proj_bin <- svm2_proj > thresholdsvm2
        svm2_proj_cut <- svm2_proj_bin * svm2_proj
        # Normaliza o modelo cut
        svm2_proj_cut <- svm2_proj_cut / maxValue(svm2_proj_cut)
        if (class(mask) == "SpatialPolygonsDataFrame") {
          svm2_proj <- crop_model(svm2_proj, mask)
          svm2_proj_bin <- crop_model(svm2_proj_bin, mask)
          svm2_proj_cut <- crop_model(svm2_proj_cut, mask)
        }
        writeRaster(x = svm2_proj,
                    filename = paste0(projection.folder, "/svm2_cont_", sp,
                                      "_", i, ".tif"), overwrite = T)
        writeRaster(x = svm2_proj_bin,
                    filename = paste0(projection.folder, "/svm2_bin_", sp,
                                      "_", i, ".tif"), overwrite = T)
        writeRaster(x = svm2_proj_cut,
                    filename = paste0(projection.folder, "/svm2_cut_", sp,
                                      "_", i, ".tif"), overwrite = T)
        rm(data2)
      }
    }
  }
  return(thsvm2)
}
