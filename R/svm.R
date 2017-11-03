#' Faz modelagem de distribuição de espécies com algotimo SVM
#'
#' @inheritParams do_bioclim
#' @return Um data.frame com metadados da modelagem (TSS, AUC, algoritmo etc.)
#' @export
do_SVM <- function(sp,
		   coordinates,
		   partitions,
		   buffer = FALSE,
		   seed = 512,
		   predictors,
		   models.dir,
		   project.model,
		   projections,
		   mask,
		   n.back = 500) {
  cat(paste("SVM", "\n"))

  if (file.exists(paste0(models.dir)) == FALSE)
    dir.create(paste0(models.dir))
  if (file.exists(paste0(models.dir, "/", sp)) == FALSE) 
    dir.create(paste0(models.dir, "/", sp))
  if (project.model == T) {
    for (proj in projections) {
      if (file.exists(paste0(models.dir, "/", sp, "/", proj)) == FALSE) 
        dir.create(paste0(models.dir, "/", sp, "/", proj))
    }
  }

  # tabela de valores
  presvals <- raster::extract(predictors, coordinates)

  if (buffer %in% c("mean", "max", "median")) {
    backgr <- createBuffer(coord = coordinates, n.back = n.back, buffer.type = buffer, 
      occs = coordinates, sp = sp, seed = seed, predictors = predictors)
  } else {
    set.seed(seed + 2)
    backgr <- dismo::randomPoints(predictors, n.back)
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
  write.table(sdmdata, file = paste0(models.dir, "/", sp, "/sdmdata.txt"))


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

    svm <- kernlab::ksvm(sdmdata_train$pa ~ ., data = envtrain, cross = partitions)  ##svm deve ser com a variável resposta binária ou contínua, eu acho que binária
    esvm <- dismo::evaluate(envtest_pre, envtest_back, svm)
    # esvm <- evaluate(pres_test,backg_test,model = svm,x = predictors)
    svm_TSS <- max(esvm@TPR + esvm@TNR) - 1
    thresholdsvm <- esvm@t[which.max(esvm@TPR + esvm@TNR)]
    thsvm <- dismo::threshold(esvm)
    thsvm$AUC <- esvm@auc
    thsvm$TSS <- svm_TSS
    thsvm$algoritmo <- "svm"
    thsvm$partition <- i
    row.names(thsvm) <- paste(sp, i, "svm")
    svm_cont <- dismo::predict(predictors, svm, progress = "text")
    svm_bin <- svm_cont > thresholdsvm
    svm_cut <- svm_bin * svm_cont
    
    # TRANSFORMA 0 A 1
    svm_cont <- svm_cont/raster::maxValue(svm_cont)
    svm_cut <- svm_cut/raster::maxValue(svm_cut)

    write.table(thsvm, file = paste0(models.dir, "/", sp, "/evaluate", 
      sp, "_", i, "_svm.txt"))

    if (class(mask) == "SpatialPolygonsDataFrame") {
      svm_cont <- cropModel(svm_cont, mask)
      svm_bin <- cropModel(svm_bin, mask)
      svm_cut <- cropModel(svm_cut, mask)
    }
    
    raster::writeRaster(x = svm_cont, filename = paste0(models.dir, "/", sp, "/svm_cont_", 
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = svm_bin, filename = paste0(models.dir, "/", sp, "/svm_bin_", 
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = svm_cut, filename = paste0(models.dir, "/", sp, "/svm_cut_", 
      sp, "_", i, ".tif"), overwrite = T)
    
    
    if (project.model == T) {
      for (proj in projections) {
        data <- list.files(paste0("./env/", proj), pattern = proj)
        data2 <- stack(data)
        svm_proj <- predict(data2, svm, progress = "text")
        svm_proj_bin <- svm_proj > thresholdsvm
        svm_proj_cut <- svm_proj_bin * svm_proj
        
        # Normaliza o modelo cut svm_proj_cut <- svm_proj_cut/maxValue(svm_proj_cut)
        if (class(mask) == "SpatialPolygonsDataFrame") {
          source("./fct/cropModel.R")
          svm_proj <- cropModel(svm_proj, mask)
          svm_proj_bin <- cropModel(svm_proj_bin, mask)
          svm_proj_cut <- cropModel(svm_proj_cut, mask)
        }
        writeRaster(x = svm_proj, filename = paste0(models.dir, "/", sp, "/", 
          proj, "/svm_cont_", sp, "_", i, ".tif"), overwrite = T)
        writeRaster(x = svm_proj_bin, filename = paste0(models.dir, "/", sp, 
          "/", proj, "/svm_bin_", sp, "_", i, ".tif"), overwrite = T)
        writeRaster(x = svm_proj_cut, filename = paste0(models.dir, "/", sp, 
          "/", proj, "/svm_cut_", sp, "_", i, ".tif"), overwrite = T)
        rm(data2)
      }
    }
  }
  return(thsvm)
}
