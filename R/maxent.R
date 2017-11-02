#' Faz modelagem de distribuição de espécies com algotimo Maxent
#'
#' @inheritParams do_bioclim
#' @return Um data.frame com metadados da modelagem (TSS, AUC, algoritmo etc.)
#' @export
do_maxent <- function(sp,
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
  cat(paste("Maxent", "\n"))

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
  
  if (buffer %in% c("mean", "max")) {
    backgr <- createBuffer(coord = coordinates, n.back = n.back, buffer.type = buffer,
                           sp = sp, seed = seed, predictors = predictors)
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
  write.table(sdmdata, file = paste0(models.dir, "/", sp, "/sdmdata.txt"))
  
#  if (! file.exists(file = paste0(models.dir, "/", sp, "/evaluate", sp, "_", i, ".txt"))) {
#    write.table(data.frame(kappa = numeric(), spec_sens = numeric(), no_omission = numeric(), prevalence = numeric(), 
#			         equal_sens_spec = numeric(), sensitivity = numeric(), AUC = numeric(), TSS = numeric(), algoritmo = character(), 
#				 partition = numeric()), file = paste0(models.dir, "/", sp, "/evaluate", sp, "_", i, ".txt"))
#  }
  
  ##### Hace los modelos
  for (i in unique(group)) {
    cat(paste(sp, "partition number", i, "\n"))
    pres_train <- coordinates[group != i, ]
    if (nrow(coordinates) == 1) 
      pres_train <- coordinates[group == i, ]
    pres_test <- coordinates[group == i, ]
    
    backg_test <- backgr[bg.grp == i, ]  #new
    
    mx <- dismo::maxent(predictors, pres_train)
    emx <- dismo::evaluate(pres_test, backg_test, mx, predictors)
    thresholdmx <- emx@t[which.max(emx@TPR + emx@TNR)]
    thmx <- dismo::threshold(emx)
    mx_TSS <- max(emx@TPR + emx@TNR) - 1
    mx_cont <- dismo::predict(predictors, mx, progress = "text")
    mx_bin <- mx_cont > thresholdmx
    mx_cut <- mx_cont * mx_bin
    thmx$AUC <- emx@auc
    thmx$TSS <- mx_TSS
    thmx$algoritmo <- "maxent"
    thmx$partition <- i
    row.names(thmx) <- paste(sp, i, "maxent")

    write.table(thmx, file = paste0(models.dir, "/", sp, "/evaluate", 
      sp, "_", i, "_maxent.txt"))
    
    if (class(mask) == "SpatialPolygonsDataFrame") {
      mx_cont <- cropModel(mx_cont, mask)
      mx_bin <- cropModel(mx_bin, mask)
      mx_cut <- cropModel(mx_cut, mask)
    }
    raster::writeRaster(x = mx_cont, filename = paste0(models.dir, "/", sp, "/maxent_cont_", 
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = mx_bin, filename = paste0(models.dir, "/", sp, "/maxent_bin_", 
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = mx_cut, filename = paste0(models.dir, "/", sp, "/maxent_cut_", 
      sp, "_", i, ".tif"), overwrite = T)
    
  
    if (project.model == T) {
      for (proj in projections) {
        data <- list.files(paste0("./env/", proj), pattern = proj)
        data2 <- stack(data)
        mx_proj <- predict(data2, mx, progress = "text")
        mx_proj_bin <- mx_proj > thresholdmx
        mx_proj_cut <- mx_proj_bin * mx_proj
        # Normaliza o modelo cut do_proj_cut <- do_proj_cut/maxValue(do_proj_cut)
        if (class(mask) == "SpatialPolygonsDataFrame") {
          source("../../fct/cropModel.R")
          mx_proj <- cropModel(mx_proj, mask)
          mx_proj_bin <- cropModel(mx_proj_bin, mask)
          mx_proj_cut <- cropModel(mx_proj_cut, mask)
        }
        writeRaster(x = mx_proj, filename = paste0(models.dir, "/", sp, "/", 
          proj, "/maxent_cont_", sp, "_", i, ".tif"), overwrite = T)
        writeRaster(x = mx_proj_bin, filename = paste0(models.dir, "/", sp, "/", 
          proj, "/maxent_bin_", sp, "_", i, ".tif"), overwrite = T)
        writeRaster(x = mx_proj_cut, filename = paste0(models.dir, "/", sp, "/", 
          proj, "/maxent_cut_", sp, "_", i, ".tif"), overwrite = T)
        rm(data2)
      }
    }
  }
  return(thmx)
}
#    eval_df <- data.frame(kappa = 1, spec_sens = 1, no_omission = 1, prevalence = 1, 
#      equal_sens_spec = 1, sensitivity = 1, AUC = 1, TSS = 1, algoritmo = "foo", 
#      partition = 1)
#      eval_df <- rbind(eval_df, thdo)
