#' Faz modelagem de distribuição de espécies com algoritmo Bioclim
#'
#' @param sp Um nome de espécie
#' @param coordinates Uma tabela com pontos de ocorrência
#' @param buffer Define se será usado buffer e de que tipo ("mean" ou "max")
#' @param seed Para reprodutibilidade
#' @param predictors Objeto do tipo RasterStack com variáveis preditoras
#' @param models.dir Path do diretório onde serão escritos os arquivos de saída
#' @param mask Objeto do tipo SpatialPolygonsDataFrame com máscara
#' @param write_png Se o png vai ser criado - defaults to F
#' @param n.back Número de pontos de background
#' @return Um data.frame com metadados da modelagem (TSS, AUC, algoritmo etc.)
#' @export
do_bioclim <- function(sp,
		       coordinates,
		       partitions,
		       buffer = FALSE,
		       seed = 512,
		       predictors,
		       models.dir,
		       project.model,
		       projections,
		       mask,
		       write_png = F,
		       n.back = 500) {
  cat(paste("Bioclim", "\n"))

  if (file.exists(paste0(models.dir)) == FALSE)
       dir.create(paste0(models.dir))
    if (file.exists(paste0(models.dir, "/", sp)) == FALSE)
     dir.create(paste0(models.dir, "/", sp))
    partition.folder <- paste0(models.dir,"/",sp,"/present","/partitions")
    if (file.exists(partition.folder) == FALSE)
        dir.create(partition.folder,recursive = T)
    
  # tabela de valores
  presvals <- raster::extract(predictors, coordinates)


  if (buffer %in% c("mean", "max", "median")) {
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
  write.table(sdmdata, file = paste0(partition.folder, "/sdmdata.txt"))


  ##### Hace los modelos
  for (i in unique(group)) {
    cat(paste(sp, "partition number", i, "\n"))
    pres_train <- coordinates[group != i, ]
    if (nrow(coordinates) == 1)
      pres_train <- coordinates[group == i, ]
    pres_test <- coordinates[group == i, ]

    backg_test <- backgr[bg.grp == i, ]  #new

    bc <- dismo::bioclim(predictors, pres_train)
    ebc <- dismo::evaluate(pres_test, backg_test, bc, predictors)
    thresholdbc <- ebc@t[which.max(ebc@TPR + ebc@TNR)]
    thbc <- dismo::threshold(ebc)
    bc_TSS <- max(ebc@TPR + ebc@TNR) - 1
    bc_cont <- dismo::predict(predictors, bc, progress = "text")
    bc_bin <- bc_cont > thresholdbc
    bc_cut <- bc_cont * bc_bin
    thbc$AUC <- ebc@auc
    thbc$TSS <- bc_TSS
    thbc$algoritmo <- "BioClim"
    thbc$partition <- i
    row.names(thbc) <- paste(sp, i, "BioClim")

    write.table(thbc, file = paste0(partition.folder, "/evaluate",
      sp, "_", i, "_bioclim.txt"))

    if (class(mask) == "SpatialPolygonsDataFrame") {
    bc_cont <- cropModel(bc_cont, mask)
    bc_bin <- cropModel(bc_bin, mask)
    bc_cut <- cropModel(bc_cut, mask)
    }
    raster::writeRaster(x = bc_cont, filename = paste0(partition.folder, "/BioClim_cont_",
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = bc_bin, filename = paste0(partition.folder, "/BioClim_bin_",
      sp, "_", i, ".tif"), overwrite = T)
    raster::writeRaster(x = bc_cut, filename = paste0(partition.folder, "/BioClim_cut_",
      sp, "_", i, ".tif"), overwrite = T)

  if (write_png == T) {
                  png(filename = paste0(partition.folder,"/Bioclim",sp,"_",i,"%03d.png"))
      raster::plot(bc_cont,main = paste("Bioclim raw","\n","AUC =", round(ebc@auc,2),'-',"TSS =",round(bc_TSS,2)))
      raster::plot(bc_bin,main = paste("Bioclim P/A","\n","AUC =", round(ebc@auc,2),'-',"TSS =",round(bc_TSS,2)))
      raster::plot(bc_cut,main = paste("Bioclim cut","\n","AUC =", round(ebc@auc,2),'-',"TSS =",round(bc_TSS,2)))
      dev.off()
      }

    if (project.model == T) {
      for (proj in projections) {
      projection.folder <- paste0(models.dir,"/",sp,"/",proj)
            if (file.exists(projection.folder) == FALSE)
                dir.create(paste0(projection.folder), recursive = T)

        data <- list.files(paste0("./env/", proj), pattern = proj)
        data2 <- stack(data)
        bc_proj <- predict(data2, do, progress = "text")
        bc_proj_bin <- bc_proj > thresholdbc
        bc_proj_cut <- bc_proj_bin * bc_proj
        # Normaliza o modelo cut
         bc_proj_cut <- bc_proj_cut/maxValue(bc_proj_cut)
        if (class(mask) == "SpatialPolygonsDataFrame") {
          source("./fct/cropModel.R")
          bc_proj <- cropModel(bc_proj, mask)
          bc_proj_bin <- cropModel(bc_proj_bin, mask)
          bc_proj_cut <- cropModel(bc_proj_cut, mask)
        }
        writeRaster(x = bc_proj, filename = paste0(projection.folder, "/BioClim_cont_", sp, "_", i, ".tif"), overwrite = T)
        writeRaster(x = bc_proj_bin, filename = paste0(projection.folder, "/BioClim_bin_", sp, "_", i, ".tif"), overwrite = T)
        writeRaster(x = bc_proj_cut, filename = paste0(projection.folder, "/BioClim_cut_", sp, "_", i, ".tif"), overwrite = T)
        rm(data2)
      }
    }
  }
  return(thbc)
}
