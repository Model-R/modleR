#' Gera um modelo final por algortimo
#'
#' @param sp Um nome de espécie
#' @param select.partitions TRUE ou FALSE
#' @param algortimos Algortimos selecionados. Maxent, RF e SVM por padrão
#' @param threshold "spec_sens" por padrão
#' @param TSS.value 0.7 por padrão
#' @param models.dir Path do diretório onde estão os modelos
#' @param final.dir Em que subdiretório escrever os modelos finais
#' @return Uma string com a data e hora do término da execução
#' @import raster
#' @export
finalModel <- function(sp, select.partitions = T, algoritmos = c("maxent", "rf", 
  "svm"), threshold = c("spec_sens"), TSS.value = 0.7, models.dir = "./models", 
  final.dir = "presfinal") {
  if (file.exists(paste0(models.dir, "/", sp, "/", final.dir)) == FALSE) 
    dir.create(paste0(models.dir, "/", sp, "/", final.dir), recursive = TRUE)
  print(date())
  
  cat(paste(sp, "\n"))
  cat(paste("Reading the evaluation files", "\n"))
  evall <- list.files(path = paste0(models.dir, "/", sp), pattern = "evaluate", 
    full.names = T)
  lista <- list()
  for (i in 1:length(evall)) {
    lista[[i]] <- read.table(file = evall[i], header = T, row.names = 1)
  }
  stats <- data.table::rbindlist(lista)
  stats <- as.data.frame(stats)
  
  # Extracts only for the selected algorithm
  if (exists("algoritmos") == F) 
    algoritmos <- unique(stats$algoritmo)
  algoritmos <- as.factor(algoritmos)
  # algo <- algoritmos[1]
  todo <- raster::stack()
  for (algo in algoritmos) {
    cat(paste("Extracting data for", algo, "\n"))
    stats2 <- stats[stats$algoritmo == algo, ]
    
    
    part <- nrow(stats2)  #How many partitions were there
    
    cat(paste("Reading models from .tif files", "\n"))
    modelos.cont <- list.files(path = paste0(models.dir, "/", sp), full.names = T, 
      pattern = paste0(algo, "_cont_"))
    
    modelos.bin <- list.files(path = paste0(models.dir, "/", sp), full.names = T, 
      pattern = paste0(algo, "_bin_"))
    mod.cont <- raster::stack(modelos.cont)  #(0)
    
    mod.bin <- raster::stack(modelos.bin)  #(0)
    
    names(mod.cont) <- paste0(sp, algo, "Partition", 1:part)
    
    names(mod.bin) <- names(mod.cont)
    
    if (select.partitions == T) {
      cat("selecting partitions for", sp, algo, "\n")
      sel.index <- which(stats2[, "TSS"] >= TSS.value)
      cont.sel <- mod.cont[[sel.index]]  #(1)
      bin.sel <- mod.bin[[sel.index]]  #(5)
      
      
      th.mean <- mean(stats2[, names(stats2) == threshold][sel.index])
      
      if (length(sel.index) == 0) 
        cat(paste("No partition was selected for", sp, algo, "\n"))
      
      # en caso de que sea solo uno varios modelos son el mismo
      if (length(sel.index) == 1) {
        cat(paste(length(sel.index), "partition was selected for", sp, algo, 
          "\n"))
        
        # bin.sel #[3] bin.sel #[7]
        final <- raster::stack(bin.sel, bin.sel)
        names(final) <- c("Final.bin.mean3", "Final.mean.bin7")
      }
      
      # en caso de que sean más aplica el mapa
      
      if (length(sel.index) > 1) {
        cat(paste(length(sel.index), "partitions were selected for", sp, 
          "\n"))
        
        final.cont.mean <- mean(cont.sel)  #(2)
        final.bin.mean <- (final.cont.mean > th.mean)  #(3)
        final.sel.bin <- mean(bin.sel)  #(7)
        
        
        final <- raster::stack(final.bin.mean, final.sel.bin)
        names(final) <- c("Final.bin.mean3", "Final.mean.bin7")
      }
      
      
      if (exists("final")) {
        # plot(final) Escribe final
	raster::rasterOptions(setfileext = F)
        raster::writeRaster(x = final, filename = paste0(models.dir, "/", sp, "/", 
          final.dir, "/", names(final), sp, algo, ".tif"), bylayer = T, overwrite = T, 
          format = "GTiff")
        
#        for (i in 1:dim(final)[[3]]) {
#          png(filename = paste0(models.dir, "/", sp, "/", final.dir, "/", 
#          names(final)[i], sp, algo, ".png"))
#          plot(final[[i]], main = names(final)[i])
#          dev.off()
#        }
        todo <- raster::addLayer(todo, final)
      }
      cat("select", sp, algo, "DONE", "\n")
      
    }
  }
  
  print(paste("DONE", algo, "\n"))
  print(date())
  
  
  print(date())
}
