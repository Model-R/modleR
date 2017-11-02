createBuffer <- function(coord, sp, occs, seed, n.back, 
  buffer.type = "mean", predictors) {
  # Transformando em spatial points
  sp::coordinates(coord) <- ~lon + lat
  
  if (buffer.type == "mean") 
    dist.buf <- mean(sp::spDists(x = coord, longlat = FALSE, segments = TRUE))
  if (buffer.type == "max") 
    dist.buf <- max(sp::spDists(x = coord, longlat = FALSE, segments = TRUE))
  if (buffer.type == "median")
    dist.buf <- median(sp::spDists(x = coord, longlat = FALSE, segments = TRUE))
  
  buffer <- raster::buffer(coord, width = dist.buf, dissolve = TRUE)
  
  # Transformando coords de novo em matriz para rodar resto script
  coord <- occs[occs$sp == sp, c("lon", "lat")]
  
  # Transformando em spatial polygon data frame
  buffer <- sp::SpatialPolygonsDataFrame(buffer, data = as.data.frame(buffer@plotOrder), 
    match.ID = FALSE)
  raster::crs(buffer) <- raster::crs(predictors)
  
  ######### TENHO CERTEZA DE QUE ISTO PODE FICAR MENOS PESADO Reference raster com mesmo
  ######### extent e resolution que predictors
  r_buffer <- raster::crop(predictors, buffer)
  r_buffer <- raster::mask(r_buffer, buffer)
  
  # r_buffer <- raster(ext=extent(predictors), resolution=res(predictors))
  
  # Rasterizando o buffer p/ geração dos ptos aleatorios r_buffer <-
  # rasterize(buffer, r_buffer, field=buffer@plotOrder) Limitando a mascara
  # ambiental r_buffer <- r_buffer*(predictors[[1]]!=0)
  
  
  # Gerando pontos aleatorios no buffer
  set.seed(seed + 2)
  backgr <- dismo::randomPoints(r_buffer, n.back)
  rm(buffer)
  gc()
  return(backgr)
}
