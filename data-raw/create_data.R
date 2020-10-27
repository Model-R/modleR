# create data for the package with the newest crs information
library(raster)
library(usethis)

#get worldclim variables
raster::getData(name = "worldclim",
                download = T,
                var = "bio",
                res = 5,
                path = "./data-raw/")

# get brazil shapefile
raster::getData(name = "GADM",
                download = T,
                country = "BRA",
                level = 0,
                path = "./data-raw/")
BRA <- readRDS("./data-raw/gadm36_BRA_0_sp.rds")
#get biome shapefile
download.file("https://geoftp.ibge.gov.br/informacoes_ambientais/estudos_ambientais/biomas/vetores/Biomas_5000mil.zip", destfile = "./data-raw/Biomas_5000mil.zip")
unzip("./data-raw/Biomas_5000mil.zip", exdir = "./data-raw/")

Biomas <- rgdal::readOGR("data-raw/Biomas5000.shp",
                         encoding = "ISO-8859-1")

BAF <- Biomas[Biomas$COD_BIOMA == "MAT",]
rgdal::writeOGR(BAF, "./data-raw/BAF.shp",
                driver = "ESRI Shapefile", layer = "")
example_mask <- rgdal::readOGR("data-raw/BAF.shp")
use_data(example_mask, overwrite = TRUE)

#create example_vars
vars <- list.files("./data-raw/wc5", pattern = "bil$", full.names = TRUE)
vars_stack <- stack(vars)
plot(vars_stack[[1]])
vars <- mask(vars_stack, BRA)
vars <- crop(vars, BRA)
source("data-raw/pcavars.R")
example_vars <- pcavars(vars)
writeRaster(example_vars, filename = "./data-raw/example_vars.tif")
use_data(example_vars, overwrite = T)

