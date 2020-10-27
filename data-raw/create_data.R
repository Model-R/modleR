# create data for the package with the newest crs information
library(raster)
library(usethis)

#get worldclim variables
raster::getData(name = "worldclim",
                download = T,
                var = "bio",
                res = 10,
                path = "./data-raw/")

# get brazil shapefile
raster::getData(name = "GADM",
                download = T,
                country = "BRA",
                level = 0,
                path = "./data-raw/")
BRA <- readRDS("./data-raw/gadm36_BRA_0_sp.rds")

#get biome shapefile
download.file(
    "https://geoftp.ibge.gov.br/informacoes_ambientais/estudos_ambientais/biomas/vetores/Biomas_5000mil.zip",
              destfile = "./data-raw/Biomas_5000mil.zip")
unzip("./data-raw/Biomas_5000mil.zip", exdir = "./data-raw/")

Biomas <- rgdal::readOGR("data-raw/Biomas5000.shp",
                         encoding = "ISO-8859-1")

# select BAF BRazilian Atlantic Forest as an example mask
BAF <- Biomas[Biomas$COD_BIOMA == "MAT",]
rgdal::writeOGR(BAF, "./data-raw/BAF.shp",
                driver = "ESRI Shapefile", layer = "", encoding = "UTF-8")
example_mask <- rgdal::readOGR("data-raw/BAF.shp")
example_mask$NOM_BIOMA <- "Mata Atl\u00e2ntica"
use_data(example_mask, overwrite = TRUE, compress = "xz")

# create example_vars
vars <- list.files("./data-raw/wc10", pattern = "bil$", full.names = TRUE)
vars_stack <- stack(vars)
plot(vars_stack[[1]])
vars <- mask(vars_stack, BRA)
vars <- crop(vars, BRA)
source("data-raw/pcavars.R")
example_vars <- pcavars(vars, proportion = 0.9)
writeRaster(example_vars, filename = "./data-raw/example_vars.tif", overwrite = TRUE)
use_data(example_vars, overwrite = T, compress = "xz")
plot(example_vars)
