#' Ocurrence points for four species in the Brazilian Atlantic Forest.
#'
#' A dataframe with 618 occurrence points from four species
#'
#' @format A dataframe with 618 rows and 3 columns:
#' \describe{
#'   \item{sp}{Species name}
#'   \item{lat}{Latitude}
#'   \item{lon}{Longitude}
#' }
"coordenadas"

#' Mask of the Brazilian Atlantic Forest, based on IBGE
#'
#' @format A SpatialPolygonsDataFrame of the Brazilian Atlantic Forest
"mascara"

#' Predictor variables
#'
#' @format A RasterStack with 6 predictor variables issued from a PCA
"variaveis_preditoras"
