#' Ocurrence points for four species in the Brazilian Atlantic Forest.
#'
#' A list with four elements containing ccurrence points from four species: Abarema langsdorfii, 
#' Eugenia florida, Leandra carassana and Ouratea semiserrata. 
#' Each element contains three variables: sp species name separated by _, lat contains Latitude information and lon contains
#' Longitude information.
#'
#' @format A list of 4: 
#' \describe{
#'   \item{Abarema_langsdorffii}{104 observations of 3 variables}
#'   \item{Eugenia_florida}{341 observations of  3 variables}
#'   \item{Leandra_carassana}{82 observations of 3 variables}
#'   \item{Ouratea_semiserrata}{90 observations of 3 variables}
#' }
"coordenadas"

#' Mask of the Brazilian Atlantic Forest, based on IBGE
#'
#' @format A SpatialPolygonsDataFrame of the Brazilian Atlantic Forest
"mascara"

#' Predictor variables
#'
#' @format A RasterStack with 6 predictor variables issued from a PCA, at 10min resolution
"example_vars"
