#' Ocurrence points for four species in the Brazilian Atlantic Forest
#'
#' A list with four elements containing occurrence points for four species:
#' \emph{Abarema langsdorfii}, \emph{Eugenia florida}, \emph{Leandra carassana}
#' and \emph{Ouratea semiserrata}.
#' Each element contains three variables: \code{sp} (species name separated by
#' "_"), \code{lat} and \code{lon} (Latitude and Longitude in decimal degrees)
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
#' @format A RasterStack with 6 predictor variables issued from a PCA of the 19
#' Worldclim 1.0 variables cut for South America, at 10min resolution
"example_vars"
