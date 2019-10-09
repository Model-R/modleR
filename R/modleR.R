#' modleR: A workflow to perform Ecological Niche Modeling based on dismo
#'
#' The modleR package wraps commonly used enm functions into a four-step
#' workflow: setup_sdmdata, do_(m)any, final_model and ensemble_model.
#'
#' @section modleR functions:
#' The three-step workflow implemented here to perform ecological niche modeling using functions
#' from \code{\link{dismo}} package.
#' 1. prepare data using \code{setup_sdmdata}
#' 2. run model(s) with \code{do_any} or \code{do_many} (for multiple
#' alforithms at time
#' 3. generate consensus model per algorithm with \code{final_model} and ensemble
#' consensus with \code{ensemble_model}
#'
#' @docType package
#' @name modleR
NULL
