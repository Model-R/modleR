#' modleR: A workflow to perform ecological niche modeling
#'
#' The modleR package wraps commonly used enm functions into a four-step
#' workflow: \code{setup_sdmdata}, \code{do_any} and \code{do_many},
#' \code{final_model} and \code{ensemble_model}.
#'
#' @section modleR functions:
#' The four-step workflow implemented here to perform ecological niche modeling using functions from \code{\link{dismo}} package.
#'
#' 1. prepare data  using \code{setup_sdmdata()}
#'
#' 2. run model(s) with \code{do_any()} or \code{do_many()} (for multiple
#' algorithms at a time
#'
#' 3. generate a final model per algorithm with \code{final_model()} and
#'
#' 4. create ensemble models (algorithm consensus) with \code{ensemble_model()}
#'
#' @docType package
#' @name modleR
NULL
