#' Model fitting, predicting and evaluating of ecological niche models using one or several algorithms
#'
#' \code{\link{do_any}} reads the output from \code{\link{setup_sdmdata}} and computes ecological niche models for a species based on an algorithm specified by the user. It fits the model, calculates the predicted values and basic statistics for model evaluation. in addition to commonly adopted metrics such as AUC and TSS, this package also calculates partial ROC (pROC) (for details on model evaluation see \insertCite{phillips_maximum_2006;textual}{modleR} and \insertCite{peterson_ecological_2011;textual}{modleR}.
#' \code{\link{do_any}} performs one algorithm at a time. \code{do_many} runs internally \code{\link{do_any}} and can be used to run multiple algorithms at a time. See \strong{Details} in \code{\link{do_any}} for a description of how each algorithm is implemented. Given that there are "\emph{no silver bullets in correlative ecological niche modeling}" \insertCite{qiao_no_2015}{modleR} the choice of which algorithm to run is on the user. See details for a description of all algorithms supported in this package.
#'
#' @inheritParams setup_sdmdata
#' @inheritParams crop_model
#' @param algo The algorithm to be fitted \code{c("bioclim", "brt", "domain",
#'                                                 "glm", "maxent", "mahal",
#'                                                 "svme", "svmk", "rf")}.
#' @param bioclim Execute bioclim algorithm from the \pkg{dismo} implementation
#' with \code{\link[dismo]{bioclim}} function.
#' @param brt Execute Boosted Regression Trees with
#' \code{\link[dismo]{gbm.step}} from \pkg{dismo}.
#' @param domain Execute domain from the \pkg{dismo} implementation with
#' \code{\link[dismo]{domain}} function.
#' @param glm Execute GLM as suggested by the \pkg{dismo} documentation with
#' \code{\link{glm}} and \code{\link{step}}.
#' @param mahal Execute Mahalanobis distance from the \pkg{dismo} implementation
#'  with \code{\link[dismo]{mahal}}.
#' @param maxent Execute Maxent algorithm from the \pkg{dismo} implementation with
#'  \code{\link[dismo]{maxent}} function.
#' @param maxnet Execute Maxent algorithm from the \pkg{maxnet} implementation with
#'  \code{\link[maxnet]{maxnet}} function.
#' @param rf Execute Random forests algorithm from \pkg{randomForest} package with
#'  function\code{\link[randomForest]{tuneRF}} as suggested by the \pkg{dismo} documentation.
#' @param svme Execute Support Vector Machines (SVM) algorithm from \pkg{e1071}
#'  package with \code{\link[e1071]{best.tune}} function.
#' @param svmk Execute Support Vector Machines (SVM) algorithm from
#' \pkg{kernlab} package with \code{\link[kernlab]{ksvm}} function.
#' @param ... Any parameter from \link{do_any}
#' @param project_model Logical, whether to perform a projection.
#' @param proj_data_folder Path to directory with projections containing one or more folders with the projection datasets (e.g. "./env/proj/proj1"). Projection diretctory should only contain raster files corresponding to the environmental variables. If more than one projection, each projection should be at one directory (e.g. "./env/proj/proj1" and "./env/proj/proj2") and equivalent raster files at diferent subdirectories must have the same names (e.g. "./env/proj/proj1/layer1.asc" and "./env/proj/proj2/layer1.asc").
#' @param write_png Logical, whether png files will be written.
#' @param write_bin_cut Logical, whether binary and cut model files(.tif, .png) should be written.
#' @param threshold Character string indicating threshold (cut-off) to transform model predictions
#' to a binary score.
#' as in \code{\link[dismo]{threshold}}: "kappa", "spec_sens", "no_omission", "prevalence",
#' "equal_sens_spec", "sensitivity". Default value is "spec_sens".
#' @param conf_mat Logical, whether confusion tables should be written in the HD.
#' @param equalize Logical, whether the number of presences and absences should be
#' equalized in randomForest and brt.
#' @param proc_threshold Numeric, value from 0 to 100 that will be used as (E) for partialROC calculations in \code{\link[kuenm]{kuenm_proc}}. Default is \code{proc_threshold = 5}.
#' @param ... Other arguments from \code{\link[kuenm]{kuenm_proc}}.
#'
#' @name model_fit
#'
NULL
