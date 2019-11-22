#' Ecological niche models fit, prediction, projection and evaluation using one or several algorithms
#'
#' \code{do_any} reads the output from \code{\link{setup_sdmdata}} and
#' computes ecological niche models for a species based on an algorithm
#' specified by the user. It fits the model, predicts it into the current
#' environmental layers and calculates basic statistics for model evaluation. In
#' addition to commonly adopted metrics such as AUC and TSS, this package also
#' calculates partial ROC
#' \insertCite{peterson_rethinking_2008,cobos_kuenm_2019}{modleR}. For details
#' on model evaluation see
#'  \insertCite{phillips_maximum_2006;textual}{modleR} and
#'  \insertCite{peterson_ecological_2011;textual}{modleR}. \code{do_any}
#'  performs one algorithm at a time. \code{do_many} runs internally
#'  \code{do_any} and can be used to run multiple algorithms at a time.
#'  Given that there are "\emph{no silver bullets in correlative ecological
#'  niche modeling}" \insertCite{qiao_no_2015}{modleR} the choice of which
#'  algorithm to run is on the user. See \strong{Details} for a description of
#'  how each algorithm supported in this package is implemented.
#'
#' @inheritParams setup_sdmdata
#' @return Returns a data frame with some key threshold values and evaluation
#' statistics of each algorithm (FNR, FPR, TSSmax, AUC, pROC, FScore,
#' Jaccard dissimilarity etc.) for the selected threshold
#' @return Writes on disk a .tif model for each partition of each algorithm
#' @return Writes in disk a .csv file with thresholds and evaluation statistics
#' of each algorithm for a given threshold
#' #' @return Writes in disk a .csv file with evaluation statistics for all
#' threshold values
#' @details See below for a description on the implementation of the algorithms
#' supported in this package.
#' \describe{
#' \item{Bioclim}{
#' Specified by \code{algo = "bioclim"} uses \code{\link[dismo]{bioclim}}
#' function in \pkg{dismo} package \insertCite{hijmans_dismo_2017}{modleR}.
#' Bioclim is the climate-envelope-model implemented by Henry Nix
#' \insertCite{nix_biogeographic_1986}{modleR}, the first species distribution
#' modelling package. It is based on climate interpolation methods and despite
#' its limitations it is still used in ecological niche modeling, specially for
#' exploration and teaching purposes
#' \insertCite{@see also @booth_bioclim_2014}{modleR}. In this package it is
#' implemented by the
#' function \code{\link[dismo]{bioclim}}, evaluated and predicted using
#' \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} also from
#' \pkg{dismo} package.
#' }
#' \item{Boosted Regression Trees (BRT)}{
#' Specified by \code{algo = "brt"}, it uses \code{\link[dismo]{gbm.step}}
#' function from \pkg{dismo} package. Runs the cross-validation procedure of
#' \insertCite{hastie_elements_2001;textual}{modleR}
#' \insertCite{@see also @elith_working_2009}{modleR}. It consists in a
#' regression modeling technique combined with the boosting method, a method for
#' combining many simple models. It is implemented by the function
#' \code{\link[dismo]{gbm.step}} as a regression with the response variable set
#' to Bernoulli distribution, evaluated and predicted using
#' \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} from
#' \pkg{dismo} package.
#' }
#' \item{Domain}{
#' Specified by \code{algo = "domain"} uses \code{\link[dismo]{domain}} function
#' from \pkg{dismo} package. Computes point-to-point similarity based on Gower
#' distance between environmental variables
#' \insertCite{carpenter_domain_1993}{modleR}.
#' \insertCite{hijmans_dismo_2017}{modleR} state that one should use it with
#' caution because it does not perform well compared to other algorithms
#' \insertCite{elith_novel_2006,hijmans_ability_2006}{modleR}. We add that it is
#'  a slow algorithm. In this package it is implemented by the function
#'  \code{\link[dismo]{domain}}, evaluated and predicted using
#'  \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} also from
#'  \pkg{dismo} package.
#' }
#' \item{Generalized Linear Model (GLM)}{
#' Specified by \code{algo = "glm"} runs a GLM with modeling presences and
#' absences as a response variable following a binomial error distribution. It
#' runs a step-wise model selection based on AIC both backward and forward
#' considering all possible combinations of predictor variables in the
#' RasterStack. In this package it is implemented using functions \code{glm} and
#'  \code{step} to fit a model and choose a model by AIC in a stepwise procedure.
#'  Model is evaluated and predicted using \code{\link[dismo]{evaluate}}
#'  function from \pkg{dismo} and \code{\link[raster]{predict}} function from
#'  \pkg{raster} package both with argument \code{type = "response"} to return
#'  values in the scale of the response variable.
#' }
#' \item{Mahalanobis}{
#' Specified by \code{algo = "mahal"} uses \code{\link[dismo]{mahal}} function
#' from \pkg{dismo} package. Corresponds to a distribution model based on
#' Mahalanobis distance, a measure of the distance between a point P and a
#' distribution D \insertCite{mahalanobis_generalized_1936}{modleR}. In this
#' package it is implemented by the function \code{\link[dismo]{mahal}},
#' evaluated and predicted using \code{\link[dismo]{evaluate}} and
#' \code{\link[dismo]{predict}} also from \pkg{dismo} package.
#' }
#' \item{Maximum Entropy (Maxent)}{
#' Specified either by \code{algo = "maxent"} or \code{algo = "maxnet"}
#' corresponding to implementation by \pkg{dismo}
#' \insertCite{hijmans_dismo_2017}{modleR} and \pkg{maxnet}
#' \insertCite{phillips_maxnet_2017}{modleR} packages respectively. Maxent is a
#' machine learning method for modeling species distributions based in
#' incomplete data allowing ENM with presence-only data
#' \insertCite{phillips_maximum_2006}{modleR}. If \code{algo = "maxent"} model
#' is fit by the function \code{\link[dismo]{maxent}}, evaluated and predicted
#' using  \code{\link[dismo]{evaluate}} and \code{\link[dismo]{predict}} also in
#'  \pkg{dismo} package. If \code{algo = "maxnet"} model is fit by the function
#'  \code{\link[maxnet]{maxnet}} from \pkg{maxnet} package, evaluated using
#'  \code{\link[dismo]{evaluate}} from \pkg{dismo} package with argument
#'  \code{type = "logistic"} and predicted using \code{\link[raster]{predict}}
#'  function from \pkg{raster} package.
#' }
#' \item{Random Forest}{
#' Specified by \code{algo = "rf"} uses \code{\link[randomForest]{tuneRF}}
#' function from \pkg{randomForest} package
#' \insertCite{liaw_classification_2002}{modleR}. Corresponds to machine
#' learning regression based on decision trees. In this package uses
#' \code{\link[randomForest]{tuneRF}} function with the optimal number of
#' variables available for splitting at each tree node (i.e. \code{mtry}) found
#' as set by parameter \code{doBest = TRUE}. Random Forest model is evaluated
#' with \code{\link[dismo]{evaluate}} function from \pkg{dismo} and predicted
#' with \code{\link[raster]{predict}} function from \pkg{raster} package.
#' }
#' \item{Support Vector Machines (SVM)}{
#' Specified either by \code{algo = "svme"} or \code{algo = "svmk"}
#' corresponding to implementation on \pkg{e1071}
#' \insertCite{meyer_e1071_2017}{modleR} and \pkg{kernlab}
#' \insertCite{karatzoglou_kernlab_2004}{modleR} packages respectively. SVM are
#'  supervised learning models that use learning algorithms for classification
#'  and regression analysis. In \pkg{e1071} package SVM is implemented through
#'  function \code{\link[e1071]{best.tune}} with method set to "\code{svm}"
#'  which uses RBF-kernel (radial basis function kernel) for classification. In
#'  \pkg{kernlab} package SVM is implemented through function
#'  \code{\link[kernlab]{ksvm}} also with RBF-kernel method (in this case the
#'  default method "\code{kbfdot}"). We expect both implementations to differ
#'  only in performance. Both \code{svme} and \code{svmk} are evaluated with
#'  \code{\link[dismo]{evaluate}} function from dismo and predicted with
#'  \code{\link[raster]{predict}} function from \pkg{raster} package.
#' }
#' }
#' @references
#'     \insertAllCited{}
#' @seealso \code{\link[dismo]{bioclim}} in \pkg{dismo} package
#' @seealso \code{\link[dismo]{domain}} in \pkg{dismo} package
#' @seealso \code{\link{do_many}}
#' @seealso \code{\link[dismo]{evaluate}} in \pkg{dismo} package
#' @seealso \code{\link[dismo]{maxent}} in \pkg{dismo} package
#' @seealso \code{\link[maxnet]{maxnet}} in \pkg{maxnet} package
#' @seealso \code{\link[dismo]{mahal}} in \pkg{dismo} package
#' @seealso \code{\link[dismo]{predict}} in \pkg{dismo} package
#' @seealso \code{\link[raster]{predict}} in \pkg{raster} package
#' @import raster
#' @import grDevices
#' @importFrom utils write.csv
#' @importFrom maxnet maxnet
#' @importFrom stats complete.cases formula glm step dist
#' @importFrom Rdpack reprompt
#' @importFrom kuenm kuenm_proc
#'
#' @param algorithm Character string of length 1 specifying the algorithm to
#' be fit: "\code{bioclim}", "\code{brt}",
#' "\code{domain}", "\code{glm}", "\code{maxent}", "\code{mahal}",
#' "\code{svme}", "\code{svmk}", "\code{rf}"
#' @param bioclim Execute bioclim algorithm from the \pkg{dismo} implementation
#' with \code{\link[dismo]{bioclim}} function
#' @param brt Execute Boosted Regression Trees with
#' \code{\link[dismo]{gbm.step}} from \pkg{dismo}
#' @param domain Execute domain from the \pkg{dismo} implementation with
#' \code{\link[dismo]{domain}} function
#' @param glm Execute GLM as suggested by the \pkg{dismo} documentation with
#' \code{\link{glm}} and \code{\link{step}}
#' @param mahal Execute Mahalanobis distance from the \pkg{dismo} implementation
#'  with \code{\link[dismo]{mahal}}
#' @param maxent Execute Maxent algorithm from the \pkg{dismo} implementation
#' with \code{\link[dismo]{maxent}} function
#' @param maxnet Execute Maxent algorithm from the \pkg{maxnet} implementation
#' with \code{\link[maxnet]{maxnet}} function
#' @param rf Execute Random forest algorithm from \pkg{randomForest} package
#' with function \code{\link[randomForest]{tuneRF}} as suggested by the
#' \pkg{dismo} documentation
#' @param svme Execute Support Vector Machines (SVM) algorithm from \pkg{e1071}
#'  package with \code{\link[e1071]{best.tune}} function
#' @param svmk Execute Support Vector Machines (SVM) algorithm from
#' \pkg{kernlab} package with \code{\link[kernlab]{ksvm}} function
#' @param project_model Logical, whether to project the models to variable sets
#' in \code{proj_data_folder} directory
#' @param proj_data_folder Path to directory with projections containing one or
#' more folders with the projection datasets (e.g. "./env/proj/proj1").
#' This directory should only contain raster files corresponding to the
#' environmental variables. If more than one projection, each projection should
#' be at one directory (e.g. "./env/proj/proj1" and "./env/proj/proj2") and
#' equivalent raster files at diferent subdirectories must have the same names
#' (e.g. "./env/proj/proj1/layer1.asc" and "./env/proj/proj2/layer1.asc")
#' @param mask A SpatialPolygonsDataFrame to be used to mask the models. This
#' mask can be used if the final area of interest is smaller than the area used
#'  for model fitting, to save disk space
#' @param write_png Logical, whether png files will be written
#' @param write_bin_cut Logical, whether binary and cut model files(.tif, .png)
#' should be written
#' @param dismo_threshold Character string indicating threshold (cut-off) to
#' transform model predictions to a binary score as in
#' \code{\link[dismo]{threshold}}:
#'  "\code{kappa}", "\code{spec_sens}", "\code{no_omission}",
#'   "\code{prevalence}", "\code{equal_sens_spec}",
#'  "\code{sensitivity}". Default value is "\code{spec_sens}"
#' @param conf_mat Logical, whether confusion tables should be written in the HD
#' @param equalize Logical, whether the number of presences and absences should be
#' equalized in randomForest and brt
#' @param proc_threshold Numeric, value from 0 to 100 that will be used as (E)
#' for partialROC calculations in \code{\link[kuenm]{kuenm_proc}}. Default is
#' \code{proc_threshold = 5}
#' @param ... Other arguments from \code{\link[kuenm]{kuenm_proc}}
#'
#' @examples
#' # run setup_sdmdata first from one species in example_occs data
#' sp <- names(example_occs)[1]
#' sp_coord <- example_occs[[1]]
#' sp_setup <- setup_sdmdata(species_name = sp,
#'                           occurrences = sp_coord,
#'                           predictors = example_vars)
#'
#' # run bioclim for one species
#' sp_any <- do_any(species_name = sp,
#'                  predictors = example_vars,
#'                  algorithm = "bioclim")
#'
#' # run do_many
#' sp_many <- do_many(species_name = sp,
#'                    predictors = example_vars,
#'                    bioclim = TRUE,
#'                    maxent = TRUE)
#' @name model_fit
#'
NULL
