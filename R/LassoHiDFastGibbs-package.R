#' @title Fast Gibbs Samplers: Blocked, Partially Collapsed, and Nested Gibbs Samplers
#'
#' @description Provides fast and scalable Gibbs sampling algorithms for
#' Bayesian penalized linear regression models in high-dimensional
#' settings. The package implements efficient partially collapsed
#' and nested Gibbs samplers for Bayesian Lasso, with a focus on
#' computational efficiency when the number of predictors is large
#' relative to the sample size.
#' @name LassoHiDFastGibbs
#' @docType package
@useDynLib LassoHiDFastGibbs, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
