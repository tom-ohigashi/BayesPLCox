#' Trace plots for posterior samples
#'
#' Generic function for trace plots of posterior samples from fitted Bayesian
#' Cox models.
#'
#' @param x A fitted model object.
#' @param ... Further arguments passed to methods.
#'
#' @return A plot object, typically a `ggplot2` or `patchwork` object.
#' @export
traceplot <- function(x, ...) {
  UseMethod("traceplot")
}
