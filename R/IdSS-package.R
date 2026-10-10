#' IdSS: Identification and Analysis of Dynamic Structural Shocks
#'
#' `IdSS` contains the procedures and datasets used in the courses
#' *Introduction to Time Series* and *The Identification of Dynamic Structural
#' Shocks*. Its public interface covers simulation, estimation, impulse-response
#' analysis, forecast-error variance decompositions, and several approaches to
#' structural-shock identification.
#'
#' @keywords internal
#' @importFrom graphics abline legend lines par plot polygon text
#' @importFrom grDevices rainbow
#' @importFrom MASS mvrnorm Null
#' @importFrom Matrix rankMatrix
#' @importFrom stats arima coef cor embed lm na.omit optim qnorm quantile
#'   residuals rnorm rt runif sd var
#' @importFrom stringr str_replace
#' @importFrom utils tail
#' @importFrom vars Acoef Bcoef VAR
"_PACKAGE"
