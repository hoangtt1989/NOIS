#####check src/NWkernel.cpp for the Rcpp versions of these functions
#' #' @keywords internal
#' gausskern <- function(input, bandwidth) {
#'     h <- bandwidth
#'     out <- 1/(h * sqrt(2 * pi)) * exp(-input^2/(2 * h^2))
#'     return(out)
#' }
#' #' @keywords internal
#' nwestimator <- function(inputval, xvals, yvals, bandwidth) {
#'     kernel_val <- dnorm(inputval - xvals, 0, bandwidth)
#'     # kernel_val <- gausskern(inputval - xvals, bandwidth)
#'     est_value <- sum(kernel_val * yvals)/sum(kernel_val)
#'     return(est_value)
#' }
#' #' @keywords internal
#' nwvector <- function(x, y, bandwidth) {
#'     nw_ests <- vapply(x, function(input) {
#'         nwestimator(input, x, y, bandwidth)
#'     }, numeric(1))
#'     return(nw_ests)
#' }
#' #' @keywords internal
#' biasnwestimator <- function(inputval, xvals, yvals, bandwidth, inputnw, nwvals, shift_sq = FALSE) {
#'     # kernel_val <- gausskern(inputval - xvals, bandwidth)
#'     kernel_val <- dnorm(inputval - xvals, 0, bandwidth)
#'     shift <- yvals - nwvals + inputnw
#'     if (shift_sq == TRUE) {
#'         shift <- shift^2
#'     }
#'     est_value <- sum(kernel_val * shift)/sum(kernel_val)
#'     return(est_value)
#' }
#' #' @keywords internal
#' biasnwvector <- function(x, y, nwvals, bandwidth) {
#'     index <- 1:length(x)
#'     nwfunc <- function(input) {
#'         out <- biasnwestimator(x[input], x, y, bandwidth, nwvals[input], nwvals)
#'         return(out)
#'     }
#'     nw_ests <- vapply(index, nwfunc, numeric(1))
#'     return(nw_ests)
#' }

