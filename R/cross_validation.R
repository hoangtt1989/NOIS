#' @keywords internal
weightfun <- function(input, samp_quant) {
    output <- rep(1, length(input))
    output[input <= samp_quant[1] | input >= samp_quant[2]] <- 0
    return(output)
}

#' LOOCV
#' @keywords internal
LOOCV <- function(x, y, bandwidth, samp_quant) {
    lenx <- length(x)
    
    LOOCV_loop <- function(input) {
        est_val <- nwestimator(x[input], x[-input], y[-input], bandwidth)
        sq_error <- (y[input] - est_val)^2
        return(sq_error)
    }
    
    cv <- vapply(1:lenx, LOOCV_loop, numeric(1))
    wts <- weightfun(x, samp_quant)
    cv <- cv * wts
    
    return(sum(cv)/lenx)
}

#' Leave-one-out cross-validation
#'
#' Evaluate leave-one-out cross-validation (LOOCV) for the Nadaraya-Watson estimator over a grid
#'
#' @param x 'x' values.
#' @param y 'y' values.
#' @param hgrid The grid of test bandwidths.
#' @param weight_probs Specification for what percentage of points are used in the weight function.
#'
#' @return The minimum bandwidth h, the minimum CV error and the grid of test bandwidths.
#' @family NOIS CV functions
#' @export
LOOCV_grid <- function(x, y, hgrid = seq(from = 0.05, to = 3, length.out = 80), weight_probs = c(0.05, 0.95)) {
    samp_quant <- stats::quantile(x, weight_probs)
    cvs <- vapply(hgrid, function(input) {
        LOOCV(x, y, input, samp_quant)
    }, numeric(1))
    min_index <- which.min(cvs)
    min_h <- hgrid[min_index]
    return(list(min_h = min_h, min_CV = cvs[min_index], hgrid = hgrid))
}


#' MCV
#' @keywords internal
MCV <- function(nextmx, nextmy, firstmx, firstmy, bandwidth, samp_quant) {
    
    lenmx <- length(nextmx)
    
    loopfun <- function(input) {
        est_val <- nwestimator(nextmx[input], firstmx, firstmy, bandwidth)
        if (is.nan(est_val)) {
            est_val <- 0
        }
        sq_error <- (nextmy[input] - est_val)^2
        return(sq_error)
    }
    
    ecv <- vapply(1:lenmx, loopfun, numeric(1))
    wts <- weightfun(nextmx, samp_quant)
    ecv <- ecv * wts
    output <- sum(ecv)/(lenmx)
    return(output)
    
}

#' Modified cross-validation
#'
#' Evaluate modified cross-validation (MCV) for the Nadaraya-Watson estimator over a grid
#'
#' @param x 'x' values.
#' @param y 'y' values.
#' @param m The fraction of observations used for training. This is the most sensitive parameter.
#' @param a Constant for tuning grid (not recommended to change).
#' @param b Constant for tuning grid (not recommended to change).
#' @param eps0 Constant for tuning grid (not recommended to change).
#' @param gridlen The number of grid points.
#' @param weight_probs Specification for what percentage of points are used in the weight function.
#'
#' @return The minimum bandwidth h, the minimum CV error and the grid of test bandwidths.
#' @family NOIS CV functions
#' @export
MCV_grid <- function(x, y, m = 0.6, a = 0.75, b = 10, eps0 = 1/175, gridlen = 100, weight_probs = c(0.05, 0.95)) {
    samp_quant <- stats::quantile(x, probs = weight_probs)
    m = floor(m * length(x))
    firstmx <- x[1:m]
    firstmy <- y[1:m]
    nextmx <- x[(m + 1):length(x)]
    nextmy <- y[(m + 1):length(y)]
    hgrid <- seq(from = a * m^(-1/5 - eps0), to = b * m^(-1/5 + eps0), length.out = gridlen)
    ecvs <- vapply(hgrid, function(input) {
        MCV(nextmx, nextmy, firstmx, firstmy, input, samp_quant)
    }, numeric(1))
    min_index <- which.min(ecvs)
    min_h <- hgrid[min_index] * (m/length(x))^(1/5)
    return(list(min_h = min_h, min_CV = ecvs[min_index], hgrid = hgrid))
}


#' PCV
#' @keywords internal
PCV <- function(data, bandwidth, g, samp_quant) {
    n <- nrow(data)
    dat_num <- g * floor(n/g)
    dat_sub <- data[1:dat_num, ]
    
    dat_grp <- list()
    cnt_vec <- seq(1, dat_num, dat_num/g)
    dat_grp <- lapply(cnt_vec, function(curr_num) {
        ret <- dat_sub[curr_num:(curr_num + (dat_num/g) - 1), ]
    })
    
    cv_err <- vapply(dat_grp, function(dat_in) {
        ret <- LOOCV(dat_in$x, dat_in$y, bandwidth, samp_quant)
    }, numeric(1))
    
    pcv_err <- mean(cv_err)
}

#' Partitioned cross-validation
#'
#' Evaluate partitioned cross-validation (PCV) for the Nadaraya-Watson estimator over a grid
#'
#' @param x 'x' values.
#' @param y 'y' values.
#' @param hgrid The grid of test bandwidths.
#' @param g The number of points used in each training subgroup. This is the most sensitive parameter.
#' @param weight_probs Specification for what percentage of points are used in the weight function.
#'
#' @return The minimum bandwidth h, the minimum CV error and the grid of test bandwidths.
#'
#' @family NOIS CV functions
#' @export
PCV_grid <- function(x, y, hgrid = seq(from = 0.05, to = 3, length.out = 100), g = floor(length(x)/5), weight_probs = c(0.05, 
    0.95)) {
    samp_quant <- stats::quantile(x, probs = weight_probs)
    data <- data.frame(x = x, y = y)
    cvs <- vapply(hgrid, function(input) {
        PCV(data, input, g, samp_quant)
    }, numeric(1))
    min_index <- which.min(cvs)
    min_h <- hgrid[min_index]
    min_h <- min_h * g^(-1/5)
    return(list(min_h = min_h, min_CV = cvs[min_index], hgrid = hgrid))
}
