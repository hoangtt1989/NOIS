
#' #' @keywords internal
#' qdet <- function(local_q = 0.1, nz) {
#'     ret_val <- max(1, ceiling(local_q * nz))
#'     return(ret_val)
#' }

#' Tuning modified BIC for NOIS
#'
#' Find the optimal number of detected outliers for NOIS using modified BIC
#'
#' @param data An input \code{data.frame}.
#' @param q_tst The grid of test values for the pooled q.
#' @param bias_correct A logical indicating bias correction.
#' @param parallel A logical indicating parallel computation. A parallel backend must be registered.
#' @param return_fit A logical indicating whether a model with the optimal q should be fit (using \code{NOIS_fit}) and returned.
#' @param ... Additional arguments passed to \code{NOIS_fit}.
#' @return A list with the following components.
#' \item{\code{q_tst}}{The grid of test values.}
#' \item{\code{min_q}}{The optimal number of detected outliers.}
#' \item{\code{opt_fit}}{A \code{NOIS_fit} using the optimal number of detected outliers.}
#' \item{\code{BIC_vals}}{The BIC for each fit.}
#' \item{\code{df_vals}}{The estimated degrees of freedom for each fit.}
#' @family NOIS BIC functions.
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @export
BIC_tuner <- function(data, q_tst = 1:floor(nrow(data)/3), bias_correct = T, parallel = F, return_fit = T, ...) {

    `%fun%` <- ifelse(parallel, `%dopar%`, `%do%`)
    i <- NULL
    BIC_fits <- foreach::foreach(i = 1:length(q_tst)) %fun% {
        fit <- NOIS_fit(data, pool_q = q_tst[i], ...)
        BIC_val <- BIC_NOIS(fit, bias_correct = T)
        return(BIC_val)
    }
    BIC_vals <- purrr::map_dbl(BIC_fits, "BIC")
    df_vals <- purrr::map_dbl(BIC_fits, "df")
    min_ind <- q_tst[which.min(BIC_vals)]
    if (min_ind == q_tst[1] | min_ind == q_tst[length(q_tst)]) {
        warning("Minimum BIC is at the edge of the grid")
    }
    if (return_fit) {
      opt_fit <- NOIS_fit(data, pool_q = min_ind, ...)
    } else {
      opt_fit <- NULL
    }
    return(list(q_tst = q_tst, min_q = min_ind, opt_fit = opt_fit, BIC_vals = BIC_vals, df_vals = df_vals))
}

#' Modified BIC for NOIS
#' Compute the modified Bayesian Information Criterion for a NOIS fit
#'
#' @param NOIS_fit A \code{NOIS_fit}.
#' @param bias_correct A logical indicating bias correction.
#'
#' @return A list with the following components.
#' \item{\code{BIC}}{The value of modified BIC.}
#' \item{\code{df}}{The estimated degrees of freedom.}
#' @family NOIS BIC functions.
#' @export
BIC_NOIS <- function(NOIS_fit, bias_correct = T) {
    if (class(NOIS_fit) != "NOIS_fit") {
        stop("Input must be a NOIS_fit")
    }
    # npts <- length(NOIS_fit$y)
    npts <- length(NOIS_fit$fit_df$y)
    smooth_df <- do.call(rbind, lapply(1:npts, function(ind) {
        # tmp <- stats::dnorm(NOIS_fit$x[ind] - NOIS_fit$x, 0, NOIS_fit$CV$pool_h)
        tmp <- stats::dnorm(NOIS_fit$fit_df$x[ind] - NOIS_fit$fit_df$x, 0, NOIS_fit$CV$pool_h)
        ret <- tmp/sum(tmp)
        return(ret)
    }))
    if (bias_correct) {
        # fit <- NOIS_fit$bias_pool_fit
        fit <- NOIS_fit$fit_df$bias_fit
        # df calculation is different for bias fit
        m <- sum(diag(2 * smooth_df - smooth_df %*% smooth_df))
    } else {
        fit <- NOIS_fit$fit_df$fit
        m <- sum(diag(smooth_df))
    }
    BIC_val <- (npts - m) * log(sum((NOIS_fit$fit_df$y_adj - fit)^2)/(npts - m)) + (NOIS_fit$pool_q + 1) * (log(npts -
        m) + 1)
    return(list(BIC = BIC_val, df = m))
}


#' Fitting NOIS to data
#'
#' Fit NOIS to data and return a \code{NOIS_fit} object
#' @param data A \code{data.frame}.
#' @param x The name of the column in \code{data} containing the 'x' values.
#' @param y The name of the column in \code{data} containing the 'y' values.
#' @param CV_method The type of cross-validation to use. Possible types are \code{c('MCV', 'LOOCV', 'PCV', 'none')}.
#' @param first_h Bandwidth for first model fit. Only used when \code{CV_method = 'none'}.
#' Default value is the theoretically optimal bandwidth of \eqn{n^{-1/5}}.
#' @param pool_h Bandwidth for pooled model fit. Only used when \code{CV_method = 'none'}.
#' Default value is the theoretically optimal bandwidth of \eqn{n^{-1/5}}.
#' @param local_q Fraction of points detected as outliers for each unpooled fit.
#' @param pool_q Pooled outlier detection. For numerics \eqn{\in (0, 1)},
#' this is the fraction of points detected as outliers.
#' For integers \eqn{\ge 1} this is the number of points detected as outliers.
#' @param tol Tolerance for each unpooled fit.
#' @param maxit Maximum number of iterations for each individual kernel smoothing fit.
#' @param ... Additional arguments passed to the \code{CV_method}.
#' @return An object of class '\code{NOIS_fit}' that is a list with the following components.
#' \item{\code{fit_df}}{A \code{data_frame} with columns specifying the original 'x' and 'y' values, the pooled adjusted 'y' values, a logical indicating whether the observation is an outlier, the pooled non-bias corrected and bias corrected fits, and the non-robust non-bias corrected and bias corrected fits.}
#' \item{\code{local_fit}}{Unpooled NOIS fits.}
#' \item{\code{local_gamma}}{The unpooled \eqn{\gamma} estimates.}
#' \item{\code{pool_gamma}}{The pooled \eqn{\gamma} estimate.}
#' \item{\code{local_q}}{Fraction of points detected as outliers for each unpooled fit.}
#' \item{\code{pool_q}}{Pooled outlier detection. For numerics \eqn{\in (0, 1)},
#' this is the fraction of points detected as outliers.
#' For integers \eqn{\ge 1} this is the number of points detected as outliers.}
#' \item{\code{pool_nonout}}{The positions of the clean observations.}
#' \item{code{CV}}{A list with cross-validation information.}
#' \item{code{conv}}{A list with convergence information.}
#' @examples
#' ###generate some random data and introduce outliers
#' set.seed(123)
#' npts <- 100
#' nout <- floor(.1*npts)
#' xt <- seq(from=0, to=2*pi, length.out=npts)
#' gaussnoise <- rnorm(npts)
#' outliers <- sample(floor(npts/2):npts, size=nout)
#' randpts <- runif(nout, min=5, max=7)
#' yt <- sin(xt) + gaussnoise
#' yt[outliers] <- yt[outliers] + randpts
#' sine_data <- data.frame(x = xt, y = yt)
#' ###fit NOIS to this data
#' sine_fit <- NOIS_fit(sine_data, 'x', 'y', pool_q = .1, CV_method = 'LOOCV')
#' @family NOIS CV functions
#' @useDynLib NOIS
#' @importFrom Rcpp sourceCpp
#' @export
NOIS_fit <- function(data, x = "x", y = "y", CV_method = "LOOCV", first_h = NULL, pool_h = NULL, local_q = 0.1,
    pool_q = 0.1, tol = 1e-07, maxit = 200, ...) {

    yy <- data[[y]]
    xx <- data[[x]]

    # checking inputs
    if (!is.data.frame(data)) {
        stop("data must be a data.frame")
    }
    if (any(!(c(x, y) %in% colnames(data)))) {
        stop("x and y must be column names in data")
    }
    if (local_q <= 0 | local_q >= 1) {
        stop("local_q must be between 0 and 1")
    }
    if (pool_q <= 0) {
        stop("pool_q must be greater than 0")
    }
    if (!all(is.finite(xx)) | !all(is.finite(yy))) {
        stop("Bad values in either x or y")
    }

    # initializing vectors, parameters
    nn <- nrow(data)
    cond_check <- rep(0, nn)
    gamma_curr <- matrix(0, nrow = nn, ncol = nn)
    local_fit <- rep(0, nn)
    qq <- rep(0, nn)
    converged <- rep(0, nn)
    iter <- rep(0, nn)
    cond_check <- rep(1, nn)

    ptm <- proc.time()
    if (CV_method %in% c("LOOCV", "MCV", "PCV")) {
        first_CV <- switch(CV_method, LOOCV = LOOCV_grid(xx, yy, ...), MCV = MCV_grid(xx, yy, ...), PCV = PCV_grid(xx, yy, ...))
        first_h <- first_CV$min_h
        first_hgrid <- first_CV$hgrid
        if (first_h == first_hgrid[1] | first_h == first_hgrid[length(first_hgrid)]) {
            warning("Optimal first bandwidth is at the end of the grid")
        }
    } else if (CV_method == "none") {
        if (is.null(first_h)) {
            warning("No CV method supplied or bandwidth supplied. Using theoretically optimal bandwidth for first h.")
            first_h <- nrow(data)^(-1/5)
        }
        first_hgrid <- NULL
    } else {
        stop("Supply valid CV_method - LOOCV, MCV, PCV or none")
    }

    nw_ests <- nwvector(xx, yy, first_h)
    bias_nw_ests <- biasnwvector(xx, yy, nw_ests, first_h)

    loop_fit <- NOIS_loop(xx, yy, first_h, local_q, tol, maxit)

    local_fit <- loop_fit$local_fit
    gamma_curr <- loop_fit$gamma_curr
    converged <- loop_fit$converged
    cond_check <- loop_fit$cond_check
    iter <- loop_fit$iter
    qq_inner <- loop_fit$qq_inner

    # get max
    gam_max <- apply(gamma_curr, 1, function(x) {
        x[which.max(abs(x))]
    })

    # index for gamma vectors
    if (pool_q < 1) {
        pool_q <- min(floor(pool_q * nn), length(which(gam_max != 0)))
    } else {
        # so that we can supply pool_q as an integer if needed
        pool_q <- min(pool_q, length(which(gam_max != 0)))
    }
    gam_val <- quantile_thresh(gam_max, pool_q)
    gam_ind <- which(gam_val != 0)

    # pooled y adj
    pool_y_adj <- yy - gam_val
    # pooled cv
    if (CV_method %in% c("LOOCV", "MCV", "PCV")) {
        pool_CV <- switch(CV_method, LOOCV = LOOCV_grid(xx, pool_y_adj, ...), MCV = MCV_grid(xx, pool_y_adj, ...), PCV = PCV_grid(xx, pool_y_adj, ...))
        pool_h <- pool_CV$min_h
        pool_hgrid <- pool_CV$hgrid
        if (pool_h == pool_hgrid[1] | pool_h == pool_hgrid[length(pool_hgrid)]) {
            warning("Optimal pooled bandwidth is at the end of the grid")
        }
    } else if (CV_method == "none") {
        if (is.null(first_h)) {
            warning("No CV method supplied or bandwidth supplied. Using theoretically optimal bandwidth for first h.")
            pool_h <- nrow(data)^(-1/5)
        }
        pool_hgrid <- NULL
    }

    # pooled fit
    pool_fit <- nwvector(xx, pool_y_adj, pool_h)
    bias_pool_fit <- biasnwvector(xx, pool_y_adj, pool_fit, pool_h)
    pool_nonout <- setdiff(1:nn, gam_ind)
    etm <- proc.time() - ptm

    # CV
    CV <- list(first_h = first_h, pool_h = pool_h, first_hgrid = first_hgrid, pool_hgrid = pool_hgrid)
    # convergence
    conv <- list(iter = iter, time = etm, converged = converged, cond_check = cond_check)

    # data_frame with relevant fits
    outlier_index <- rep(FALSE, length(xx))
    outlier_index[gam_ind] <- TRUE
    out_df <- tibble::data_frame(x = xx, y = yy, y_adj = pool_y_adj, outlier = outlier_index, fit = pool_fit,
                         bias_fit = bias_pool_fit, nr_fit = nw_ests, nr_bias_fit = bias_nw_ests)

    model_output <- list(fit_df = out_df, local_fit = local_fit, local_gamma = gamma_curr, pool_gamma = gam_val, local_q = qq_inner, pool_q = pool_q, pool_nonut = pool_nonout, CV = CV, conv = conv)
    # model_output <- list(local_fit = local_fit, pool_fit = pool_fit, bias_pool_fit = bias_pool_fit, first_fit = nw_ests,
    #     bias_first_fit = bias_nw_ests, local_gamma = gamma_curr, pool_gamma = gam_val, pool_outlier = gam_ind,
    #     local_q = qq_inner, pool_q = pool_q, pool_nonout = pool_nonout, x = xx, y_adj = pool_y_adj, y = yy, CV = CV, conv = conv)
    class(model_output) <- "NOIS_fit"
    model_output
}


#' Print method for a NOIS fit
#' @param x A \code{NOIS_fit}.
#' @param ... Not used.
#' @export
print.NOIS_fit <- function(x, ...) {
    cat("Number of detected outliers =", length(which(x$fit_df$outlier)), "\nNumber of observations =", length(x$fit_df$y), "\nTime =",
        x$conv$time[3], "\nConvergence =", all(as.logical(x$conv$converged)), "\nMSE =", mean((x$fit_df$y_adj - x$fit_df$fit)^2),
        "\nBias corrected MSE =", mean((x$fit_df$y_adj - x$fit_df$bias_fit)^2), "\nFirst optimal bandwidth =", x$CV$first_h,
        "\nPooled optimal bandwidth =", x$CV$pool_h)
}


#' Create a plot highlighting the outliers and a NOIS curve.
#'
#' Use a \code{NOIS_fit} to create a plot of the outliers and the NOIS curve.
#'
#' @param NOIS_fit A \code{NOIS_fit}
#' @param color An aesthetic specification for \code{geom_point}.
#' This must be a string correpsponding to a factor variable from a \code{NOIS_df} (usually \code{'outlier'}).
#' @param ... Additional parameters passed to \code{aes_string} in \code{geom_point}.
#' @param bias_correct Logical specfiying bias correction for the NOIS curve.
#' @return A \code{ggplot}.
#' @export
outlier_plot <- function(NOIS_fit, color = "outlier", ..., bias_correct = T) {
    if (class(NOIS_fit) != "NOIS_fit") {
        stop("Input must be a NOIS_fit")
    }
    fit_type <- ifelse(bias_correct, "bias_fit", "fit")
    # df <- NOIS_df(NOIS_fit)
    df <- NOIS_fit$fit_df
    pt <- ggplot2::ggplot(df) + ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y", color = color, ...)) +
        ggplot2::geom_line(ggplot2::aes_string(x = "x", y = fit_type))
    pt
}

