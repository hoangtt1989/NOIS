
#' @keywords internal
qdet <- function(local_q = 0.1, nz) {
    ret_val <- max(1, ceiling(local_q * length(nz)))
    return(ret_val)
}

#' Tuning modified BIC for NOIS
#'
#' Find the optimal number of detected outliers for NOIS using modified BIC
#'
#' @param data An input \code{data.frame}.
#' @param q_tst The grid of test values for the pooled q.
#' @param bias_correct A logical indicating bias correction.
#' @param parallel A logical indicating parallel computation. A parallel backend must be registered.
#' @param ... Additional arguments passed to \code{NOIS_fit}.
#' @return A list with the following components.
#' \item{\code{q_tst}}{The grid of test values.}
#' \item{\code{min_q}}{The optimal number of detected outliers.}
#' \item{\code{BIC_vals}}{The individual outputs from \code{BIC.NOIS_fit}.}
#' @family NOIS BIC functions.
#' @importFrom foreach %dopar%
#' @importFrom foreach %do%
#' @export
BIC_tuner <- function(data, q_tst = 1:floor(nrow(data)/3), bias_correct = T, parallel = F, ...) {

    `%fun%` <- ifelse(parallel, `%dopar%`, `%do%`)

    BIC_fits <- foreach::foreach(i = 1:length(q_tst)) %fun% {
        fit <- NOIS_fit(data, pool_q = q_tst[i], ...)
        BIC_val <- BIC_NOIS(fit, bias_correct = T)
        return(BIC_val)
    }

    min_ind <- q_tst[which.min(purrr::map_dbl(BIC_fits, "BIC"))]
    if (min_ind == q_tst[1] | min_ind == q_tst[length(q_tst)]) {
        warning("Minimum BIC is at the edge of the grid")
    }
    return(list(q_tst = q_tst, min_q = min_ind, BIC_vals = BIC_fits))
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
    npts <- length(NOIS_fit$y)
    smooth_df <- do.call(rbind, lapply(1:npts, function(ind) {
        # tmp <- gausskern(NOIS_fit$x[ind] - NOIS_fit$x, NOIS_fit$pool_h)
        tmp <- qnorm(NOIS_fit$x[ind] - NOIS_fit$x, 0, NOIS_fit$pool_h)
        ret <- tmp/sum(tmp)
        return(ret)
    }))
    if (bias_correct) {
        fit <- NOIS_fit$bias_pool_fit
        # df calculation is different for bias fit
        m <- sum(diag(2 * smooth_df - smooth_df %*% smooth_df))
    } else {
        fit <- NOIS_fit$local_fit
        m <- sum(diag(smooth_df))
    }
    BIC_val <- (npts - m) * log(sum((NOIS_fit$y_adj - fit)^2)/(npts - m)) + (NOIS_fit$pool_q + 1) * (log(npts -
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
#' \item{\code{local_fit}}{Unpooled NOIS fits.}
#' \item{\code{pool_fit}}{Pooled NOIS fit.}
#' \item{\code{bias_pool_fit}}{Bias corrected pooled NOIS fit.}
#' \item{\code{first_fit}}{The first (non-robust) kernel smoothing fit.}
#' \item{\code{bias_first_fit}}{The first (non-robust) bias corrected kernel smoothing fit.}
#' \item{\code{local_gamma}}{The unpooled \eqn{\gamma} estimates.}
#' \item{\code{pool_gamma}}{The pooled \eqn{\gamma} estimate.}
#' \item{\code{pool_outlier}}{The positions of the pooled outlier estimates.}
#' \item{\code{local_q}}{Fraction of points detected as outliers for each unpooled fit.}
#' \item{\code{pool_q}}{Pooled outlier detection. For numerics \eqn{\in (0, 1)},
#' this is the fraction of points detected as outliers.
#' For integers \eqn{\ge 1} this is the number of points detected as outliers.}
#' \item{\code{pool_nonout}}{The positions of the clean observations.}
#' \item{\code{x}}{Original 'x' values.}
#' \item{\code{y_adj}}{Pooled adjusted 'y' values.}
#' \item{\code{y}}{Original 'x' values.}
#' \item{\code{first_h}}{Bandwidth used for first (non-robust) kernel smoothing fit.}
#' \item{\code{pool_h}}{Bandwidth used for pooled NOIS fit.}
#' \item{\code{first_hgrid}}{Grid of bandwidths used for first (non-robust) kernel smoothing fit.}
#' \item{\code{pool_hgrid}}{Grid of bandwidths used for pooled NOIS fit.}
#' \item{\code{iter}}{Number of iterations for each unpooled fit.}
#' \item{\code{time}}{Elapsed time.}
#' \item{\code{converged}}{Convergence of each unpooled fit. 0 means the fit did not converge in \code{maxit}, 1 means it did.}
#' \item{\code{cond_check}}{The infinity norm of each unpooled \eqn{\gamma} estimate.}
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


    if (CV_method %in% c("LOOCV", "MCV", "PCV")) {
        if (CV_method == "LOOCV") {
            first_CV <- LOOCV_grid(xx, yy, ...)
        } else if (CV_method == "MCV") {
            first_CV <- MCV_grid(xx, yy, ...)
        } else if (CV_method == "PCV") {
            first_CV <- PCV_grid(xx, yy, ...)
        }
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

    ptm <- proc.time()
    for (jj in 1:nn) {

        xx_inner <- xx[jj]

        kernlist <- qnorm(xx_inner - xx, 0, first_h)
        # kernlist <- gausskern(xx_inner - xx, first_h)
        nz_ind <- which(kernlist != 0 & kernlist >= 1e-20)
        kern_nz <- kernlist[nz_ind]
        kern_nzsqrt <- sqrt(kern_nz)
        kern_nzsqrtinv <- 1/kern_nzsqrt
        qq[jj] <- qdet(local_q, kern_nz)

        gamma_inner <- rep(0, nn)
        qq_inner <- qq[jj]

        for (ii in 1:maxit) {
            gamma_next <- rep(0, nn)
            yy_adj <- yy - gamma_inner
            local_inner <- nwestimator(xx_inner, xx, yy_adj, first_h)
            rr <- kern_nzsqrt * ((yy - local_inner))[nz_ind]
            gamma_next[nz_ind] <- kern_nzsqrtinv * quantile_thresh(rr, qq_inner)
            cond_inner <- max(abs(gamma_next - gamma_inner))
            gamma_inner <- gamma_next

            if (max(abs(cond_inner)) <= tol) {
                converge_inner <- T
                break
            } else {
                converge_inner = F
            }
        }

        local_fit[jj] <- local_inner
        gamma_curr[, jj] <- gamma_inner
        converged[jj] <- converge_inner
        cond_check[jj] <- cond_inner
        iter[jj] <- ii

        if (ii == maxit) {
            warning(paste("Model did not converge at j =", jj))
        }
    }
    etm <- proc.time() - ptm

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
        if (CV_method == "LOOCV") {
            pool_CV <- LOOCV_grid(xx, pool_y_adj, ...)
        } else if (CV_method == "MCV") {
            pool_CV <- MCV_grid(xx, pool_y_adj, ...)
        } else if (CV_method == "PCV") {
            pool_CV <- PCV_grid(xx, pool_y_adj, ...)
        }
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

    model_output <- list(local_fit = local_fit, pool_fit = pool_fit, bias_pool_fit = bias_pool_fit, first_fit = nw_ests,
        bias_first_fit = bias_nw_ests, local_gamma = gamma_curr, pool_gamma = gam_val, pool_outlier = gam_ind,
        local_q = qq, pool_q = pool_q, pool_nonout = pool_nonout, x = xx, y_adj = pool_y_adj, y = yy, first_h = first_h,
        pool_h = pool_h, first_hgrid = first_hgrid, pool_hgrid = pool_hgrid, iter = iter, time = etm, converged = converged,
        cond_check = cond_check)
    class(model_output) <- "NOIS_fit"
    model_output
}

#' Print method for a NOIS fit
#' @param x A \code{NOIS_fit}.
#' @param ... Not used.
#' @export
print.NOIS_fit <- function(x, ...) {
    cat("Number of detected outliers =", length(x$pool_outlier), "\nNumber of observations =", length(x$y), "\nTime =",
        x$time[3], "\nConvergence =", all(as.logical(x$converged)), "\nMSE =", mean((x$y_adj - x$pool_fit)^2),
        "\nBias corrected MSE =", mean((x$y_adj - x$bias_pool_fit)^2), "\nFirst optimal bandwidth =", x$first_h,
        "\nPooled optimal bandwidth =", x$pool_h)
}


#' Construct a \code{data_frame} from a \code{NOIS_fit}.
#'
#' Construct a \code{data_frame} using output from a \code{NOIS_fit}.
#'
#' @param NOIS_fit A \code{NOIS_fit}
#' @return A \code{data_frame} with the following columns.
#' \item{\code{x}}{'x' values.}
#' \item{\code{y}}{'y' values.}
#' \item{\code{y_adj}}{Pooled adjusted 'y' values.}
#' \item{\code{fit}}{Pooled NOIS fit without bias correction.}
#' \item{\code{bias}}{Pooled NOIS fit with bias correction.}
#' \item{\code{outlier}}{A logical indicating whether this point is an outlier.}
#' @export
NOIS_df <- function(NOIS_fit) {
    if (class(NOIS_fit) != "NOIS_fit") {
        stop("Input must be a NOIS_fit")
    }
    outlier_index <- rep(FALSE, length(NOIS_fit$x))
    outlier_index[NOIS_fit$pool_outlier] <- TRUE
    df <- with(NOIS_fit, tibble::data_frame(x = x, y = y, y_adj = y_adj, fit = pool_fit, bias_fit = bias_pool_fit,
        outlier = outlier_index))
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
    df <- NOIS_df(NOIS_fit)
    pt <- ggplot2::ggplot(df) + ggplot2::geom_point(ggplot2::aes_string(x = "x", y = "y", color = color, ...)) +
        ggplot2::geom_line(ggplot2::aes_string(x = "x", y = fit_type))
    pt
}

