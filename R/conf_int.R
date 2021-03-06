#' @keywords internal
pred_resids_func <- function(x, theta, first_resid, sd_est, bandwidth, bias_correct) {
  npts <- length(x)
  resid.resamp <- sample(first_resid, size = npts, replace = TRUE)
  newy <- theta + sd_est * resid.resamp
  newfit <- nwvector(x, newy, bandwidth)
  if (bias_correct == TRUE) {
    newfit <- biasnwvector(x, newy, newfit, bandwidth)
  }
  return(newfit)
}

#' @keywords internal
cvloop <- function(input, x, y, bandwidth) {
  est_val <- nwestimator(x[input], x[-input], y[-input], bandwidth)
  return(est_val)
}

#' @keywords internal
biascvloop <- function(input, x, y, nwvals, bandwidth, shift_sq = FALSE) {
  est_val <- biasnwestimator(x[input], x[-input], y[-input], bandwidth, nwvals[input], nwvals[-input], shift_sq = shift_sq)
  return(est_val)
}

#' Predictive residuals bootstrap pointwise confidence bands for '\code{NOIS_fit}' objects
#'
#' @param NOIS_fit A \code{NOIS_fit}.
#' @param conf_level The significance level.
#' @param fit_type The type of fit to use for confidence bands.
#' Valid types are \code{c('NOIS', 'regular')}, where regular is the non-robust fit..
#' @param bias_correct A logical indicating usage of bias correction.
#' @param parallel A logical indicating parallel computation. A backend must be registered first.
#' @param B Number of bootstrap replicates.
#' @return A list with the following components.
#' \item{\code{up_predicted}}{The upper band.}
#' \item{\code{low_predicted}}{The lower band.}
#' @family NOIS confidence bands
#' @importFrom foreach %dopar%
#' @export
pred_resid_BS_confint <- function(NOIS_fit, conf_level = 0.05, fit_type = "NOIS", bias_correct = T, parallel = F,
    B = 500) {
    x <- NOIS_fit$fit_df$x
    if (fit_type == "NOIS") {
        y <- NOIS_fit$fit_df$y_adj
        bandwidth <- NOIS_fit$CV$pool_h
        theta <- NOIS_fit$fit_df$fit
        if (bias_correct == T) {
            theta <- NOIS_fit$fit_df$bias_fit
            nw_est <- NOIS_fit$fit_df$fit
        }
    } else if (fit_type == "regular") {
        y <- NOIS_fit$y
        bandwidth <- NOIS_fit$CV$first_h
        theta <- NOIS_fit$fit_df$nr_fit
        if (bias_correct == T) {
            theta <- NOIS_fit$fit_df$nr_bias_fit
            nw_est <- NOIS_fit$fit_df$nr_fit
        }
    } else {
        stop("Must supply valid fit_type: NOIS or regular.")
    }

    loo_est <- sapply(1:length(x), cvloop, x, y, bandwidth)

    if (bias_correct == TRUE) {
        loo_est <- sapply(1:length(x), biascvloop, x, y, nw_est, bandwidth)
    }

    sd_estimator <- function(x, y, h, kernel_fit) {
        M_x <- sapply(1:length(x), cvloop, x, y^2, h)
        if (bias_correct == TRUE) {
            M_x <- sapply(1:length(x), biascvloop, x, y, nw_est, bandwidth, shift_sq = TRUE)
        }
        sd_est <- sqrt(M_x - kernel_fit^2)
    }

    sd_est <- sd_estimator(x, y, bandwidth, loo_est)

    fitted_resid <- (y - loo_est)/sd_est

    center_resid <- fitted_resid - mean(fitted_resid)


    if (parallel == T) {
        `%fun%` <- doRNG::`%dorng%`
        i <- NULL
        ret <- foreach::foreach(i = 1:B, .combine = cbind) %fun% {
            loop_ret <- pred_resids_func(x, theta, center_resid, sd_est, bandwidth, bias_correct)
            return(loop_ret)
        }
        resids_boots <- ret
    } else {
        resids_boots <- replicate(B, pred_resids_func(x, theta, center_resid, sd_est, bandwidth, bias_correct))
    }

    roots <- apply(resids_boots, 2, function(x) {
        theta - x
    })

    q_lower <- apply(roots, 1, stats::quantile, probs = conf_level/2)
    q_upper <- apply(roots, 1, stats::quantile, probs = (1 - conf_level/2))


    cis.lower <- theta + q_lower
    cis.upper <- theta + q_upper

    return(list(low_predicted = cis.lower, up_predicted = cis.upper, sd = sd_est, fitted_resid = fitted_resid,
        center_resid = center_resid))
}

#' @keywords internal
resids_func <- function(x, func_est, first_resid, bandwidth, bias_correct) {
  npts <- length(x)
  resid.resamp <- sample(first_resid, size = npts, replace = TRUE)
  newy <- func_est + resid.resamp
  newfit <- nwvector(x, newy, bandwidth)
  if (bias_correct == TRUE) {
    newfit <- biasnwvector(x, newy, newfit, bandwidth)
  }
  return(newfit)
}


#' Residual resampling bootstrap pointwise confidence bands for '\code{NOIS_fit}' objects
#'
#' @inheritParams pred_resid_BS_confint
#' @return A list with the following components.
#' \item{\code{up_predicted}}{The upper band.}
#' \item{\code{low_predicted}}{The lower band.}
#' @family NOIS confidence bands
#' @importFrom foreach %dopar%
#' @export
resid_BS_confint <- function(NOIS_fit, conf_level = 0.05, fit_type = "NOIS", bias_correct = T, parallel = F, B = 500) {
    x <- NOIS_fit$fit_df$x
    if (fit_type == "NOIS") {
        y <- NOIS_fit$fit_df$y_adj
        bandwidth <- NOIS_fit$CV$pool_h
        func_est <- NOIS_fit$fit_df$fit
        if (bias_correct == T) {
            func_est <- NOIS_fit$fit_df$bias_fit
        }
    } else if (fit_type == "regular") {
        y <- NOIS_fit$fit_df$y
        bandwidth <- NOIS_fit$CV$first_h
        func_est <- NOIS_fit$fit_df$nr_fit
        if (bias_correct == T) {
            func_est <- NOIS_fit$fit_df$nr_bias_fit
        }
    } else {
        stop("Must supply valid fit_type: NOIS or regular.")
    }



    resid <- y - func_est
    center_resid <- resid - mean(resid)


    if (parallel == T) {
        `%fun%` <- doRNG::`%dorng%`
        i <- NULL
        ret <- foreach::foreach(i = 1:B, .combine = cbind) %fun% {
            loop_ret <- resids_func(x, func_est, center_resid, bandwidth, bias_correct)
            return(loop_ret)
        }
        resids_boots <- ret
    } else {
        resids_boots <- replicate(B, resids_func(x, func_est, center_resid, bandwidth, bias_correct))
    }

    roots <- apply(resids_boots, 2, function(x) {
        func_est - x
    })

    q_lower <- apply(roots, 1, stats::quantile, probs = conf_level/2)
    q_upper <- apply(roots, 1, stats::quantile, probs = (1 - conf_level/2))


    cis.lower <- func_est + q_lower
    cis.upper <- func_est + q_upper
    return(list(up_predicted = cis.upper, low_predicted = cis.lower))
}


#' @keywords internal
NOIS_logelr_root <- function(yvals, hyp_theta, gkcalc, conf_level = 0.05, nwest_val = NULL, index = NULL,
    calib_type = "F") {
    if (is.null(nwest_val) & is.null(index)) {
        score_vec <- gkcalc * (yvals - hyp_theta)
    } else if (!is.null(nwest_val) & !is.null(index)) {
        score_vec <- gkcalc * (yvals - hyp_theta - nwest_val + nwest_val[index])
    } else {
        stop("Missing a parameter")
    }
    npts <- dim(score_vec)[1]
    elm <- emplik(score_vec)
    if (calib_type == "chisq") {
        thresh <- stats::qchisq(1 - conf_level, df = 1)
    } else if (calib_type == "F") {
        thresh <- stats::qf(1 - conf_level, 1, length(score_vec) - 1)
    } else {
        stop("Must supply a valid calibration type (F or chisq)")
    }
    funcoutput <- elm$logelr - thresh
    return(funcoutput)
}


EL_rootfun <- function(hyp_val, y, gkcalc, calib_type, conf_level, bias_correct, nwfit, i) {
  if (bias_correct == FALSE) {
    llr <- NOIS_logelr_root(y, hyp_val, gkcalc, calib_type = calib_type, conf_level = conf_level)
  } else {
    llr <- NOIS_logelr_root(y, hyp_val, gkcalc, calib_type = calib_type, conf_level = conf_level, nwest_val = nwfit,
                            index = i)
  }
  return(llr)
}



#' Empirical likelihood pointwise confidence bands for '\code{NOIS_fit}' objects
#'
#' Returns EL confidence bands.
#'
#' @param NOIS_fit A \code{NOIS_fit}.
#' @param conf_level The significance level.
#' @param fit_type The type of fit to use for confidence bands.
#' Valid types are \code{c('NOIS', 'regular')}, where regular is the non-robust fit..
#' @param bias_correct A logical indicating usage of bias correction.
#' @param parallel A logical indicating parallel computation. A backend must be registered first.
#' @param left The left summand for the root finding procedure.
#' @param right The right summand for the root finding procedure.
#' @param maxit The maximum number of iterations for each root finding procedure.
#' @param calib_type The distribution for calibrating Wilks' theorem. Valid types are \code{c('F', 'chisq')}.
#' @return A list with the following components.
#' \item{\code{up_predicted}}{The upper band.}
#' \item{\code{low_predicted}}{The lower band.}
#' \item{\code{time}}{Elapsed time.}
#' \item{\code{up_iter}}{Number of iterations for the upper band.}
#' \item{\code{low_iter}}{Number of iterations for the lower band.}
#' @family NOIS confidence bands
#' @export
EL_confint <- function(NOIS_fit, conf_level = 0.05, fit_type = "NOIS", bias_correct = T, parallel = F, calib_type = "F",
    left = 0, right = 20, maxit = 50) {
    x <- NOIS_fit$fit_df$x
    if (fit_type == "NOIS") {
        bandwidth <- NOIS_fit$CV$pool_h
        y <- NOIS_fit$fit_df$y_adj
        if (bias_correct == TRUE) {
            nwfit <- NOIS_fit$fit_df$fit
            theta <- NOIS_fit$fit_df$bias_fit
        } else {
            nwfit <- NULL
            theta <- NOIS_fit$fit_df$fit
        }
    } else if (fit_type == "regular") {
        bandwidth <- NOIS_fit$CV$first_h
        y <- NOIS_fit$fit_df$y
        if (bias_correct == TRUE) {
            nwfit <- NOIS_fit$fit_df$nr_fit
            theta <- NOIS_fit$fit_df$nr_bias_fit
        } else {
            nwfit <- NULL
            theta <- NOIS_fit$fit_df$nr_fit
        }
    } else {
        stop("Must supply valid fit_type: NOIS or regular.")
    }

    if (parallel == F) {
        `%fun%` <- foreach::`%do%`
    } else {
        `%fun%` <- foreach::`%dopar%`
    }
    i <- NULL

    ptm <- proc.time()
    loop_obj <- foreach::foreach(i = 1:length(x)) %fun% {
        gkcalc <- as.matrix(stats::dnorm(x - x[i], 0, bandwidth))
        # gkcalc <- as.matrix(gausskern(x - x[i], bandwidth))

        # EL_rootfun <- function(hyp_val) {
        #   if (bias_correct == FALSE) {
        #     llr <- NOIS_logelr_root(y, hyp_val, gkcalc, calib_type = calib_type, conf_level = conf_level)
        #   } else {
        #     llr <- NOIS_logelr_root(y, hyp_val, gkcalc, calib_type = calib_type, conf_level = conf_level, nwest_val = nwfit,
        #                             index = i)
        #   }
        #   return(llr)
        # }


        mletheta <- theta[i]
        left <- left
        right <- right
        nwupsoln <- stats::uniroot(EL_rootfun, interval = c(mletheta + left, mletheta + right), y = y, gkcalc = gkcalc, calib_type = calib_type, conf_level = conf_level, bias_correct = bias_correct, nwfit = nwfit, i = i, check.conv = TRUE, maxiter = maxit)
        # nwupsoln <- stats::uniroot(EL_rootfun, c(mletheta + left, mletheta + right), check.conv = TRUE, maxiter = maxit)
        up_val <- nwupsoln$root
        nwlowsoln <- stats::uniroot(EL_rootfun, interval = c(mletheta - left, mletheta - right), y = y, gkcalc = gkcalc, calib_type = calib_type, conf_level = conf_level, bias_correct = bias_correct, nwfit = nwfit, i = i, check.conv = TRUE, maxiter = maxit)
        # nwlowsoln <- stats::uniroot(EL_rootfun, c(mletheta - left, mletheta - right), check.conv = TRUE, maxiter = maxit)
        low_val <- nwlowsoln$root
        upiter <- nwupsoln$iter
        lowiter <- nwlowsoln$iter
        return(list(up_val = up_val, low_val = low_val, up_iter = upiter, low_iter = lowiter))
    }
    looptm <- proc.time() - ptm

    up_val <- purrr::map_dbl(loop_obj, "up_val")
    low_val <- purrr::map_dbl(loop_obj, "low_val")

    up_iter <- purrr::map_int(loop_obj, "up_iter")
    low_iter <- purrr::map_int(loop_obj, "low_iter")

    return(list(up_predicted = up_val, low_predicted = low_val, time = looptm, up_iter = up_iter, low_iter = low_iter))
}

#' Pointwise confidence bands for '\code{NOIS_fit}' objects
#'
#' Returns pointwise confidence bands.
#'
#' This method calls three types of confidence intervals -
#' predictive residuals bootstrap, residual resampling bootstrap or empirical likelihood.
#'
#' @param NOIS_fit A \code{NOIS_fit}.
#' @param conf_type The type of confidence interval construction. Valid types are \code{c('pred_BS', 'resid_BS', 'EL')}.
#' @param ... Additional parameters passed to
#' \code{\link{pred_resid_BS_confint}}, \code{\link{resid_BS_confint}}, \code{\link{EL_confint}}.
#' @family NOIS confidence bands
#' @export
NOIS_confint <- function(NOIS_fit, conf_type = "pred_BS", ...) {
    if (class(NOIS_fit) != "NOIS_fit") {
        stop("Input must be a NOIS_fit")
    }
    switch(conf_type, pred_BS = pred_resid_BS_confint(NOIS_fit, ...), resid_BS = resid_BS_confint(NOIS_fit, ...), EL = EL_confint(NOIS_fit, ...))
}


