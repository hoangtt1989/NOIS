#' #' @keywords internal
#' NOIS_inner <- function(xx_inp, xx, yy, nn, first_h, local_q, tol, maxit) {
#'
#'   xx_inner <- xx_inp
#'   kernlist <- stats::dnorm(xx_inner - xx, 0, first_h)
#'   nz_ind <- which(kernlist != 0 & kernlist >= 1e-20)
#'   kern_nz <- rep(0, nn)
#'   kern_nz[nz_ind] <- kernlist[nz_ind]
#'   kern_nzsqrt <- rep(0, nn)
#'   kern_nzsqrt[nz_ind] <- sqrt(kern_nz)[nz_ind]
#'   kern_nzsqrtinv <- rep(0, nn)
#'   kern_nzsqrtinv[nz_ind] <- 1/kern_nzsqrt[nz_ind]
#'   qq_inner <- qdet(local_q, length(nz_ind))
#'
#'   gamma_inner <- rep(0, nn)
#'
#'   for (ii in 1:maxit) {
#'     yy_adj <- yy - gamma_inner
#'     local_inner <- nwestimator(xx_inner, xx, yy_adj, first_h)
#'     rr <- kern_nzsqrt * (yy - local_inner)
#'     gamma_next <- kern_nzsqrtinv * quantile_thresh(rr, qq_inner)
#'     cond_inner <- max(abs(gamma_next - gamma_inner))
#'     gamma_inner <- gamma_next
#'
#'     if (max(abs(cond_inner)) <= tol) {
#'       converge_inner <- T
#'       break
#'     } else {
#'       converge_inner = F
#'     }
#'   }
#'   if(ii == maxit) {
#'     warning('One of the points did not converge. Check status in conv.')
#'   }
#'   return(list(local_fit = local_inner, gamma_curr = gamma_inner, qq_inner = qq_inner, converged = converge_inner,
#'               cond_check = cond_inner, iter = ii))
#' }
#
# NOIS_fit_lapply <- function(data, x = "x", y = "y", CV_method = "LOOCV", first_h = NULL, pool_h = NULL, local_q = 0.1,
#                             pool_q = 0.1, tol = 1e-07, maxit = 200, ...) {
#
#   yy <- data[[y]]
#   xx <- data[[x]]
#
#   # checking inputs
#   if (!is.data.frame(data)) {
#     stop("data must be a data.frame")
#   }
#   if (any(!(c(x, y) %in% colnames(data)))) {
#     stop("x and y must be column names in data")
#   }
#   if (local_q <= 0 | local_q >= 1) {
#     stop("local_q must be between 0 and 1")
#   }
#   if (pool_q <= 0) {
#     stop("pool_q must be greater than 0")
#   }
#   if (!all(is.finite(xx)) | !all(is.finite(yy))) {
#     stop("Bad values in either x or y")
#   }
#
#   # initializing vectors, parameters
#   nn <- nrow(data)
#   # cond_check <- rep(0, nn)
#   # gamma_curr <- matrix(0, nrow = nn, ncol = nn)
#   # local_fit <- rep(0, nn)
#   # qq <- rep(0, nn)
#   # converged <- rep(0, nn)
#   # iter <- rep(0, nn)
#   # cond_check <- rep(1, nn)
#
#   if (CV_method %in% c("LOOCV", "MCV", "PCV")) {
#     first_CV <- switch(CV_method, LOOCV = LOOCV_grid(xx, yy, ...), MCV = MCV_grid(xx, yy, ...), PCV = PCV_grid(xx, yy, ...))
#     first_h <- first_CV$min_h
#     first_hgrid <- first_CV$hgrid
#     if (first_h == first_hgrid[1] | first_h == first_hgrid[length(first_hgrid)]) {
#       warning("Optimal first bandwidth is at the end of the grid")
#     }
#   } else if (CV_method == "none") {
#     if (is.null(first_h)) {
#       warning("No CV method supplied or bandwidth supplied. Using theoretically optimal bandwidth for first h.")
#       first_h <- nrow(data)^(-1/5)
#     }
#     first_hgrid <- NULL
#   } else {
#     stop("Supply valid CV_method - LOOCV, MCV, PCV or none")
#   }
#
#   nw_ests <- nwvector(xx, yy, first_h)
#   bias_nw_ests <- biasnwvector(xx, yy, nw_ests, first_h)
#
#   ptm <- proc.time()
#   inner_fit <- lapply(xx, NOIS_inner, xx = xx, yy = yy, nn = nn, first_h = first_h, local_q = local_q, tol = tol, maxit = maxit)
#   local_fit <- purrr::map_dbl(inner_fit, 'local_fit')
#   qq_inner <- purrr::map_dbl(inner_fit, 'qq_inner')
#   gamma_curr <- do.call(cbind, lapply(inner_fit, function(x){x$gamma_curr}))
#   converged <- purrr::map_lgl(inner_fit, 'converged')
#   cond_check <- purrr::map_dbl(inner_fit, 'cond_check')
#   iter <- purrr::map_int(inner_fit, 'iter')
#   # for (jj in 1:nn) {
#   #
#   #     xx_inner <- xx[jj]
#   #
#   #     # kernlist <- dnorm(xx_inner - xx, 0, first_h)
#   #     # nz_ind <- which(kernlist != 0 & kernlist >= 1e-20)
#   #     # kern_nz <- kernlist[nz_ind]
#   #     # kern_nzsqrt <- sqrt(kern_nz)
#   #     # kern_nzsqrtinv <- 1/kern_nzsqrt
#   #     kernlist <- dnorm(xx_inner - xx, 0, first_h)
#   #     nz_ind <- which(kernlist != 0 & kernlist >= 1e-20)
#   #     kern_nz <- rep(0, nn)
#   #     kern_nz[nz_ind] <- kernlist[nz_ind]
#   #     kern_nzsqrt <- rep(0, nn)
#   #     kern_nzsqrt[nz_ind] <- sqrt(kern_nz)[nz_ind]
#   #     kern_nzsqrtinv <- rep(0, nn)
#   #     kern_nzsqrtinv[nz_ind] <- 1/kern_nzsqrt[nz_ind]
#   #     qq[jj] <- qdet(local_q, kern_nz)
#   #
#   #     gamma_inner <- rep(0, nn)
#   #     qq_inner <- qq[jj]
#   #
#   #     for (ii in 1:maxit) {
#   #         # gamma_next <- rep(0, nn)
#   #         yy_adj <- yy - gamma_inner
#   #         local_inner <- nwestimator(xx_inner, xx, yy_adj, first_h)
#   #         rr <- kern_nzsqrt * (yy - local_inner)
#   #         # rr <- kern_nzsqrt * ((yy - local_inner))[nz_ind]
#   #         gamma_next <- kern_nzsqrtinv * quantile_thresh(rr, qq_inner)
#   #         # gamma_next[nz_ind] <- kern_nzsqrtinv * quantile_thresh(rr, qq_inner)
#   #         cond_inner <- max(abs(gamma_next - gamma_inner))
#   #         gamma_inner <- gamma_next
#   #
#   #         if (max(abs(cond_inner)) <= tol) {
#   #             converge_inner <- T
#   #             break
#   #         } else {
#   #             converge_inner = F
#   #         }
#   #     }
#   #
#   #     local_fit[jj] <- local_inner
#   #     gamma_curr[, jj] <- gamma_inner
#   #     converged[jj] <- converge_inner
#   #     cond_check[jj] <- cond_inner
#   #     iter[jj] <- ii
#   #
#   #     if (ii == maxit) {
#   #         warning(paste("Model did not converge at j =", jj))
#   #     }
#   # }
#   etm <- proc.time() - ptm
#
#   # get max
#   gam_max <- apply(gamma_curr, 1, function(x) {
#     x[which.max(abs(x))]
#   })
#
#   # index for gamma vectors
#   if (pool_q < 1) {
#     pool_q <- min(floor(pool_q * nn), length(which(gam_max != 0)))
#   } else {
#     # so that we can supply pool_q as an integer if needed
#     pool_q <- min(pool_q, length(which(gam_max != 0)))
#   }
#   gam_val <- quantile_thresh(gam_max, pool_q)
#   gam_ind <- which(gam_val != 0)
#
#   # pooled y adj
#   pool_y_adj <- yy - gam_val
#   # pooled cv
#   if (CV_method %in% c("LOOCV", "MCV", "PCV")) {
#     pool_CV <- switch(CV_method, LOOCV = LOOCV_grid(xx, pool_y_adj, ...), MCV = MCV_grid(xx, pool_y_adj, ...), PCV = PCV_grid(xx, pool_y_adj, ...))
#     pool_h <- pool_CV$min_h
#     pool_hgrid <- pool_CV$hgrid
#     if (pool_h == pool_hgrid[1] | pool_h == pool_hgrid[length(pool_hgrid)]) {
#       warning("Optimal pooled bandwidth is at the end of the grid")
#     }
#   } else if (CV_method == "none") {
#     if (is.null(first_h)) {
#       warning("No CV method supplied or bandwidth supplied. Using theoretically optimal bandwidth for first h.")
#       pool_h <- nrow(data)^(-1/5)
#     }
#     pool_hgrid <- NULL
#   }
#
#   # pooled fit
#   pool_fit <- nwvector(xx, pool_y_adj, pool_h)
#   bias_pool_fit <- biasnwvector(xx, pool_y_adj, pool_fit, pool_h)
#   pool_nonout <- setdiff(1:nn, gam_ind)
#
#   # CV
#   CV <- list(first_h = first_h, pool_h = pool_h, first_hgrid = first_hgrid, pool_hgrid = pool_hgrid)
#   # convergence
#   conv <- list(iter = iter, time = etm, converged = converged, cond_check = cond_check)
#
#   model_output <- list(local_fit = local_fit, pool_fit = pool_fit, bias_pool_fit = bias_pool_fit, first_fit = nw_ests,
#                        bias_first_fit = bias_nw_ests, local_gamma = gamma_curr, pool_gamma = gam_val, pool_outlier = gam_ind,
#                        local_q = qq_inner, pool_q = pool_q, pool_nonout = pool_nonout, x = xx, y_adj = pool_y_adj, y = yy, CV = CV, conv = conv)
#   class(model_output) <- "NOIS_fit"
#   model_output
# }


###loop inside NOIS_fit
# for (jj in 1:nn) {
#
#     xx_inner <- xx[jj]
#     kernlist <- stats::dnorm(xx_inner - xx, 0, first_h)
#     nz_ind <- which(kernlist != 0 & kernlist >= 1e-20)
#     kern_nz <- rep(0, nn)
#     kern_nz[nz_ind] <- kernlist[nz_ind]
#     kern_nzsqrt <- rep(0, nn)
#     kern_nzsqrt[nz_ind] <- sqrt(kern_nz[nz_ind])
#     kern_nzsqrtinv <- rep(0, nn)
#     kern_nzsqrtinv[nz_ind] <- 1/(kern_nzsqrt[nz_ind])
#     qq[jj] <- qdet(local_q, length(nz_ind))
#
#     gamma_inner <- rep(0, nn)
#     qq_inner <- qq[jj]
#
#     for (ii in 1:maxit) {
#         yy_adj <- yy - gamma_inner
#         local_inner <- nwestimator(xx_inner, xx, yy_adj, first_h)
#         rr <- kern_nzsqrt * (yy - local_inner)
#         gamma_next <- kern_nzsqrtinv * quantile_thresh(rr, qq_inner)
#         cond_inner <- max(abs(gamma_next - gamma_inner))
#         gamma_inner <- gamma_next
#
#         if (max(abs(cond_inner)) <= tol) {
#             converge_inner <- T
#             break
#         } else {
#             converge_inner = F
#         }
#     }
#
#     local_fit[jj] <- local_inner
#     gamma_curr[, jj] <- gamma_inner
#     converged[jj] <- converge_inner
#     cond_check[jj] <- cond_inner
#     iter[jj] <- ii
#
#     if (ii == maxit) {
#         warning(paste("Model did not converge at j =", jj))
#     }
# }



#' #' Construct a \code{data_frame} from a \code{NOIS_fit}.
#' #'
#' #' Construct a \code{data_frame} using output from a \code{NOIS_fit}.
#' #'
#' #' @param NOIS_fit A \code{NOIS_fit}
#' #' @return A \code{data_frame} with the following columns.
#' #' \item{\code{x}}{'x' values.}
#' #' \item{\code{y}}{'y' values.}
#' #' \item{\code{y_adj}}{Pooled adjusted 'y' values.}
#' #' \item{\code{fit}}{Pooled NOIS fit without bias correction.}
#' #' \item{\code{bias}}{Pooled NOIS fit with bias correction.}
#' #' \item{\code{outlier}}{A logical indicating whether this point is an outlier.}
#' #' @export
#' NOIS_df <- function(NOIS_fit) {
#'   if (class(NOIS_fit) != "NOIS_fit") {
#'     stop("Input must be a NOIS_fit")
#'   }
#'   outlier_index <- rep(FALSE, length(NOIS_fit$x))
#'   outlier_index[NOIS_fit$pool_outlier] <- TRUE
#'   df <- with(NOIS_fit, tibble::data_frame(x = x, y = y, y_adj = y_adj, fit = pool_fit, bias_fit = bias_pool_fit,
#'                                           outlier = outlier_index))
#' }

