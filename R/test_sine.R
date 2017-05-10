# library(NOIS) set.seed(123) npts <- 100 nout <- floor(.1*npts) xt <- seq(from=0, to=2*pi, length.out=npts)
# gaussnoise <- rnorm(npts) outliers <- sample(floor(npts/2):npts, size=nout) randpts <- runif(nout, min=5,
# max=7) yt <- sin(xt) + gaussnoise yt[outliers] <- yt[outliers] + randpts xt_outliers = xt[outliers] orig_func
# <- sin(xt) data <- data.frame(x=xt, y=yt) test_fit <- NOIS_fit(data, CV_method = 'LOOCV')
# test_fit$pool_outlier outlier_plot(test_fit) length(intersect(test_fit$pool_outlier, outliers)) test_fit <-
# NOIS_fit(data, CV_method = 'MCV') outlier_plot(test_fit) length(intersect(test_fit$pool_outlier, outliers))
# test_fit <- NOIS_fit(data, CV_method = 'PCV') outlier_plot(test_fit) length(intersect(test_fit$pool_outlier,
# outliers)) BIC_test <- BIC_tuner(data, CV_method = 'MCV') BIC_test$min_q BIC_test$q_tst
