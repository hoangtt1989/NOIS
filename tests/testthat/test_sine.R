library(NOIS)
context('Test on a sine curve')

set.seed(123)
npts <- 100
nout <- floor(.1*npts)

xt <- seq(from=0, to=2*pi, length.out=npts)
gaussnoise <- rnorm(npts)

outliers <- sample(floor(npts/2):npts, size=nout)
randpts <- runif(nout, min=5, max=7)
yt <- sin(xt) + gaussnoise
yt[outliers] <- yt[outliers] + randpts
xt_outliers = xt[outliers]

orig_func <- sin(xt)
data <- data.frame(x=xt, y=yt)

test_that('Detected outliers match true outliers', {
  expect_equal(length(intersect(NOIS_fit(data, CV_method = 'MCV')$pool_outlier, outliers)), nout)
  expect_equal(length(intersect(NOIS_fit(data, CV_method = 'PCV')$pool_outlier, outliers)), nout)
  expect_equal(length(intersect(NOIS_fit(data, CV_method = 'LOOCV')$pool_outlier, outliers)), nout)
})

test_that('BIC tuning', {
  expect_equal(BIC_tuner(data, CV_method = 'MCV')$min_q, 13)
})
