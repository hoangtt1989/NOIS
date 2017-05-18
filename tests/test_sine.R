library(NOIS)

set.seed(123)
npts <- 1000
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

fit <- NOIS_fit(data, CV_method = 'MCV')
fit$conv$time
fit$local_q
outlier_plot(fit)
length(intersect(fit$pool_outlier, outliers))
# EL <- NOIS_confint(fit, conf_type = "EL", right = 5)
EL$up_predicted > EL$low_predicted
