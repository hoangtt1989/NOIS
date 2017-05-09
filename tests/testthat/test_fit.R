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


sine_fit <- NOIS_fit(data, 'x', 'y', pool_q = .1, CV_method = 'LOOCV')
class(sine_fit)
sine_fit
length(intersect(sine_fit$pool_outlier, outliers))

BIC(sine_fit)

sine_df <- NOIS_df(sine_fit)

sine_plot <- outlier_plot(sine_fit)
sine_plot

sine_BS <- confint(sine_fit, conf_type = 'resid_BS', B = 100, parallel = F)
all(sine_BS$up_predicted > sine_BS$low_predicted)

cl <- parallel::makeForkCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
sine_BS <- confint(sine_fit, conf_type = 'resid_BS', B = 100, parallel = T)
parallel::stopCluster(cl)
all(sine_BS$up_predicted > sine_BS$low_predicted)


cl <- parallel::makeForkCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
sine_BS <- confint(sine_fit, conf_type = 'pred_BS', B = 100, parallel = T)
parallel::stopCluster(cl)
all(sine_BS$up_predicted > sine_BS$low_predicted)

cl <- parallel::makeForkCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
sine_EL <- confint(sine_fit, conf_type = 'EL', right = 5, parallel = T)
parallel::stopCluster(cl)
sine_EL$time
all(sine_EL$up_predicted > sine_EL$low_predicted)


# cl <- parallel::makeForkCluster(parallel::detectCores())
# doParallel::registerDoParallel(cl)
# sine_EL <- EL_confint(sine_fit, right = 10, parallel = T)
# parallel::stopCluster(cl)
# sine_EL$time
#
# sine_EL <- EL_confint(sine_fit, right = 10, parallel = F)
# sine_EL$time
