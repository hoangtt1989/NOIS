NOIS
================
Hoang Tran
5/9/2017

Non-Parametric Outlier Identification and Smoothing.
----------------------------------------------------

This package implements outlier identification and smoothing for non-linear and time series data.

Installation
------------

``` r
install.packages('devtools')
devtools::install_github('hoangtt1989/NOIS')
```

Simulated Data Example
----------------------

We generate a random sine curve and perturb it with outliers.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(NOIS)
set.seed(123)
npts <- 200
nout <- floor(.1*npts)

xt <- seq(from=0, to=2*pi, length.out=npts)
gaussnoise <- rnorm(npts)

outliers <- sample(floor(npts/2):npts, size=nout)
randpts <- runif(nout, min=5, max=7)
yt <- sin(xt) + gaussnoise
yt[outliers] <- yt[outliers] + randpts
xt_outliers = xt[outliers]

orig_func <- sin(xt)
data <- data.frame(x = xt, y = yt, true_outlier = xt %in% xt_outliers)

ggplot(data) + geom_point(aes(x, y, color = true_outlier))
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

We use NOIS to detect outliers and estimate the underlying curve. Inputs to this function include a `data.frame` as well as strings specifying which columns contain "x" and "y" values. See the help documentation for `NOIS_fit` for additional information about function inputs.

``` r
sine_fit <- NOIS_fit(data, x = 'x', y = 'y', CV_method = 'LOOCV', pool_q = nout)
sine_fit
```

    ## Number of detected outliers = 20 
    ## Number of observations = 200 
    ## Convergence = TRUE 
    ## MSE = 0.7681321 
    ## Bias corrected MSE = 0.7573429 
    ## First optimal bandwidth = 1.767722 
    ## Pooled optimal bandwidth = 0.4234177

We plot the estimated curve and mark the detected outliers. NOIS correctly detects all of the outliers.

``` r
pt <- outlier_plot(sine_fit)
pt
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)

We use the predictive residuals bootstrap to construct confidence bands.

``` r
sine_BS <- NOIS_confint(sine_fit, conf_type = 'pred_BS', B = 500, conf_level = .05)
sine_df <- NOIS_df(sine_fit) %>%
  mutate(low_conf = sine_BS$low_predicted,
         up_conf = sine_BS$up_predicted)
pt2 <- ggplot(sine_df) + geom_point(aes(x, y, color = outlier)) + geom_line(aes(x, bias_fit)) + geom_ribbon(aes(x = x, ymin = low_conf, ymax = up_conf), alpha = .5)
pt2
```

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)
