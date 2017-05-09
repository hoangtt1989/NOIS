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

We use NOIS to detect outliers and estimate the underlying curve.

``` r
sine_fit <- NOIS_fit(data, CV_method = 'LOOCV', pool_q = nout)
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
outlier_plot(sine_fit)
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)
