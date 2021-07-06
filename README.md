# FasanoMDKS

This package contains functions for performing multi-dimensional Kolmogorov-Smirnov (MDKS) tests as defined by [Fasano and Franceschini 1987](https://academic.oup.com/mnras/article/225/1/155/1007281). For now, this package only performs two-dimensional tests. Three-dimensional tests may be added at a later time.

Before use of the package, please install `dplyr`, `progress`, and `purrr`.

<br>

To download the package, run the following:

```
devtools::install_github("kaitlynstrickfaden/FasanoMDKS")
```

<br>

To use the package, simply input a data frame with an "x" and "y" column of paired observations you want to analyze. You can also set the number of randomizations to perform and the desired alpha level. A progress bar will keep track of progress through the randomizations.

```
xy <- data.frame(x = sample(c(1:100),50), y = sample(c(1:100),50))
Fasano_2DKS(xy, rands = 5000, alpha = 0.05)
```

The `Fasano_2DKS` function will output a data frame of summary values and statistics containing the maximum D-value, the x and y values at the maximum D-value, the minimum and maximum x and y values for the statistically-significant observations, and the p-value.

