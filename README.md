# FasanoMDKS

This package contains functions for performing multi-dimensional Kolmogorov-Smirnov (MDKS) tests as defined by [Fasano and Franceschini 1987](https://academic.oup.com/mnras/article/225/1/155/1007281). This package can perform both two-dimensional (2D) tests and three-dimensional (3D) tests.

Before use of the package, please install `dplyr`, `progress`, and `purrr`.

<br>

To download the package, run the following:

```
devtools::install_github("kaitlynstrickfaden/FasanoMDKS")
```

<br>

To use the package, simply input "x" and "y" vectors (and a "z" vector for a 3D test) of data you want to analyze. You can also set the number of randomizations to perform and the desired alpha level. A progress bar will keep track of progress through the randomizations.

```
xcol <- sample(c(1:100),50)
ycol <- sample(c(1:100),50)
zcol <- sample(c(1:100),50)

Fasano_2DKS(xcol, ycol, rands = 5000, alpha = 0.05)
Fasano_3DKS(xcol, ycol, zcol, rands = 5000, alpha = 0.05)
```

The `Fasano_2DKS` and `Fasano_3DKS` functions will output a data frame of summary values and statistics containing the maximum D-value; the x, y, z, and p values at the maximum D-value; and the minimum and maximum x, y, and z values for the statistically-significant observations.

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/kaitlynstrickfaden/FasanoMDKS/workflows/R-CMD-check/badge.svg)](https://github.com/kaitlynstrickfaden/FasanoMDKS/actions)
  <!-- badges: end -->

