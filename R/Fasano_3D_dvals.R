
#' Calculate D-values for a Three-dimensional Kolmogorov-Smirnov Test
#'
#' Find D-values of observations for a three-dimensional Kolmogorov-Smirnov test as defined by Fasano and Franceschini (1987).
#'
#' @importFrom purrr map2_dbl pmap_dbl
#' @param xcol a vector of x values
#' @param ycol a vector of y values
#' @param zcol a vector of z values
#' @return The D-values for each set of observations as a vector
#' @examples
#' xcol <- rnorm(10)
#' ycol <- rnorm(10)
#' zcol <- rnorm(10)
#'
#' Fasano_3D_dvals(xcol, ycol, zcol)
#'
#' \dontrun{
#' xcol <- rnorm(10)
#' ycol <- rnorm(10)
#' zcol <- rnorm(5)
#'
#' Fasano_3D_dvals(xcol, ycol, zcol)
#' }
#' @export


## A function to determine D-values
Fasano_3D_dvals <- function(xcol, ycol, zcol) {

  if (length(xcol) != length(ycol) |
      length(xcol) != length(zcol) |
      length(ycol) != length(zcol) ) {
    stop("xcol, ycol, and zcol must be the same length")
  }

  n <- length(xcol)  # sample size

  # store number of values less/greater than a given x/y value
  xl <- sapply(xcol, function(x) {length(xcol[xcol < x])}) # all x-values less than the X
  xg <- sapply(xcol, function(x) {length(xcol[xcol > x])}) # all x-values greater than the X
  yl <- sapply(ycol, function(y) {length(ycol[ycol < y])}) # all y-values less than the Y
  yg <- sapply(ycol, function(y) {length(ycol[ycol > y])}) # all y-values greater than the Y
  zl <- sapply(zcol, function(z) {length(zcol[zcol < z])}) # all z-values less than the Z
  zg <- sapply(zcol, function(z) {length(zcol[zcol > z])}) # all z-values greater than the Z

  # store number of values for a given quadrant
  lll <- pmap_dbl(list(xcol,ycol,zcol), # <x, <y, <z
                 function(x,y,z) {length(xcol[xcol < x & ycol < y & zcol < z])})
  lgl <- pmap_dbl(list(xcol,ycol,zcol), # <x, >y, <z
                 function(x,y,z) {length(xcol[xcol < x & ycol > y & zcol < z])})
  gll <- pmap_dbl(list(xcol,ycol,zcol), # >x, <y, <z
                 function(x,y,z) {length(xcol[xcol > x & ycol < y & zcol < z])})
  ggl <- pmap_dbl(list(xcol,ycol,zcol), # >x, >y, <z
                 function(x,y,z) {length(xcol[xcol > x & ycol > y & zcol < z])})
  llg <- pmap_dbl(list(xcol,ycol,zcol), # <x, <y, >z
                 function(x,y,z) {length(xcol[xcol < x & ycol < y & zcol > z])})
  lgg <- pmap_dbl(list(xcol,ycol,zcol), # <x, >y, >z
                 function(x,y,z) {length(xcol[xcol < x & ycol > y & zcol > z])})
  glg <- pmap_dbl(list(xcol,ycol,zcol), # >x, <y, >z
                 function(x,y,z) {length(xcol[xcol > x & ycol < y & zcol > z])})
  ggg <- pmap_dbl(list(xcol,ycol,zcol), # >x, >y, >z
                 function(x,y,z) {length(xcol[xcol > x & ycol > y & zcol > z])})

  # the EXPECTED proportion of observations in a given quadrant:
  # We'd 'expect' it to be a product of the proportion of observations on the left/right side of x, and top/bottom side of Y
  expected1 <- pmap_dbl(list(xl,yl,zl), function(x,y,z) {(x/n) * (y/n) * (z/n)}) # bottom lower left
  expected2 <- pmap_dbl(list(xl,yg,zl), function(x,y,z) {(x/n) * (y/n) * (z/n)}) # bottom upper left
  expected3 <- pmap_dbl(list(xg,yl,zl), function(x,y,z) {(x/n) * (y/n) * (z/n)}) # bottom lower right
  expected4 <- pmap_dbl(list(xg,yg,zl), function(x,y,z) {(x/n) * (y/n) * (z/n)}) # bottom upper right
  expected5 <- pmap_dbl(list(xl,yl,zg), function(x,y,z) {(x/n) * (y/n) * (z/n)}) # top lower left
  expected6 <- pmap_dbl(list(xl,yg,zg), function(x,y,z) {(x/n) * (y/n) * (z/n)}) # top upper left
  expected7 <- pmap_dbl(list(xg,yl,zg), function(x,y,z) {(x/n) * (y/n) * (z/n)}) # top lower right
  expected8 <- pmap_dbl(list(xg,yg,zg), function(x,y,z) {(x/n) * (y/n) * (z/n)}) # top upper right

  # the OBSERVED proportion of observations in a given quadrant
  observed1 <- sapply(lll, function(x) {x/n}) # bottom lower left
  observed2 <- sapply(lgl, function(x) {x/n}) # bottom upper left
  observed3 <- sapply(gll, function(x) {x/n}) # bottom lower right
  observed4 <- sapply(ggl, function(x) {x/n}) # bottom upper right
  observed5 <- sapply(llg, function(x) {x/n}) # top lower left
  observed6 <- sapply(lgg, function(x) {x/n}) # top upper left
  observed7 <- sapply(glg, function(x) {x/n}) # top lower right
  observed8 <- sapply(ggg, function(x) {x/n}) # top upper right

  # compute differences in observed vs. expected proportions of observations within each quadrant
  d1 <- map2_dbl(observed1, expected1, function(x,y) {abs(x - y)}) # bottom lower left
  d2 <- map2_dbl(observed2, expected2, function(x,y) {abs(x - y)}) # bottom upper left
  d3 <- map2_dbl(observed3, expected3, function(x,y) {abs(x - y)}) # bottom lower right
  d4 <- map2_dbl(observed4, expected4, function(x,y) {abs(x - y)}) # bottom upper right
  d5 <- map2_dbl(observed5, expected5, function(x,y) {abs(x - y)}) # top lower left
  d6 <- map2_dbl(observed6, expected6, function(x,y) {abs(x - y)}) # top upper left
  d7 <- map2_dbl(observed7, expected7, function(x,y) {abs(x - y)}) # top lower right
  d8 <- map2_dbl(observed8, expected8, function(x,y) {abs(x - y)}) # top upper right

  # D-value: the maximum difference among all eight quadrants
  maxdiff <- pmap_dbl(list(d1, d2, d3, d4, d5, d6, d7, d8), max)

  return(maxdiff) # Function output

} # end function
