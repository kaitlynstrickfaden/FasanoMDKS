
#' Calculate D-values for a Two-dimensional Kolmogorov-Smirnov Test
#'
#' Find D-values of observations for a two-dimensional Kolmogorov-Smirnov test as defined by [Fasano and Franceschini 1987](https://academic.oup.com/mnras/article/225/1/155/1007281).
#'
#' @importFrom dplyr arrange group_by mutate summarize
#' @param xcol a vector of x values
#' @param ycol a vector of y values
#' @return The D-values for each pair of observations as a vector.
#' @export


## A function to determine D-values
Fasano_2D_dvals <- function(xcol, ycol) {

  if (length(xcol) != length(ycol)) {
    stop("xcol and ycol must be the same length")
  }

  n <- length(xcol)  # sample size

  # store number of values less/greater than a given x/y value
  xl = sapply(xcol, function(x) {length(xcol[xcol < x])}) # all x-values less than the X
  xg = sapply(xcol, function(x) {length(xcol[xcol > x])})# all x-values greater than the X
  yl = sapply(ycol, function(y) {length(ycol[ycol < y])}) # all y-values less than the Y
  yg = sapply(ycol, function(y) {length(ycol[ycol > y])})# all y-values greater than the Y

  # store number of values for a given quadrant
  ll = map2_dbl(xcol,ycol, function(x,y) {length(xcol[xcol < x & ycol < y])}) # <x, <y
  lg = map2_dbl(xcol,ycol, function(x,y) {length(xcol[xcol < x & ycol > y])}) # <x, >y
  gl = map2_dbl(xcol,ycol, function(x,y) {length(xcol[xcol > x & ycol < y])}) # >x, <y
  gg = map2_dbl(xcol,ycol, function(x,y) {length(xcol[xcol > x & ycol > y])}) # >x, >y

  # the EXPECTED proportion of observations in a given quadrant:
  # We'd 'expect' it to be a product of the proportion of observations on the left/right side of x, and top/bottom side of Y
  expected1 = map2_dbl(xl,yl, function(x,y) {(x/n) * (y/n)}) # lower left
  expected2 = map2_dbl(xl,yg, function(x,y) {(x/n) * (y/n)}) # upper left
  expected3 = map2_dbl(xg,yl, function(x,y) {(x/n) * (y/n)}) # lower right
  expected4 = map2_dbl(xg,yg, function(x,y) {(x/n) * (y/n)}) # upper right

  # the OBSERVED proportion of observations in a given quadrant
  observed1 = sapply(ll, function(x) {x/n}) # lower left
  observed2 = sapply(lg, function(x) {x/n}) # upper left
  observed3 = sapply(gl, function(x) {x/n}) # lower right
  observed4 = sapply(gg, function(x) {x/n}) # upper right

  # compute differences in observed vs. expected proportions of observations within each quadrant
  d1 = map2_dbl(observed1, expected1, function(x,y) {abs(x - y)}) # lower left
  d2 = map2_dbl(observed2, expected2, function(x,y) {abs(x - y)}) # upper left
  d3 = map2_dbl(observed3, expected3, function(x,y) {abs(x - y)}) # lower right
  d4 = map2_dbl(observed4, expected4, function(x,y) {abs(x - y)}) # upper right

  # D-value: the maximum difference among all four quadrants
  maxdiff = pmap_dbl(list(d1, d2, d3, d4), max)

  return(maxdiff) # Function output

} # end quadrants function
