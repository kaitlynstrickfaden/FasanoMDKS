
#' Calculate D-values for a Three-dimensional Kolmogorov-Smirnov Test
#'
#' Find D-values of observations for a three-dimensional Kolmogorov-Smirnov test as defined by [Fasano and Franceschini 1987](https://academic.oup.com/mnras/article/225/1/155/1007281).
#'
#' @importFrom dplyr arrange group_by mutate summarize
#' @param xyz The D-values for each observation as a vector.
#' @return The original data frame with a new "maxdiff" column containing D-values.
#' @export



## A function to determine D-values
Fasano_3D_dvals <- function(xyz) {

  if ("x" %in% colnames(xyz) == FALSE |
      "y" %in% colnames(xyz) == FALSE |
      "z" %in% colnames(xyz) == FALSE) {
    stop("xyz must contain an 'x', 'y', and 'z' column")
  }

  n <- nrow(xyz)  # sample size
  xth <- xyz$x    # isolate x column
  yth <- xyz$y    # isolate y column
  zth <- xyz$z    # isolate z column

  # manipulate input data to find D values
  xyznew <- xyz %>%

    arrange(x) %>%

    mutate(i = 1:n) %>% # add row index

    group_by(i) %>% # group by row index (need to summarize each row individually rather than the whole data frame)

    summarize(

      # store of values less/greater than a given x/y value
      xl = length(xth[xth < .data$x]), # all x-values less than the X
      xg = length(xth[xth > .data$x]), # all x-values greater than the X
      yl = length(yth[yth < .data$y]), # all y-values less than the Y
      yg = length(yth[yth > .data$y]), # all x-values greater than the Y
      zl = length(zth[zth < .data$z]), # all z-values less than the Z
      zg = length(zth[zth > .data$z]), # all z-values greater than the Z

      # store of values for a given quadrant
      lll = length(xth[xth < .data$x & yth < .data$y & zth < .data$z]), # all obs with lower x's, lower y's, AND lower z's (bottom lower left quadrant)
      lgl = length(xth[xth < .data$x & yth > .data$y & zth < .data$z]), # all obs with lower x's, greater y's, AND lower z's (bottom upper left quadrant)
      gll = length(xth[xth > .data$x & yth < .data$y & zth < .data$z]), # all obs with greater x's, lower y's, AND lower z's (bottom lower right quadrant)
      ggl = length(xth[xth > .data$x & yth > .data$y & zth < .data$z]), # all obs with greater x's, greater y's, AND lower z's (bottom upper right quadrant)
      llg = length(xth[xth < .data$x & yth < .data$y & zth > .data$z]), # all obs with lower x's, lower y's, AND greater z's (top lower left quadrant)
      lgg = length(xth[xth < .data$x & yth > .data$y & zth > .data$z]), # all obs with lower x's, greater y's, AND greater z's (top upper left quadrant)
      glg = length(xth[xth > .data$x & yth < .data$y & zth > .data$z]), # all obs with greater x's, lower y's, AND greater z's (top lower right quadrant)
      ggg = length(xth[xth > .data$x & yth > .data$y & zth > .data$z]), # all obs with greater x's, greater y's, AND greater z's (top upper right quadrant)

      # the EXPECTED proportion of observations in a given quadrant:
      # We'd 'expect' it to be a product of the proportion of observations on the left/right side of x, and top/bottom side of Y
      expected1 = (xl/n) * (yl/n) * (zl/n), # bottom lower left
      expected2 = (xl/n) * (yg/n) * (zl/n), # bottom upper left
      expected3 = (xg/n) * (yl/n) * (zl/n), # bottom lower right
      expected4 = (xg/n) * (yg/n) * (zl/n), # bottom upper right
      expected5 = (xl/n) * (yl/n) * (zg/n), # top lower left
      expected6 = (xl/n) * (yg/n) * (zg/n), # top upper left
      expected7 = (xg/n) * (yl/n) * (zg/n), # top lower right
      expected8 = (xg/n) * (yg/n) * (zg/n), # top upper right

      # the OBSERVED proportion of observations in a given quadrant
      observed1 = lll/n, # bottom lower left
      observed2 = lgl/n, # bottom upper left
      observed3 = gll/n, # bottom lower right
      observed4 = ggl/n, # bottom upper right
      observed5 = llg/n, # top lower left
      observed6 = lgg/n, # top upper left
      observed7 = glg/n, # top lower right
      observed8 = ggg/n, # top upper right

      # compute differences in observed vs. expected proportions of observations within each quadrant
      d1 = abs(observed1 - expected1), # bottom lower left
      d2 = abs(observed2 - expected2), # bottom upper left
      d3 = abs(observed3 - expected3), # bottom lower right
      d4 = abs(observed4 - expected4), # bottom upper right
      d5 = abs(observed5 - expected5), # top lower left
      d6 = abs(observed6 - expected6), # top upper left
      d7 = abs(observed7 - expected7), # top lower right
      d8 = abs(observed8 - expected8), # top upper right

      # D-value: the maximum difference among all four quadrants
      maxdiff = max(d1, d2, d3, d4, d5, d6, d7, d8)

    ) # end dplyr data manipulation

  return(xyznew$maxdiff) # Function output

} # end quadrants function
