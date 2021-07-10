
#' Calculate D-values for a Two-dimensional Kolmogorov-Smirnov Test
#'
#' Find D-values of observations for a two-dimensional Kolmogorov-Smirnov test as defined by [Fasano and Franceschini 1987](https://academic.oup.com/mnras/article/225/1/155/1007281).
#'
#' @importFrom dplyr arrange group_by mutate summarize
#' @param xy The D-values for each observation as a vector.
#' @return The original data frame with a new "maxdiff" column containing D-values.
#' @export



## A function to determine D-values
Fasano_2D_dvals <- function(xy) {

  if ("x" %in% colnames(xy) == FALSE | "y" %in% colnames(xy) == FALSE) {
    stop("xy must contain an 'x' and 'y' column")
  }

  n <- nrow(xy)  # sample size (repeated here to allow quadrants to run as its own function)
  xth <- xy$x    # isolate x column
  yth <- xy$y    # isolate y column

  # manipulate input data to find D values
  xynew <- xy %>%

    arrange(x) %>%

    mutate(i = 1:n) %>% # add row index

    group_by(i) %>% # group by row index (need to summarize each row individually rather than the whole data frame)

    summarize(

      # store of values less/greater than a given x/y value
      xl = length(xth[xth < .data$x]), # all x-values less than the X
      xg = length(xth[xth > .data$x]), # all x-values greater than the X
      yl = length(yth[yth < .data$y]), # all y-values less than the Y
      yg = length(yth[yth > .data$y]), # all x-values greater than the Y

      # store of values for a given quadrant
      ll = length(xth[xth < .data$x & yth < .data$y]), # all obs with lower x's AND lower y's (lower left quadrant)
      lg = length(xth[xth < .data$x & yth > .data$y]), # all obs with lower x's AND greater y's (upper left quadrant)
      gl = length(xth[xth > .data$x & yth < .data$y]), # all obs with greater x's AND lower y's (lower right quadrant)
      gg = length(xth[xth > .data$x & yth > .data$y]), # all obs with greater x's AND greater y's (upper right quadrant)

      # the EXPECTED proportion of observations in a given quadrant:
      # We'd 'expect' it to be a product of the proportion of observations on the left/right side of x, and top/bottom side of Y
      expected1 = (xl/n) * (yl/n), # lower left
      expected2 = (xl/n) * (yg/n), # upper left
      expected3 = (xg/n) * (yl/n), # lower right
      expected4 = (xg/n) * (yg/n), # upper right

      # the OBSERVED proportion of observations in a given quadrant
      observed1 = ll/n, # lower left
      observed2 = lg/n, # upper left
      observed3 = gl/n, # lower right
      observed4 = gg/n, # upper right

      # compute differences in observed vs. expected proportions of observations within each quadrant
      d1 = abs(observed1 - expected1), # lower left
      d2 = abs(observed2 - expected2), # upper left
      d3 = abs(observed3 - expected3), # lower right
      d4 = abs(observed4 - expected4), # upper right

      # D-value: the maximum difference among all four quadrants
      maxdiff = max(d1, d2, d3, d4)

    ) # end dplyr data manipulation

  return(xynew$maxdiff) # Function output

} # end quadrants function
