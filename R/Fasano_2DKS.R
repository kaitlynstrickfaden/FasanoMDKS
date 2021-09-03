
#' Run a Two-dimensional Kolmogorov-Smirnov Test
#'
#' Run a two-dimensional Kolmogorov-Smirnov test as defined by \href{https://academic.oup.com/mnras/article/225/1/155/1007281}{Fasano and Franceschini 1987}. Allows for specification of the number of randomization to perform and the desired alpha level.
#'
#' @importFrom dplyr arrange group_by mutate summarize
#' @importFrom purrr map2_dbl pmap_dbl
#' @import progress
#' @param xcol a vector of x values
#' @param ycol a vector of y values
#' @param rands A numeric indicating how many randomizations to perform. Default is 5000.
#' @param alpha A numeric between 0 and 1 for the desired alpha level. Default is 0.05.
#' @return A data frame of summary values and statistics (the maximum D-value; the x, y, and p values at the maximum D-value; and the minimum and maximum x and y values for the statistically-significant observations)
#' @examples
#' xcol <- rnorm(10)
#' ycol <- rnorm(10)
#'
#' Fasano_2DKS(xcol, ycol, rands = 2000, alpha = 0.05)
#'
#' \dontrun{
#' xcol <- rnorm(10)
#' ycol <- rnorm(5)
#'
#' Fasano_3DKS(xcol, ycol, rands = 2000, alpha = 0.05)
#' }
#' @export
#'

Fasano_2DKS <- function(xcol, ycol, rands = 5000, alpha = 0.05) {

    if (length(xcol) != length(ycol)) {
        stop("xcol and ycol must be the same length")
    }

    if (rands < 2000) {
        warning("rands < 2000 may not converge", immediate. = T)
    }

    if (alpha < 0 | alpha > 1) {
        stop("alpha must be between 0 and 1")
    }


    ##### PART 1: D-values for each observation (i.e., potential threshold) #####
    ## Set observations as potential thresholds
    n <- length(xcol)                     # sample size

    ## Compute D-values for observed data
    d_obs <- Fasano_2D_dvals(xcol, ycol)  # calculate D-values for observations


    ##### PART 2: p-values for each observation  #####

    # Progress bar for keeping track of analysis
    pb <- progress_bar$new(format = "Performing randomizations  Estimated completion: :eta  [:bar]  :percent",
                                  total = rands, clear = F)

    d_rand_sum <- 0 # empty list for storing randomizations of data

    for (r in 1:rands) { # for each randomization (total # set by input value 'rands'):

        xrand <- xcol[sample(c(1:n), n, replace = T)] # randomly index X-values (with replacement)
        yrand <- ycol[sample(c(1:n), n, replace = T)] # randomly index Y-values (with replacement)

        d_rand <- Fasano_2D_dvals(xrand, yrand)           # generate D-values for each 'observation'

        d_rand_max <- max(d_rand)                    # find the maximum D-value for this randomization (the observation that generated this would be the 'threshold')
        d_rand_bigger <- d_rand_max > d_obs    # (functionally) a count of how many observed D-values were >D-max from this randomly generated dataset
        d_rand_sum <- d_rand_sum + d_rand_bigger # accumulates the above counts from each loop

        pb$tick()

    } # end "for" loop

    # compute p-values for each observation
    p <- d_rand_sum / rands # the proportion of randomly-generated Dmaxes (n = rands) that exceeded the observed Dmax



    ## Combine results of all potential thresholds
    xy_thresh <- data.frame(x = xcol,       # true X-values
                            y = ycol,       # true y-values
                            d_obs = d_obs,  # true D-values
                            p = p           # true p-values
    )

    ## Find "best" XY threshold and associated P-value in the observed dataset
    dmax    <- max(xy_thresh$d_obs)              # find Dmax
    drow    <- which(xy_thresh$d_obs == dmax)[1] # find the observation that generated Dmax (i.e., the threshold)
    best_x  <- xy_thresh$x[drow]                 # its x value
    best_y  <- xy_thresh$y[drow]                 # its y value
    p_value <- xy_thresh$p[drow]                 # its p value

    ## Find smallest and largest X threshold that is significant
    sig_xy   <- xy_thresh[xy_thresh$p <= alpha,] # select rows with significant p-values

    if (nrow(sig_xy) == 0) {

        lo_x_val <- NA                       # lowest 'significant' x value
        hi_x_val <- NA                       # highest 'significant' y value
        lo_y_val <- NA                       # lowest 'significant' x value
        hi_y_val <- NA                       # highest 'significant' y value

    } else {
        lo_x_val <- min(sig_xy$x)            # lowest 'significant' x value
        hi_x_val <- max(sig_xy$x)            # highest 'significant' y value
        lo_y_val <- min(sig_xy$y)            # lowest 'significant' x value
        hi_y_val <- max(sig_xy$y)            # highest 'significant' y value
    }



    ## Combine results
    return(data.frame(D_Max     = dmax,
                      Best_X    = best_x,
                      Best_Y    = best_y,
                      P_Value   = format(round(p_value, 4), nsmall = 4),
                      Min_Sig_X = lo_x_val,
                      Max_Sig_X = hi_x_val,
                      Min_Sig_Y = lo_y_val,
                      Max_Sig_Y = hi_y_val)
    )

} # end function




