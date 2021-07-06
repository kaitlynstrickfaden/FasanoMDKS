#' Run a Two-dimensional Kolmogorov-Smirnov Test
#'
#' Run a two-dimensional Kolmogorov-Smirnov test as defined by [https://watermark.silverchair.com/mnras225-0155.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAs0wggLJBgkqhkiG9w0BBwagggK6MIICtgIBADCCAq8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMgHp1AmsYKnHXNEOdAgEQgIICgM8aZSoVKTeXKk_3u2Q16WV38bdHKZl47PPAZRnydI6VywiYEDK-o9G4LCYrsBtp9W7CQdG-sCouNnAg31KLOaX34sD0cYkMKIc0SELP0rLEkJIMTSsMBRdyn01EmO5hrH9_kGe24VPb9RIQq2LEB31sbzQmAOzwnGAOaPD7rlw0hn8cXxMgYjOpUDkJE9SDb5fhMmrZSOodnBcRAP3pdc_egvoYTnAmUqDSI-5AKHj9pErWPhVrAJ0HkdNTfI1CQLyXAqDSgoD8l109TgAfQAKk7WZ-Z0oP_wNT4ClVKV6cbMEIX1Hszd1JNShLs6qYg38RrZAi8lEAiONfvRnfBFdS7i_easW9jBrt0QWTcno-QopLWUXze6rE7qvyFjiqyB3g-UvKv0P0HOv4DgYw8STmvIZ_0KDh5YMn3y6uPG2nvP8ktudoaY1wZ2WZwE_IoxLRt6WZ0G1o9FOmnLw_tXRosMGKnCgV1YKSbshPC5rKZ8aoHOsmhhUfRxzwLv6ew9Am8M_2wINVaBn9HnGbu8WVRzRjvOlXTEaNk-ZUbTVgsPxCU2H4isCjrWbHXoCVGf4r8lM3a6Pc6oqXl3fIvZwvQmgPt3KgwFZxyqJHsXMJZWXPAZrU8c81EzdQ4816AyRGo8p0MlWmk7ArmBJ_uaOGu3jMGXnLJzJrzTZ3B-nPiLP7ctrfZsFcWipRuhGlWtL9PcKK3FBWa_gSprHUtEW_VHLihVHJOSU59qwoBao2sxsGQbI-3fKVSDofPsn-qgKKgPWtyM_K6CA_fDFjxkyUy7fy1hIm-ix50f9pMVsA8HpSPz4CDC7qywL5yrG8a3AeHzq8H8vesd33rfTuKDg](Fasano and Franceschini 1987). Allows for specification of the number of randomization to perform and the desired alpha level.
#'
#' @importFrom dplyr arrange group_by mutate summarize
#' @importFrom purrr map
#' @import progress
#' @param xy A data frame with an "x" and "y" column of paired observations.
#' @param rands A numeric indicating how many randomizations to perform. Default is 5000.
#' @param alpha A numeric between 0 and 1 for the desired alpha level. Default is 0.05.
#' @return A data frame of summary values and statistics (the maximum D-value, the x and y values at the maximum D-value, the minimum and maximum x and y values for the statistically-significant observations, and the p-value)
#' @export
#'

Fasano_2DKS <- function(xy, rands = 5000, alpha = 0.05) {

    if ("x" %in% colnames(xy) == FALSE | "y" %in% colnames(xy) == FALSE) {
        stop("xy must contain an 'x' and 'y' column")
    }

    if (rands < 2000) {
        warning("rands < 2000 may not converge", immediate. = T)
    }

    if (alpha < 0 | alpha > 1) {
        stop("alpha must be between 0 and 1")
    }


    ##### PART 1: D-values for each observation (i.e., potential threshold) #####
    ## Set observations as potential thresholds
    xy <- arrange(xy, x)    # sort rows by X in ascending order
    n <- nrow(xy)           # sample size
    xth <- xy$x             # isolate x column
    yth <- xy$y             # isolate y column


    ## Compute D-values for observed data
    d_obs <- Fasano_2D_dvals(xy)  # calculate D-values for observations


    ##### PART 2: p-values for each observation  #####

    # Progress bar for keeping track of analysis
    pb <- progress_bar$new(format = "Performing randomizations  Estimated completion: :eta  [:bar]  :percent",
                                  total = rands, clear = F)

    d_rand_sum <- 0 # empty list for storing randomizations of data

    for (r in 1:rands) { # for each randomization (total # set by input value 'rands'):

        xy_rand <- data.frame(x = rep(0,n), y = rep(0,n)) # empty data frame for re-shuffled X and Y values
        xy_rand$x <- xy$x[sample(c(1:n), n, replace = T)] # randomly index X-values (with replacement)
        xy_rand$y <- xy$y[sample(c(1:n), n, replace = T)] # randomly index Y-values (with replacement)

        xy_rand <- arrange(xy_rand, x)               # sort rows by X in ascending order and append to xy_rand list

        d_rand <- Fasano_2D_dvals(xy_rand)           # generate D-values for each 'observation'

        d_rand_max <- max(d_rand)                    # find the maximum D-value for this randomization (the observation that generated this would be the 'threshold')
        d_rand_bigger <- d_rand_max > d_obs    # (functionally) a count of how many observed D-values were >D max from this randomly generated dataset
        d_rand_sum <- d_rand_sum + d_rand_bigger # this list accumulates the above counts from each loop

        pb$tick()

    } # end "for" loop

    # compute p-values for each observation
    p <- d_rand_sum / rands # the proportion of randomly-generated Dmaxes (n = rands) that exceeded the observed Dmax



    ## Combine results of all potential thresholds
    xy_thresh <- data.frame(x = xth,       # true X-values
                            y = yth,       # true y-values
                            d_obs = d_obs, # true D-values
                            p = p          # true p-values
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
                      Min_sig_X = lo_x_val,
                      Max_sig_X = hi_x_val,
                      Min_sig_Y = lo_y_val,
                      Max_sig_Y = hi_y_val,
                      P_Value   = format(round(p_value, 4), nsmall = 4))
    )

} # end function




