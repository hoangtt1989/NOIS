#' Quantile Thresholding
#'
#' Perform quantile thresholding
#' @param input A numeric or integer vector.
#' @param thresh_val An integer specifying the number of values to keep (all others set to zero).
#' @return The thresholded input.
#' @export
quantile_thresh <- function(input, thresh_val) {
    sortobj <- sort(abs(input), decreasing = TRUE, index.return = TRUE)
    outvec <- rep(0, length(input))
    thresh_ind <- sortobj$ix[1:thresh_val]
    outvec[thresh_ind] <- input[thresh_ind]
    return(outvec)
}
