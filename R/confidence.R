#' confidence
#'
#' This function annotated a confidence level to the score
#' @param s the score
#' @param type coding or noncoding
#'
#' @return returns a confidence level or remark/error message
#' @export
#' @examples
#'
#' remark <- confidence(0.8, type='Coding')


confidence <- function(s,
                       type) {
    # Returns the confidence level for (neutral, low-confidence or
    # high-confidence) for the given p-value string using the given threshold.
    # If the string is invalid, returns "no prediction".

    # Initialization:
    # Confidence thresholds
    DEFAULT_THRESH <- 0.5
    CODING_THRESH <- 0.89
    NONCODING_THRESH <- 0.7
    if (type == "Coding") {
        thresh <- CODING_THRESH
    } else {
    thresh <- NONCODING_THRESH
    }

    # Function:
    x <- as.numeric(s)
    if (x >= thresh) {
        return("High-confidence")
    } else if (x >= DEFAULT_THRESH) {
        return("Low-confidence")
    } else {
        return("Neutral")
    }
}
