# Dependencies for the original cross-validation and energy-table functions.
#' @importFrom dplyr arrange bind_rows filter
#' @importFrom magrittr %>%
NULL

# Silence false-positive global variable checks triggered by dplyr data masking.
utils::globalVariables(c("P", "Q", "R", "CV_MSE", "CV_SE",
                         "ModelDF", "RelativeEnergy"))
