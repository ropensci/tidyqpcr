#' Remove baseline from amplification curves (BETA)
#'
#' Remove baseline from qPCR amplification curves
#' by subtracting median of initial cycles.
#'
#' BETA function version because:
#'
#' - assumes Roche Lightcycler format,
#' we should ideally replace "program_no == 2" by something more generic?
#'
#' - the rule-of thumb "baseline is median of initial 10 cycles"
#' has not been tested robustly
#'
#' @param plateamp data frame with plate amplification data, including variables
#'   well, cycle, fluor_raw (raw fluorescence value), and program_no. Assume
#'   program 2 for amplification curves from Roche Lightcycler format data.
#' @param maxcycle maximum cycle value to use for baseline, before
#'   amplification.
#'
#' @return platemap with additional columns per well:
#' \tabular{ll}{
#'   fluor_base      \tab baseline /background value  \cr
#'   fluor_signal    \tab normalized fluorescence signal, 
#'   i.e. fluor_raw - fluor_base
#'   }
#'
#' @export
#' @importFrom tidyr %>%
#' @importFrom rlang .data
#'
#' @examples
#' # create simple dataset of raw fluorescence
#' # with two samples over 15 cycles
#' raw_fluor_tibble <- tibble(sample_id = rep(c("S1", "S2"), each = 15),
#'                           target_id = "T1",
#'                           well_row = "A",
#'                           well_col = rep(c(1, 2), each = 15),
#'                           well = rep(c("A1", "A2"), each = 15),
#'                           cycle = rep(1:15,2),
#'                           fluor_raw = c(1:15, 6:20),
#'                           program_no = 2)
#'
#' # remove base fluorescence from dataset
#' raw_fluor_tibble %>%
#'     debaseline()
#'     
debaseline <- function(plateamp, maxcycle = 10) {
    assertthat::assert_that(
        assertthat::has_name(
            plateamp,
            c("well", "program_no", "cycle", "fluor_raw"))
    )
    baseline <-
        plateamp %>%
        dplyr::group_by(.data$well) %>%
        dplyr::filter(.data$program_no == 2,
                      .data$cycle <= maxcycle) %>%
        dplyr::summarize(fluor_base = stats::median(.data$fluor_raw))
    plateamp %>%
        dplyr::left_join(baseline) %>%
        dplyr::mutate(fluor_signal = .data$fluor_raw - .data$fluor_base)
}


#' Calculate dy/dx vector from vectors y and x
#'
#' Used in tidyqpcr to calculate dR/dT for a melt curve of fluorescence signal R
#' vs temperature T.
#'
#' @param x input variable, numeric vector, assumed to be temperature
#' @param y output variable, numeric vector of same length as x, assumed
#'   to be fluorescence signal.
#' @param method to use for smoothing:
#'
#'   "spline" default, uses smoothing spline stats::smooth.spline.
#'
#'   "diff" base::diff for lagged difference
#'
#' @param ... other arguments to pass to smoothing method.
#'
#' @return estimated first derivative of y with respect to x, numeric vector of
#'   same length as y.
#'
#' @family melt_curve_functions
#'
#' @export
#' 
#' @examples
#' # create simple curve
#' x = 1:5
#' y = x^2
#'
#' # calculate gradient of curve
#' #----- use case 1 : using splines
#' calculate_dydx_1(x, y)
#' 
#' # optional arguments are passed to smooth.splines function
#' calculate_dydx_1(x, y, spar = 0.5)
#' 
#' #----- use case 2 : using difference between adjacent points
#' calculate_dydx_1(x, y, method = "diff")
#'
calculate_dydx_1 <- function(x, y, method = "spline", ...) {
    assertthat::assert_that(is.numeric(x))
    assertthat::assert_that(is.numeric(y))
    assertthat::assert_that(length(x) == length(y))
    if (method == "diff") {
       return(-c(diff(y) / diff(x), NA))
    } else if (method == "spline") {
        fit <- stats::smooth.spline(x = x, y = y, ...)
        return(-1 * stats::predict(object = fit, x = x, deriv = 1)$y)
    }
}


#' Calculate dR/dT of melt curves for of every well in a plate.
#'
#' dR/dT, the derivative of the melt curve (of fluorescence signal R vs
#' temperature T), has a maximum at the melting temperature Tm. A single peak in
#' this suggests a single-length PCR product is present in the well.
#' 
#' Note that this function does not group by plate, only by well.
#' The function will give strange results if you pass it data from 
#' more than one plate. Avoid this by analysing one plate at a time.
#'
#' @param platemelt data frame describing melt curves, including variables
#'   well, temperature, fluor_raw (raw fluorescence value).
#' @param method to use for smoothing:
#'
#'   "spline" default, uses smoothing spline stats::smooth.spline.
#'
#'   "diff" base::diff for lagged difference
#'
#' @param ... other arguments to pass to smoothing method.
#'
#' @return platemelt with additional column dRdT.
#'
#' @family melt_curve_functions
#'
#' @export
#' @importFrom tidyr %>%
#'
#' @examples
#' # create simple curve
#' # create simple dataset of raw fluorescence with two samples
#' temp_tibble <- tibble(sample_id = rep(c("S1", "S2"), each = 10),
#'                           target_id = "T1",
#'                           well_row = "A",
#'                           well_col = rep(c(1, 2), each = 10),
#'                           well = rep(c("A1", "A2"), each = 10),
#'                           temperature = rep(56:65,2),
#'                           fluor_raw = c(1:10, 6:15))
#'
#' # calculate drdt of all melt curves
#' #----- use case 1 : using splines
#' temp_tibble %>%
#'     calculate_drdt_plate()
#' 
#' # optional arguments are passed to smooth.splines function
#' temp_tibble %>%
#'     calculate_drdt_plate(spar = 0.5)
#' 
#' #----- use case 2 : using difference between adjacent points
#' temp_tibble %>%
#'     calculate_drdt_plate(method = "diff")
#'
calculate_drdt_plate <- function(platemelt, method = "spline", ...) {
    assertthat::assert_that(
        assertthat::has_name(
            platemelt,
            c("well", "temperature", "fluor_raw"))
    )
    if (assertthat::has_name(platemelt,"plate") ) {
        warning("platemelt has a plate column, but calculate_drdt_plate works only for single plates.")
    }
    platemelt %>%
        dplyr::arrange(.data$well, .data$temperature) %>%
        dplyr::group_by(.data$well) %>%
        dplyr::mutate(dRdT =
                          calculate_dydx_1(x = .data$temperature,
                                           y = .data$fluor_raw,
                                           method = method,
                                           ...)
                      ) %>%
        dplyr::ungroup()
}
