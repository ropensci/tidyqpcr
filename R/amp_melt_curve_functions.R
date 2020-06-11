

#' Remove baseline from amplification curves (BETA)
#'
#' Remove baseline from qPCR amplification curves
#' by subtracting median of initial cycles.
#'
#' BETA function version because:
#'
#' - assumes Roche Lightcycler format,
#' we should ideally replace "program==2" by something more generic?
#'
#' - the rule-of thumb "baseline is median of initial 10 cycles"
#' has not been tested robustly
#'
#' @param plateamp data frame with plate amplification data, including variables
#'   Well, Cycle, Fluor (fluorescence value), and Program. Assume program 2 for
#'   cmplification curves from Roche Lightcycler format data.
#' @param maxcycle maximum cycle value to use for baseline, before
#'   amplification.
#'
#' @return platemap with additional columns per Well:
#' \tabular{ll}{
#'   Base      \tab baseline /background value  \cr
#'   Signal    \tab normalized fluorescence signal, i.e. Fluor - Base
#'   }
#'
#'
#' @export
#' @importFrom magrittr %>%
#'
debaseline <- function(plateamp, maxcycle = 10) {
    baseline <-
        plateamp %>%
        dplyr::group_by(Well) %>%
        dplyr::filter(Program == 2, Cycle <= maxcycle) %>%
        dplyr::summarize(Base = stats::median(Fluor))
    plateamp %>%
        dplyr::left_join(baseline) %>%
        dplyr::mutate(Signal = Fluor - Base)
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
#' this suggests a single-liength PCR product is present in the well.
#'
#' @param platemelt data frame describing melt curves, including variables
#'   Well, Temperature, Fluor (fluorescence value).
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
#' @importFrom magrittr %>%
#'
calculate_drdt_plate <- function(platemelt, method = "spline", ...) {
    platemelt %>%
        dplyr::arrange(Well, Temperature) %>%
        # @ewallace: doesn't group by plate, only by well,
        # so will fail strangely if used on data from multiple plates
        dplyr::group_by(Well) %>%
        dplyr::mutate(dRdT =
                          calculate_dydx_1(x = Temperature,
                                           y = Fluor,
                                           method = method,
                                           ...)
                      ) %>%
        dplyr::ungroup()
}

#' @describeIn calculate_drdt_plate
#'
#' @export
#'
getdRdTall <- function(platemelt, method = "spline") {
    lifecycle::deprecate_warn("0.2", "getdRdTall()",
                              "calculate_drdt_plate()",
        details = "Replaced with more specific name")
    calculate_drdt_plate(platemelt = platemelt,
                         method = method)
}
