

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


#' Calculate dR/dT for a melt curve of signal R vs temperature T
#'
#' @param TT Temperature.
#' @param RR Signal, assumed fluorescence signal matched to TT.
#' @param method to use for smoothing:
#' 
#'   "spline" default, uses smoothing spline stats::smooth.spline.
#'   
#'   "diff" base::diff for lagged difference
#'
#' @param ... other arguments to pass to smoothing method.
#'
#' @return estimated first derivative of RR with respect to TT,
#' numeric vector of same length as RR.
#'
#' @family melt_curve_functions
#'
#' @export
#'
getdRdT <- function(TT, RR, method = "spline", ...) {
    if (method == "diff") {
       return(-c(diff(RR) / diff(TT), NA))
    } else if (method == "spline") {
        fit <- stats::smooth.spline(x = TT, y = RR, ...)
        return(-1 * stats::predict(object = fit, x = TT, deriv = 1)$y)
    }
}


#' Calculate dR/dT for melt curves of every well in a plate
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
getdRdTall <- function(platemelt, method = "spline", ...) {
    platemelt %>%
        dplyr::arrange(Well, Temperature) %>%
        # @ewallace: doesn't group by plate, only by well,
        # so will fail strangely if used on data from multiple plates
        dplyr::group_by(Well) %>%
        dplyr::mutate(dRdT = getdRdT(Temperature, 
                                     Fluor,
                                     method=method, 
                                     ...)
                      ) %>%
        dplyr::ungroup()
}
