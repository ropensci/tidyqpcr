
#' Calibrate probes by calculating detection efficiency and R squared
#' 
#' Note efficiency is given in ratio, not per cent; multiply by 100 for that.
#'
#' @param ct_df_1 a data frame with Ct data, 1 row per well.
#' 
#' Must have columns Ct, Dilution. 
#' 
#' Assumes data are only for 1 probe, i.e. all values in ct_df_1 are fit with the same slope.
#' 
#' @param formula formula to use for log-log regression fit. 
#' 
#' Default value assumes multiple biological replicates Ct ~ log2(Dilution) + BioRep.
#' If only a single Biological Replicate, change to Ct ~ log2(Dilution).
#'
#' @return data frame with 1 single row, and columns: efficiency, efficiency.sd, r.squared
#'
#' @seealso est_efficiency
#' 
#' @export
#' @importFrom tibble tibble
#' 
est_efficiency_1 <- function(ct_df_1,formula = Ct ~ log2(Dilution) + BioRep) {
    slopefit <-  stats::lm(formula = formula,
                    data    = ct_df_1)
    slopefitsummary <- summary(slopefit)
    tibble(efficiency = - slopefit$coefficients[2],
           efficiency.sd = slopefitsummary$coefficients[2,2],
           r.squared = slopefitsummary$r.squared)
}

#' Calibrate multiple probes by calculating detection efficiency and R squared
#' 
#' See calibration vignette for example of usage.
#' 
#' Note efficiency is given in ratio, not per cent; multiply by 100 for that.
#'
#' @param ct_df a data frame with Ct data, 1 row per well
#' 
#' Must have columns Type, Probe, Ct, Dilution. Only Type=="+RT" columns are used.
#' 
#' @param formula formula to use for log-log regression fit. 
#' 
#' Default value assumes multiple biological replicates Ct ~ log2(Dilution) + BioRep.
#' If only a single Biological Replicate, change to Ct ~ log2(Dilution).
#' If multiple SampleIDs, change to Ct ~ log2(Dilution) + SampleID. See ?formula for background
#'
#' @return data frame with columns: Probe, efficiency, efficiency.sd, r.squared
#' 
#' @seealso est_efficiency_1
#' 
#' @export
#' @importFrom magrittr %>%
#' 
est_efficiency <- function(ct_df,formula = Ct ~ log2(Dilution) + BioRep) {
    ct_df %>%
        dplyr::filter(Type=="+RT") %>%
        dplyr::group_by(Probe) %>%
        dplyr::do(est_efficiency_1(.,formula=formula))
}

