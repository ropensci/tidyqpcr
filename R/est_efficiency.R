
#' Calibrate primer sets / probes by calculating detection efficiency and
#' R squared
#'
#' Note efficiency is given in ratio, not per cent; multiply by 100 for that.
#'
#' @param cq_df_1 data frame with Cq (quantification cycle) data,
#' 1 row per well.
#'
#' Must have columns Cq, Dilution.
#'
#' Assumes data are only for 1 probe/primer set/target_id, i.e. all values in
#' cq_df_1 are fit with the same slope.
#'
#' @param formula formula to use for log-log regression fit.
#'
#' Default value assumes multiple biological replicates,
#' Cq ~ log2(Dilution) + BioRep.
#'
#' If only a single Biological Replicate, change to Cq ~ log2(Dilution).
#'
#' @return data frame with 1 single row, and columns:
#' efficiency, efficiency.sd, r.squared.
#'
#' @seealso est_efficiency_bytargetid
#'
#' @export
#' @importFrom tibble tibble
#'
est_efficiency <- function(cq_df_1, formula = Cq ~ log2(Dilution) + BioRep) {
    if (length(unique(cq_df_1$target_id)) > 1 ) {
            warning("multiple target_ids, did you mean to run est_efficiency_plate?")
    }
    slopefit <- stats::lm(formula = formula, data = cq_df_1)
    slopefitsummary <- summary(slopefit)
    tibble(efficiency = -slopefit$coefficients[2],
           efficiency.sd = slopefitsummary$coefficients[2, 2],
           r.squared = slopefitsummary$r.squared)
}

#' Calibrate multiple probes by calculating detection efficiency and R squared
#'
#' See calibration vignette for example of usage.
#'
#' Note efficiency is given in ratio, not per cent; multiply by 100 for that.
#'
#' @param cq_df a data frame with Cq (quantification cycle) data, 1 row per well
#'
#' Must have columns Type, target_id, Cq, Dilution.
#' Only Type=="+RT" columns are used.
#'
#' @param formula formula to use for log-log regression fit.
#'
#' Default value assumes multiple biological replicates,
#' Cq ~ log2(Dilution) + BioRep.
#'
#' If only a single Biological Replicate, change to Cq ~ log2(Dilution).
#' If multiple sample_ids, change to Cq ~ log2(Dilution) + sample_id.
#'
#' See ?formula for background and help.
#'
#' @param usetypes Type column values to use, default "+RT" for RT-qPCR.
#'
#' By default, this includes only reverse-transcribed values in the efficiency
#' estimation, so excludes negative controls such as no-template and no-RT.
#'
#' To skip this filtering step, set usetypes=NA.
#'
#'
#' @return data frame with columns:
#' target_id, efficiency, efficiency.sd, r.squared.
#'
#' @seealso est_efficiency
#'
#' @export
#' @importFrom magrittr %>%
#'
est_efficiency_bytargetid <- function(cq_df,
                           formula = Cq ~ log2(Dilution) + BioRep,
                           usetypes="+RT") {
    if (!is.na(usetypes)) {
        cq_df <- dplyr::filter(cq_df, Type %in% usetypes)
    }
    cq_df %>%
        dplyr::group_by(target_id) %>%
        dplyr::do(est_efficiency(., formula = formula))
}
