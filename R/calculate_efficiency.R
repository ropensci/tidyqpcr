#' Calibrate primer sets / probes by calculating detection efficiency and
#' R squared
#'
#' Note efficiency is given in ratio, not per cent; multiply by 100 for that.
#'
#' @param cq_df_1 data frame with cq (quantification cycle) data,
#' 1 row per well.
#'
#' Must have columns cq, dilution.
#'
#' Assumes data are only for 1 probe/primer set/target_id, i.e. all values in
#' cq_df_1 are fit with the same slope.
#'
#' @param formula formula to use for log-log regression fit.
#'
#' Default value assumes multiple biological replicates,
#' cq ~ log2(dilution) + biol_rep.
#'
#' If only a single Biological Replicate, change to cq ~ log2(dilution).
#'
#' @return data frame with 1 single row, and columns:
#' efficiency, efficiency.sd, r.squared.
#'
#' @seealso calculate_efficiency_bytargetid
#'
#' @export
#' @importFrom tibble tibble
#' 
#' @examples
#' # create simple dilution dataset
#' dilution_tibble <- tibble(dilution = rep(c(1, 0.1, 0.001, 0.0001), 2),
#'                      cq = c(1, 3, 4, 6,
#'                             4, 5, 6, 7),
#'                      biol_rep = rep(c(1,2), each = 4),
#'                      target_id = "T1")
#'                      
#' # calculate primer efficiency
#' 
#' #----- use case 1: include difference across replicates in model
#' dilution_tibble %>%
#'     calculate_efficiency()
#' 
#' #----- use case 2: ignore difference across replicates
#' dilution_tibble %>%
#'     calculate_efficiency(formula = cq ~ log2(dilution))
#'
calculate_efficiency <- function(cq_df_1, formula = cq ~ log2(dilution) + biol_rep) {
    if (length(unique(cq_df_1$target_id)) > 1) {
            warning("multiple target_ids, did you mean calculate_efficiency_bytargetid?")
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
#' @param cq_df a data frame with cq (quantification cycle) data, 1 row per well
#'
#' Must have columns prep_type, target_id, cq, dilution.
#' Only prep_type=="+RT" columns are used.
#'
#' @param formula formula to use for log-log regression fit.
#'
#' Default value assumes multiple biological replicates,
#' cq ~ log2(dilution) + biol_rep.
#'
#' If only a single Biological Replicate, change to cq ~ log2(dilution).
#' If multiple sample_ids, change to cq ~ log2(dilution) + sample_id.
#'
#' See ?formula for background and help.
#'
#' @param use_prep_types prep_type column values to use, default "+RT" for RT-qPCR.
#'
#' By default, this includes only reverse-transcribed values in the efficiency
#' estimation, so excludes negative controls such as no-template and no-RT.
#'
#' To skip this filtering step, set use_prep_types=NA.
#'
#'
#' @return data frame with columns:
#' target_id, efficiency, efficiency.sd, r.squared.
#'
#' @seealso calculate_efficiency
#'
#' @export
#' @importFrom tidyr %>%
#' @importFrom rlang .data
#'
#' @examples
#' 
#' # create simple dilution dataset for two target_ids with two biological reps
#' dilution_tibble <- tibble(target_id = rep(c("T_1",
#'                                             "T_2"), each = 8),
#'                           well_row = rep(c("A",
#'                                            "B"), each = 8),
#'                           well_col = rep(1:8, 2),
#'                           well = paste0(well_row, well_col),
#'                           dilution = rep(c(1, 0.1, 0.001, 0.0001), 4),
#'                           cq = c(1, 3, 4, 6, 1, 3, 5, 7,
#'                                  4, 5, 6, 7, 3, 7, 8, 9),
#'                           biol_rep = rep(c(1, 1, 1, 1, 2, 2, 2, 2), 2),
#'                           prep_type = "+RT")
#'                      
#' # calculate primer efficiency for multiple targets
#' 
#' #----- use case 1: include difference across replicates in model
#' dilution_tibble %>%
#'     calculate_efficiency_bytargetid()
#' 
#' #----- use case 2: ignore difference across replicates
#' dilution_tibble %>%
#'     calculate_efficiency_bytargetid(formula = cq ~ log2(dilution))
#'
calculate_efficiency_bytargetid <- function(cq_df,
                           formula = cq ~ log2(dilution) + biol_rep,
                           use_prep_types="+RT") {
    if (!is.na(use_prep_types)) {
        cq_df <- dplyr::filter(cq_df, .data$prep_type %in% use_prep_types)
    }
    cq_df %>%
        dplyr::group_by(.data$target_id) %>%
        dplyr::do(calculate_efficiency(.data, formula = formula)) %>%
        dplyr::ungroup()
}