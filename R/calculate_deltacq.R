#' Calculate a normalized value for a subset of reference ids
#'
#' This is used to calculate the normalized `cq` values for reference
#' `target_ids` (e.g. genes), to use in `delta_cq` calculation for each
#' `sample_id`.
#'
#' Also used to calculate the normalized `delta_cq` values for reference
#' `sample_ids`, to use in `deltadelta_cq` calculation for each `target_id`.
#'
#' @param value_df data frame containing relevant columns, those named in
#'   `value_name` and `id_name` parameters.
#' @param ref_ids values of reference ids, that are used to calculate
#'   normalized reference value.
#' @param value_name name of column containing values. This column should be
#'   numeric.
#' @param id_name name of column containing ids.
#' @param norm_function Function to use to calculate the value to normalize by.
#'   Default function is median, alternatively could use mean, geometric mean,
#'   etc.
#'
#' @export
#' @importFrom tidyr %>%
#' @importFrom stats median
#'
#' @examples
#' # create simple cq dataset with one sample, two targets  and 3 reps
#' cq_tibble <- tibble(sample_id = "S_1",
#'                      target_id = rep(c("T_1",
#'                                        "T_norm"), each = 3),
#'                      tech_rep = rep(1:3, 2),
#'                      well_row = rep(c("A",
#'                                       "B"), each = 3),
#'                      well_col = 1,
#'                      well = paste0(well_row, well_col),
#'                      cq = c(10, 10, 10,
#'                             12, 12, 11))
#'                      
#' # normalise cq to reference target_id called 'T_norm'
#' 
#' #----- use case 1: median reference target_id value
#' cq_tibble %>%
#'     calculate_normvalue(ref_ids = "T_norm",
#'                         value_name = "cq",
#'                         id_name = "target_id")
#' 
#' #----- use case 2: mean reference target_id value 
#' cq_tibble %>%
#'     calculate_normvalue(ref_ids = "T_norm",
#'                         value_name = "cq",
#'                         id_name = "target_id",
#'                         norm_function = mean)
#'
calculate_normvalue <- function(value_df,
                      ref_ids,
                      value_name = "value",
                      id_name = "id",
                      norm_function = median) {
    
    assertthat::assert_that(
        assertthat::has_name(value_df, 
                             c(value_name,id_name)))
    
    # make subset of value_df where gene is one or more ref_ids
    value_to_norm_by <- dplyr::filter(value_df,
                             .data[[id_name]] %in% ref_ids) %>%
        dplyr::pull({{value_name}}) %>%
        norm_function(na.rm = TRUE)
    
    # assign summary (median) value to value_df$value_to_norm_by
    dplyr::mutate(value_df, value_to_norm_by = value_to_norm_by)
}

#' Calculate delta cq (\eqn{\Delta Cq}) to normalize quantification cycle
#' (log2-fold) data within sample_id.
#'
#' This function implements relative quantification by the delta Cq method. For
#' each sample, the Cq values of all targets (e.g. genes, probes, primer sets)
#' are compared to one or more reference target ids specified in
#' `ref_target_ids`.
#'
#' @param cq_df a data frame containing columns `sample_id`, value_name (default
#'   `cq`) and tid_name (default `target_id`). Crucially, sample_id should be
#'   the same for different technical replicates measuring identical reactions
#'   in different wells of the plate, but differ for different biological and
#'   experimental replicates. See tidyqpcr vignettes for examples.
#' @param ref_target_ids names of targetss to normalize by, i.e. reference
#'   genes, hydrolysis probes, or primer sets. This can be one reference target
#'   id, a selection of multiple target ids, or even all measured target ids. In
#'   the case of all of them, the delta Cq value would be calculated relative to
#'   the median (or other `norm_function`) of all measured targets.
#' @param norm_function Function to use to calculate the value to normalize by
#'   on given scale. Default is median, alternatively could use mean.
#'
#' @return data frame like cq_df with three additional columns:
#'
#'   \tabular{ll}{ ref_cq    \tab summary (median/mean) cq value for reference
#'   target ids \cr delta_cq  \tab normalized value, \eqn{\Delta Cq} \cr
#'   rel_abund \tab normalized ratio, \eqn{2^(-\Delta Cq)} }
#'
#' @export
#' @importFrom tidyr %>%
#' @importFrom stats median
#' @importFrom rlang .data
#' 
#' @examples
#' # create simple cq dataset with two samples, two targets  and 3 reps
#' cq_tibble <- tibble(sample_id = rep(c("S_1","S_1","S_1", "S_2", "S_2", "S_2"), 2),
#'                      target_id = rep(c("T_1",
#'                                        "T_norm"), each = 6),
#'                      tech_rep = rep(1:3, 4),
#'                      well_row = rep(c("A",
#'                                       "B"), each = 6),
#'                      well_col = rep(1:6, 2),
#'                      well = paste0(well_row, well_col),
#'                      cq = c(10, 10, 10, 12, 12, 11,
#'                              9,  9,  9,  9,  9,  9))
#'                      
#' # calculate deltacq using reference target_id called 'T_norm'
#' 
#' #----- use case 1: median reference target_id value
#' cq_tibble %>%
#'     calculate_deltacq_bysampleid(ref_target_ids = "T_norm")
#' 
#' #----- use case 2: mean reference target_id value 
#' cq_tibble %>%
#'     calculate_deltacq_bysampleid(ref_target_ids = "T_norm",
#'                                  norm_function = mean)
#'
calculate_deltacq_bysampleid <- function(cq_df,
                                         ref_target_ids,
                                         norm_function = median) {
    
    assertthat::assert_that(
        assertthat::has_name(cq_df, 
                             c("target_id", "sample_id","cq")))
    
    cq_df %>%
        dplyr::group_by(.data$sample_id) %>%
        dplyr::do(calculate_normvalue(.data,
                                   ref_ids = ref_target_ids,
                                   value_name = "cq",
                                   id_name = "target_id",
                                   norm_function = norm_function)) %>%
        dplyr::rename(ref_cq = .data$value_to_norm_by) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
               delta_cq    = .data$cq - .data$ref_cq,
               rel_abund   = 2^-.data$delta_cq)
}


#' Calculate delta delta cq (\eqn{\Delta \Delta Cq}) to globally normalize
#' quantification cycle (log2-fold) data across sample_id.
#'
#' By default, \eqn{\Delta \Delta Cq} is positive if a target is more highly
#' detected in the relevant sample, compared to reference samples. This can be
#' flipped by setting the parameter `ddcq_positive` to `FALSE`. In either case,
#' The fold change, \eqn{2^{\Delta \Delta Cq}}, is also reported.
#'
#' This function does a global normalization, where all samples are compared to
#' one or more reference samples specified in `ref_sample_ids`. There are other
#' experimental designs that require comparing samples in pairs or small groups,
#' e.g. a time course comparing `delta_cq` values against a reference strain at
#' each time point. For those situations, instead we recommend adapting code
#' from this function, changing the grouping variables used in to
#' `dplyr::group_by` to draw the contrasts appropriate for the experiment.
#'
#' @param deltacq_df a data frame containing columns `sample_id`, value_name
#'   (default `delta_cq`) and tid_name (default `target_id`). Crucially,
#'   sample_id should be the same for different technical replicates measuring
#'   identical reactions in different wells of the plate, but differ for
#'   different biological and experimental replicates.
#'
#'   Usually this will be a data frame that was output by
#'   `calculate_deltacq_bysampleid`.
#'
#' @param ref_sample_ids reference sample_ids to normalize by
#' @param norm_function Function to use to calculate the value to normalize by
#'   on given scale. Default is median, alternatively could use mean.
#' @param ddcq_positive (default TRUE) output \eqn{\Delta \Delta Cq} as positive
#'   if a target is more highly detected in the relevant sample, compared to
#'   reference samples.
#'
#' @return data frame like cq_df with three additional columns:
#'
#'   \tabular{ll}{ ref_delta_cq  \tab summary (median/mean) \eqn{\Delta Cq}
#'   value for target_id in reference sample ids \cr deltadelta_cq \tab the
#'   normalized value, \eqn{\Delta \Delta Cq} \cr fold_change   \tab the
#'   normalized fold-change ratio, \eqn{2^(-\Delta \Delta Cq)} }
#'
#' @export
#' @importFrom tidyr %>%
#' @importFrom stats median
#'
#' @examples
#' # create simple deltacq dataset with two samples, two targets and 3 reps
#' deltacq_tibble <- tibble(sample_id = rep(c("S_1","S_1","S_1", "S_norm", "S_norm", "S_norm"), 2),
#'                      target_id = rep(c("T_1",
#'                                        "T_2"), each = 6),
#'                      tech_rep = rep(1:3, 4),
#'                      well_row = rep(c("A",
#'                                       "B"), each = 6),
#'                      well_col = rep(1:6, 2),
#'                      well = paste0(well_row,well_col),
#'                      delta_cq = c(1, 1, 1, 3, 3, 2,
#'                                   4, 5, 4, 5, 5, 5))
#'                      
#' # calculate deltadeltacq using reference target_id called 'S_norm'
#' 
#' #----- use case 1: median reference sample_id value
#' deltacq_tibble %>%
#'     calculate_deltadeltacq_bytargetid(ref_sample_ids = "S_norm")
#' 
#' #----- use case 2: mean reference sample_id value 
#' deltacq_tibble %>%
#'     calculate_deltadeltacq_bytargetid(ref_sample_ids = "S_norm",
#'                                  norm_function = mean)
#'
calculate_deltadeltacq_bytargetid <- function(deltacq_df,
                                         ref_sample_ids,
                                         norm_function = median,
                                         ddcq_positive = TRUE) {
    
    assertthat::assert_that(
        assertthat::has_name(deltacq_df, 
                             c("target_id", "sample_id","delta_cq")))
    
    ddcq_factor <- (-1) ^ ddcq_positive
    
    deltacq_df %>%
        dplyr::group_by(.data$target_id) %>%
        dplyr::do(calculate_normvalue(.data,
                                   ref_ids = ref_sample_ids,
                                   value_name = "delta_cq",
                                   id_name = "sample_id",
                                   norm_function = norm_function)) %>%
        dplyr::rename(ref_delta_cq = .data$value_to_norm_by) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
               deltadelta_cq = ddcq_factor *
                   (.data$delta_cq - .data$ref_delta_cq),
               fold_change   = 2 ^ .data$deltadelta_cq)
}