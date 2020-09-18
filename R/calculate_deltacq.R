
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
calculate_normvalue <- function(value_df,
                      ref_ids,
                      value_name = "value",
                      id_name = "id",
                      norm_function = median) {
    # make subset of value_df where gene is one or more ref_ids
    value_to_norm_by <- dplyr::filter(value_df,
                             !!dplyr::sym(id_name) %in% ref_ids) %>%
        dplyr::pull(!!dplyr::sym(value_name)) %>%
        norm_function(na.rm = TRUE)
    #
    # assign summary (median) value to value_df$value_to_norm_by
    # note this is the same value for every row, a waste of space technically
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
calculate_deltacq_bysampleid <- function(cq_df,
                                         ref_target_ids,
                                         norm_function = median) {
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
calculate_deltadeltacq_bytargetid <- function(deltacq_df,
                                         ref_sample_ids,
                                         norm_function = median) {
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
               deltadelta_cq = .data$delta_cq - .data$ref_delta_cq,
               fold_change   = 2^-.data$deltadelta_cq)
}