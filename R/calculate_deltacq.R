
#' @describeIn calculate_deltacq_bysampleid get the median value of a set of
#'   normalization (reference) probes, for a single sample.
#'
#' @param norm_function Function to use to calculate the value to
#' normalise by on log2/cq scale.
#' Default function is median, alternatively could use mean.
#'
#' @export
#' @importFrom tidyr %>%
#' @importFrom stats median
#'
calculate_normcq <- function(cq_df,
                      value_name = "cq",
                      norm_target_ids = "ALG9",
                      tid_name = "target_id",
                      norm_function = median) {
    # make subset of cq_df where gene is one or more norm_target_ids
    value_to_norm_by <- dplyr::filter(cq_df,
                             !!dplyr::sym(tid_name) %in% norm_target_ids) %>%
        dplyr::pull(!!dplyr::sym(value_name)) %>%
        norm_function(na.rm = TRUE)
    #
    # assign summary (median) value to cq_df$value_to_norm_by
    # note this is the same value for every row, a waste of space technically
    dplyr::mutate(cq_df, value_to_norm_by = value_to_norm_by)
}

#' Calculate delta cq to normalize quantification cycle (log2-fold) data within
#' sample_id.
#'
#' @param cq_df a data frame containing columns `sample_id`, value_name (default
#'   `cq`) and tid_name (default `target_id`). Crucially, sample_id should be
#'   the same for different technical replicates measuring identical reactions
#'   in different wells of the plate, but differ for different biological and
#'   experimental replicates.
#' @param value_name the column name of the value that will be normalized
#' @param norm_target_ids names of PCR probes (or primer sets) to normalize by,
#'   i.e. reference genes
#' @param tid_name the column name for probe sets
#'
#' @return data frame like cq_df with three additional columns:
#'
#'   \tabular{ll}{ value_to_norm_by       \tab the median value of the reference
#'   probes \cr value_norm    \tab the normalized value, \eqn{\Delta Cq} \cr
#'   value_normexp \tab the normalized ratio, \eqn{2^(-\Delta Cq)} }
#'
#' @export
#' @importFrom tidyr %>%
#'
calculate_deltacq_bysampleid <- function(cq_df,
                                         norm_target_ids,
                                         value_name = "cq",
                                         tid_name = "target_id") {
    cq_df %>%
        dplyr::group_by(sample_id) %>%
        dplyr::do(calculate_normcq(.,
                                   value_name,
                                   norm_target_ids,
                                   tid_name)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(.value = !!dplyr::sym(value_name), # a tidyeval trick
               value_norm    = .value - value_to_norm_by,
               value_normexp = 2^-value_norm) %>%
        dplyr::select(-.value) %>%
        return()
}

#' @describeIn calculate_deltacq_bysampleid Synonym for
#'   calculate_deltacq_plates.
#'
#' @export
#'
normalizeqPCR <- function(cq_df,
                          value_name = "cq",
                          norm_target_ids = "ALG9",
                          tid_name = "target_id") {
    lifecycle::deprecate_warn("0.2", "normalizeqPCR()",
                              "calculate_deltacq_bysampleid()",
        details = "Replaced with more descriptive name")
    calculate_deltacq_bysampleid(cq_df = cq_df,
                  norm_target_ids = norm_target_ids,
                  value_name = value_name,
                  tid_name = tid_name)
}
