
#' @describeIn calculate_deltacq_bysampleid get the median value of a set of normalization
#'   (reference) probes, for a single sample.
#'
#' @param normby.function Function to use to calculate the value to
#' normalise by on log2/Cq scale.
#' Default value is median, alternatively could use mean.
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom stats median
#'
calculate_normcq <- function(cq_df,
                      value = "Cq",
                      normTargetIDs = "ALG9",
                      probename = "TargetID",
                      normby.function = median) {
    # make subset of cq_df where gene is one or more normTargetIDs
    norm.by <- dplyr::filter(cq_df,
                             !!dplyr::sym(probename) %in% normTargetIDs) %>%
        .[[value]] %>%
        normby.function(na.rm = TRUE)
    #
    # assign summary (median) value to cq_df$norm.by
    # note this is the same value for every row, a waste of space technically
    dplyr::mutate(cq_df, norm.by = norm.by)
}

#' Calculate delta Cq to normalize quantification cycle (log2-fold) data within
#' SampleID.
#'
#' @param cq_df a data frame containing columns `SampleID`, value (default `Cq`)
#'   and probename (default `TargetID`). Crucially, SampleID should be the same
#'   for different technical replicates measuring identical reactions in
#'   different wells of the plate, but differ for different biological and
#'   experimental replicates.
#' @param value the column name of the value that will be normalized
#' @param normTargetIDs names of PCR probes (or primer sets) to normalize by,
#'   i.e. reference genes
#' @param probename the column name for probe sets
#'
#' @return data frame like cq_df with three additional columns:
#'
#'   \tabular{ll}{ norm.by       \tab the median value of the reference probes
#'   \cr Value.norm    \tab the normalized value, \eqn{\Delta Cq} \cr
#'   Value.normexp \tab the normalized ratio, \eqn{2^(-\Delta Cq)} }
#'
#' @export
#' @importFrom magrittr %>%
#'   
calculate_deltacq_bysampleid <- function(cq_df,
                                         normTargetIDs,
                                         value = "Cq",
                                         probename = "TargetID") {
    cq_df %>%
        dplyr::group_by(SampleID) %>%
        dplyr::do(calculate_normcq(., value, normTargetIDs, probename)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(.Value = !!dplyr::sym(value), # a tidyeval trick
               Value.norm    = .Value - norm.by,
               Value.normexp = 2^-Value.norm) %>%
        dplyr::select(-.Value) %>%
        return()
}

#' @describeIn calculate_deltacq_bysampleid Synonym for calculate_deltacq_plates.
#'
#' @export
#'
normalizeqPCR <- function(cq_df,
                          value = "Cq",
                          normTargetIDs = "ALG9",
                          probename = "TargetID") {
    lifecycle::deprecate_warn("0.2", "normalizeqPCR()",
                              "calculate_deltacq_bysampleid()",
        details = "Replaced with more descriptive name")
    calculate_deltacq_bysampleid(cq_df = cq_df,
                  normTargetIDs = normTargetIDs,
                  value = value,
                  probename = probename)
}
