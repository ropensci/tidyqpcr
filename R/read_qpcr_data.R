#' Reads raw text-format fluorescence data in 1 colour from Roche Lightcyclers
#' 
#' This is the data from "export in text format" from the Lightcycler software.
#' The other data format, .ixo, can be converted to .txt format by the
#' Lightcycler software.
#'
#' This function is a thin wrapper around readr::read_tsv.
#'
#' @param filename file name
#' @param skip number of lines to skip, defaults to 2
#' @param col_names names to give to columns
#' @param col_types data types of columns
#' @param ... other arguments to pass to read_tsv, if needed
#'
#' @return tibble containing raw data, with default column names:
#'
#'   well: the well of the plate, e.g. A1
#'
#'   sample_info: this is the "Sample" field entered in lightcycler software,
#'   defaults to "Sample X"
#'
#'   program_no: the number of the cycler program, for 2-step PCR defaults to 1
#'   = melt, 2 = amplify, 3 = melt.
#'
#'   segment_no: the number of the segment of the cycler program, e.g.
#'   hold/raise/lower temperature
#'
#'   cycle: the cycle number, for programs with repeated cycles (i.e.
#'   amplification)
#'
#'   time: the time of fluorescence reading acquisition (in what units???)
#'
#'   temperature: the temperature of the block at fluorescence acquisition
#'
#'   fluor_raw: the raw fluorescence reading in "arbitrary units". For SYBR safe, this
#'   would be 483nm excitation, 533nm emission.
#'   
#' @export
#' @seealso read_lightcycler_1colour_cq
#'
#' @examples read_lightcycler_1colour_raw(system.file("extdata/Edward_qPCR_Nrd1_calibration_2019-02-02.txt.gz", 
#'                                                  package = "tidyqpcr"))
#'
read_lightcycler_1colour_raw <- function(
    filename,
    skip = 2,
    col_names = c(
        "well", "sample_info", "program_no", "segment_no",
        "cycle", "time", "temperature", "fluor_raw"
    ), 
    col_types = "ccffinnn",
    ...) {
    readr::read_tsv(file = filename,
                    skip = skip,
                    col_names = col_names,
                    col_types = col_types,
                    ...)
}

#' Reads quantification cycle (cq) data in 1 colour from Roche Lightcyclers
#'
#' This is the data from "export in text format" from the analysis tab in the
#' Lightcycler software. That software calls cq, "Cp".
#'
#' This function is a thin wrapper around readr::read_tsv.
#'
#' @param filename file name
#' @param skip number of lines to skip, defaults to 2
#' @param col_names names to give to columns
#' @param col_types data types of columns
#' @param ... other arguments to pass to read_tsv, if needed
#'
#' @return tibble containing cq data
#' @export
#' @seealso read_lightcycler_1colour_raw
#'
#' @examples read_lightcycler_1colour_cq(system.file("extdata/Edward_qPCR_Nrd1_calibration_2019-02-02_Cq.txt.gz", 
#'                                                  package = "tidyqpcr"))
#'
read_lightcycler_1colour_cq <- function(
    filename,
    skip = 2,
    col_names = c(
        "include", "color", "well", "sample_info",
        "cq", "concentration", "standard", "status"
    ), 
    col_types = "liccddil",
    ...) {
    readr::read_tsv(file = filename,
                    skip = 2,
                    col_names=col_names,
                    col_types=col_types,
                    ...)
}