
#' Create a blank plate template as a tibble
#'
#' For more help, examples and explanations, see the plate setup vignette:
#' \code{vignette("platesetup_vignette", package = "tidyqpcr")}
#'
#' @param WellR Vector of Row labels, usually LETTERS
#' @param WellC Vector of Column labels, usually numbers
#' @return tibble (data frame) with columns WellR, WellC, Well. This contains
#'   all pairwise combinations of WellR and WellC, as well as individual Well
#'   names. Both WellR and WellC are coerced to factors (even if WellC
#'   is supplied as numbers), to ensure order is consistent.
#'
#'   However, Well is a character vector as that is the default behaviour of
#'   "unite", and display order doesn't matter.
#'
#'   Default value describes a full 384-well plate.
#'
#' @examples
#' create_blank_plate(WellR=LETTERS[1:2],WellC=1:3)
#' create_blank_plate_96well()
#'
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom forcats as_factor
#'
create_blank_plate <- function(WellR = LETTERS[1:16], WellC = 1:24) {
    tidyr::crossing(WellR = as_factor(WellR),
                             WellC = as_factor(WellC)) %>%
        as_tibble() %>%
        tidyr::unite(Well, WellR, WellC, sep = "", remove = FALSE)
}

#' @describeIn create_blank_plate create blank 96-well plate
#' @export
#'
create_blank_plate_96well <- function() {
    create_blank_plate(WellR = LETTERS[1:8], WellC = 1:12)
}

#' @describeIn create_blank_plate Row names for 1536-well plates on
#' Roche Lightcycler (tm) 1536 Aa,Ab,Ac,Ad,Ba,...,Hd.
#' @export
#'
make_row_names_lc1536 <- function() {
    paste0(rep(LETTERS[1:8], each = 4), letters[1:4])
}

#' @describeIn create_blank_plate Row names for 1536-well plates on
#' Labcyte Echo A,B,...,Z,AA,AB,...,AF.
#' @export
#'
make_row_names_echo1536 <- function() {
    c(LETTERS[1:26], paste0("A", LETTERS[1:6]))
}

#' @describeIn create_blank_plate create blank 1536-well plate
#' @export
#'
create_blank_plate_1536well <- function(
    WellR = make_row_names_lc1536(),
    WellC = 1:48) {
    create_blank_plate(WellR, WellC)
}

#' Create a 6-value, 24-column key for plates
#'
#' Create a 24-column key with 6 values repeated over 24 plate columns.
#' Each of the 6 values is repeated over 3x +RT Techreps and 1x -RT.
#'
#' @param ... Vectors of length 6 describing well contents,
#' e.g. sample or probe.
#' @return tibble (data frame) with 24 rows, and columns
#' WellC, Type, TechRep, and supplied values.
#'
#' @examples
#' create_colkey_6in24(sample_id=LETTERS[1:6])
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom magrittr %>%
#'
create_colkey_6in24 <- function(...) {
    colkey <- tibble(WellC   = factor(1:24),
                     Type    = c(rep("+RT", 18), rep("-RT", 6)) %>%
                         factor(levels = c("+RT", "-RT")),
                     TechRep = rep(c(1, 2, 3, 1), each = 6) %>%
                         factor(levels = 1:3)
                     )
    if (!missing(...)) {
        pieces6 <- list(...) %>% as_tibble()
        stopifnot(nrow(pieces6) == 6)
        pieces24 <- dplyr::bind_rows(pieces6, pieces6, pieces6, pieces6)
        colkey <- dplyr::bind_cols(colkey, pieces24)
    }
    return(colkey)
}

#' Create a 4-dilution column key for primer calibration
#'
#' Creates a 24-column key for primer calibration, with 2x BioReps and 2x
#' TechReps, and 5-fold dilution until 5^4 of +RT; then -RT (no reverse
#' transcriptase), NT (no template) negative controls. That is a total of 6
#' versions of each sample replicate.
#'
#' @param Dilution Numeric vector of length 6 describing sample dilutions
#' @param DilutionNice Character vector of length 6 with nice labels for sample
#'   dilutions
#' @param Type Character vector of length 6 describing type of sample (+RT, -RT,
#'   NT)
#' @param BioRep Character vector of length 6 describing biological replicates
#' @param TechRep Character vector of length 6 describing technical replicates
#' @return tibble (data frame) with 24 rows, and columns WellC, Dilution,
#'   DilutionNice, Type, BioRep, TechRep.
#' @examples
#' create_colkey_4diln_2ctrl_in24()
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble
#'
create_colkey_4diln_2ctrl_in24 <- function(
                     Dilution     = c(5 ^ (0:-3), 1, 1),
                     DilutionNice = c("1x", "5x", "25x", "125x", "-RT", "NT"),
                     Type         = c(rep("+RT", 4), "-RT", "NT"),
                     BioRep       = rep(c("A", "B"), each = 12,
                                        length.out = 24),
                     TechRep      = rep(1:2, each = 6,
                                        length.out = 24)
                     ) {
    tibble(WellC = factor(1:24),
           Dilution = rep(Dilution, 4),
           DilutionNice = rep(DilutionNice, 4),
           Type = factor(rep(Type, 4),
                         levels = c("+RT", "-RT", "NT")),
           BioRep = factor(BioRep),
           TechRep = factor(TechRep)
    )
}

#' Create a 6-dilution column key for primer calibration
#'
#' Creates a 24-column key for primer calibration, with 1x BioReps and 3x
#' TechReps, and 5-fold dilution until 5^6 of +RT; then -RT (no reverse
#' transcriptase), NT (no template) negative controls. That is a total of 8
#' versions of each replicate.
#'
#' @param Dilution Numeric vector of length 8 describing sample dilutions
#' @param DilutionNice Character vector of length 8 with nice labels for sample
#'   dilutions
#' @param Type Character vector of length 8 describing type of sample (+RT, -RT,
#'   NT)
#' @param TechRep Character vector of length 8 describing technical replicates
#' @return tibble (data frame) with 24 rows, and variables WellC, Dilution,
#'   DilutionNice, Type, BioRep, TechRep.
#' @examples
#' create_colkey_6diln_2ctrl_in24()
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble
#'
create_colkey_6diln_2ctrl_in24 <- function(
                     Dilution=c(5 ^ (0:-5), 1, 1),
                     DilutionNice=c("1x", "5x", "25x", "125x",
                                    "625x", "3125x", "-RT", "NT"),
                     Type=c(rep("+RT", 6), "-RT", "NT"),
                     TechRep = rep(1:3, each = 8, length.out = 24)
                     ) {
    tibble(WellC = factor(1:24),
           Dilution = rep(Dilution, 3),
           DilutionNice = rep(DilutionNice, 3),
           Type = factor(rep(Type, 3),
                         levels = c("+RT", "-RT", "NT")),
           TechRep = factor(TechRep))
}

#' Create a 4-value, 16-row key for plates
#'
#' Create a 16-row key with 4 values repeated over 16 plate rows. Each of the 4
#' values is repeated over 3x +RT Techreps and 1x -RT.
#'
#' @param ... Vectors of length 4 describing well contents, e.g. sample or
#'   probe.
#' @return tibble (data frame) with 16 rows, and variables WellR, Type, TechRep,
#'   and supplied values.
#' @examples
#' create_rowkey_4in16(sample_id=c("sheep","goat","cow","chicken"))
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom magrittr %>%
#'
create_rowkey_4in16 <- function(...) {
    rowkey <- tibble(WellR = LETTERS[1:16],
                     Type = factor(c(rep("+RT", 12), rep("-RT", 4)),
                                   levels = c("+RT", "-RT")),
                     TechRep = factor(rep(c(1, 2, 3, 1), each = 4),
                                      levels = 1:3))
    if (!missing(...)) {
        pieces4 <- list(...) %>% as_tibble()
        stopifnot(nrow(pieces4) == 4)
        pieces16 <- dplyr::bind_rows(pieces4, pieces4, pieces4, pieces4)
        rowkey <- dplyr::bind_cols(rowkey, pieces16)
    }
    return(rowkey)
}

#' Create a plain 8-value, 16-row key for plates
#'
#' Create a 16-row key with 8 values repeated over 16 plate rows. No other
#' information is included by default, hence "plain".
#'
#' @param ... Vectors of length 8 describing well contents, e.g. sample or
#'   probe.
#' @return tibble (data frame) with 16 rows, and variables WellC, and supplied
#'   values.
#' @examples
#' create_rowkey_8in16_plain(sample_id=c("me","you","them","him",
#'                                    "her","dog","cat","monkey"))
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#'
create_rowkey_8in16_plain <- function(...) {
    rowkey <- tibble(WellR = LETTERS[1:16])
    if (!missing(...)) {
        pieces8 <- list(...) %>% as_tibble()
        stopifnot(nrow(pieces8) == 8)
        pieces16 <- dplyr::bind_rows(pieces8, pieces8)
        rowkey <- dplyr::bind_cols(rowkey, pieces16)
    }
    return(rowkey)
}

#' Label a plate with sample and probe information
#'
#' See vignettes for further examples
#'
#' @param plate tibble (data frame) with variables WellR, WellC, Well. This
#'   would usually be produced by create_blank_plate(). It is possible to
#'   include other information in additional variables.
#' @param rowkey tibble (data frame) describing plate rows, with variables WellR
#'   and others.
#' @param colkey tibble (data frame) describing plate columns, with variables
#'   WellC and others.
#' @param coercefactors if TRUE, coerce WellR in rowkey and WellC in colkey
#' to factors
#'
#' @return tibble (data frame) with variables WellR, WellC, Well, and others.
#'
#'   This tibble contains all combinations of WellR and WellC found in the input
#'   plate, and all information supplied in rowkey and colkey distributed across
#'   every well of the plate. Return plate is ordered by row WellR then column
#'   WellC.
#'
#'   Note this ordering may cause a problem if WellC is supplied as a character
#'   (1,10,11,...), instead of a factor or integer (1,2,3,...). For this reason,
#'   the function by default converts WellR in `rowkey`, and WellC in `colkey`,
#'   to factors, taking factor levels from `plate`, and warns the user.
#'
#'   Other tidyqpcr functions require plate plans to contain variables
#'   sample_id, target_id, and Type, so `label_plate_rowcol` will warn if any of
#'   these are missing. This is a warning, not an error, because these variables
#'   can be added by users later.
#'
#' @examples
#' label_plate_rowcol(plate = create_blank_plate()) # returns blank plate
#' @family plate creation functions
#'
#' @export
#'
label_plate_rowcol <- function(plate,
                               rowkey = NULL,
                               colkey = NULL,
                               coercefactors = TRUE) {
    if (!is.null(colkey)) {
        assertthat::assert_that(assertthat::has_name(colkey, "WellC"))
        # Note: should this if clause be a freestanding function?
        # coerce_column_to_factor(df, col, warn=TRUE)?
        if (!is.factor(colkey$WellC) & coercefactors) {
            warning("coercing WellC to a factor with levels from plate$WellC")
            colkey <- dplyr::mutate(colkey,
                                    WellC = factor(WellC,
                                                   levels = levels(plate$WellC))
                                    )
        }
        plate <- dplyr::left_join(plate, colkey, by = "WellC")
    }
    if (!is.null(rowkey)) {
        assertthat::assert_that(assertthat::has_name(rowkey, "WellR"))
        if (!is.factor(rowkey$WellR) & coercefactors) {
            warning("coercing WellR to a factor with levels from plate$WellR")
            rowkey <- dplyr::mutate(rowkey,
                                    WellR = factor(WellR,
                                                   levels = levels(plate$WellR))
                                    )
        }
        plate <- dplyr::left_join(plate, rowkey, by = "WellR")
    }
    # check that plate contains sample_id, target_id, Type, warn if not
    if (! "sample_id" %in% names(plate)) {
        warning("plate does not contain variable sample_id")
    }
    if (! "target_id" %in% names(plate)) {
        warning("plate does not have variable target_id")
    }
    if (! "Type" %in% names(plate)) {
        warning("plate does not have variable Type")
    }
    return(dplyr::arrange(plate, WellR, WellC))
}


#' Display plate plan with sample_id, target_id, Type per Well
#'
#' @param plate tibble with variables WellC, WellR, sample_id, target_id, Type.
#'   Output from label_plate_rowcol.
#'
#' @return ggplot object; major output is to plot it
#'
#' @examples # !Needs a labeled example from label_plate_rowcol...
#' @family plate creation functions
#'
#' @export
#' @importFrom forcats as_factor
#'
display_plate <- function(plate) {
    rowlevels <- plate %>%
        dplyr::pull(WellR) %>%
        as_factor() %>%
        levels()
    #
    ggplot2::ggplot(data = plate,
                    aes(x = as_factor(WellC),
                        y = as_factor(WellR))) +
        ggplot2::geom_tile(aes(fill = target_id), alpha = 0.3) +
        ggplot2::geom_text(aes(label = paste(target_id,
                                             sample_id,
                                             Type,
                                             sep = "\n")),
                           size = 2.5, lineheight = 1) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0),
                                  limits = rev(rowlevels)) +
        ggplot2::coord_equal() +
        ggplot2::theme_void() +
        ggplot2::theme(axis.text = ggplot2::element_text(angle = 0),
                       panel.grid.major = ggplot2::element_blank(),
                       legend.position = "none",
                       plot.margin = grid::unit(rep(0.01, 4), "npc"),
                       panel.border = ggplot2::element_blank())
}
