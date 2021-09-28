#' Create a blank plate template as a tibble
#'
#' For more help, examples and explanations, see the plate setup vignette:
#' \code{vignette("platesetup_vignette", package = "tidyqpcr")}
#'
#' @param well_row Vector of Row labels, usually LETTERS
#' @param well_col Vector of Column labels, usually numbers
#' @return tibble (data frame) with columns well_row, well_col, well. This
#'   contains all pairwise combinations of well_row and well_col, as well as
#'   individual well names. Both well_row and well_col are coerced to factors
#'   (even if well_col is supplied as numbers), to ensure order is consistent.
#'
#'   However, well is a character vector as that is the default behaviour of
#'   "unite", and display order doesn't matter.
#'
#'   Default value describes a full 384-well plate.
#'
#' @examples
#' create_blank_plate(well_row=LETTERS[1:2],well_col=1:3)
#' create_blank_plate_96well()
#'
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr %>%
#' @importFrom forcats as_factor
#' @importFrom rlang .data
#'
create_blank_plate <- function(well_row = LETTERS[1:16], well_col = 1:24) {
    tidyr::crossing(well_row = as_factor(well_row),
                    well_col = as_factor(well_col)) %>%
        as_tibble() %>%
        tidyr::unite("well", .data$well_row, .data$well_col, 
                     sep = "", remove = FALSE)
}

#' @describeIn create_blank_plate create blank 96-well plate
#' @export
#'
create_blank_plate_96well <- function() {
    create_blank_plate(well_row = LETTERS[1:8], well_col = 1:12)
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
    well_row = make_row_names_lc1536(),
    well_col = 1:48) {
    create_blank_plate(well_row, well_col)
}

#' Create a 6-value, 24-column key for plates
#'
#' Create a 24-column key with 6 values repeated over 24 plate columns.
#' Each of the 6 values is repeated over 3x +RT Techreps and 1x -RT.
#'
#' @param ... Vectors of length 6 describing well contents,
#' e.g. sample_id or target_id
#' @return tibble (data frame) with 24 rows, and columns
#' well_col, prep_type, tech_rep, and supplied values.
#'
#' @examples
#' create_colkey_6_in_24(sample_id=LETTERS[1:6])
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble as_tibble

#' @importFrom tidyr %>%
#'
create_colkey_6_in_24 <- function(...) {
    colkey <- tibble(well_col   = factor(1:24),
                     prep_type    = c(rep("+RT", 18), rep("-RT", 6)) %>%
                         factor(levels = c("+RT", "-RT")),
                     tech_rep = rep(c(1, 2, 3, 1), each = 6) %>%
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
#' Creates a 24-column key for primer calibration, with 2x biol_reps and 2x
#' tech_reps, and 5-fold dilution until 5^4 of +RT; then -RT (no reverse
#' transcriptase), NT (no template) negative controls. That is a total of 6
#' versions of each sample replicate.
#'
#' @param dilution Numeric vector of length 6 describing sample dilutions
#' @param dilution_nice Character vector of length 6 with nice labels for sample
#'   dilutions
#' @param prep_type Character vector of length 6 describing type of sample (+RT,
#'   -RT, NT)
#' @param biol_rep Character vector of length 6 describing biological replicates
#' @param tech_rep Character vector of length 6 describing technical replicates
#' @return tibble (data frame) with 24 rows, and columns well_col, dilution,
#'   dilution_nice, prep_type, biol_rep, tech_rep.
#' @examples
#' create_colkey_4diln_2ctrl_in_24()
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble
#'
create_colkey_4diln_2ctrl_in_24 <- function(
                     dilution      = c(5 ^ (0:-3), 1, 1),
                     dilution_nice = c("1x", "5x", "25x", "125x", "-RT", "NT"),
                     prep_type         = c(rep("+RT", 4), "-RT", "NT"),
                     biol_rep       = rep(c("A", "B"), each = 12,
                                        length.out = 24),
                     tech_rep      = rep(1:2, each = 6,
                                        length.out = 24)
                     ) {
    tibble(well_col = factor(1:24),
           dilution = rep(dilution, 4),
           dilution_nice = rep(dilution_nice, 4),
           prep_type = factor(rep(prep_type, 4),
                         levels = c("+RT", "-RT", "NT")),
           biol_rep = factor(biol_rep),
           tech_rep = factor(tech_rep)
    )
}

#' Create a 6-dilution column key for primer calibration
#'
#' Creates a 24-column key for primer calibration, with 1x biol_reps and 3x
#' tech_reps, and 5-fold dilution until 5^6 of +RT; then -RT (no reverse
#' transcriptase), NT (no template) negative controls. That is a total of 8
#' versions of each replicate.
#'
#' @param dilution Numeric vector of length 8 describing sample dilutions
#' @param dilution_nice Character vector of length 8 with nice labels for sample
#'   dilutions
#' @param prep_type Character vector of length 8 describing type of sample (+RT,
#'   -RT, NT)
#' @param tech_rep Character vector of length 8 describing technical replicates
#' @return tibble (data frame) with 24 rows, and variables well_col, dilution,
#'   dilution_nice, prep_type, biol_rep, tech_rep.
#' @examples
#' create_colkey_6diln_2ctrl_in_24()
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble
#'
create_colkey_6diln_2ctrl_in_24 <- function(
                     dilution = c(5 ^ (0:-5), 1, 1),
                     dilution_nice = c("1x", "5x", "25x", "125x",
                                    "625x", "3125x", "-RT", "NT"),
                     prep_type=c(rep("+RT", 6), "-RT", "NT"),
                     tech_rep = rep(1:3, each = 8, length.out = 24)
                     ) {
    tibble(well_col = factor(1:24),
           dilution = rep(dilution, 3),
           dilution_nice = rep(dilution_nice, 3),
           prep_type = factor(rep(prep_type, 3),
                         levels = c("+RT", "-RT", "NT")),
           tech_rep = factor(tech_rep))
}

#' Create a 4-value, 16-row key for plates
#'
#' Create a 16-row key with 4 values repeated over 16 plate rows. Each of the 4
#' values is repeated over 3x +RT Techreps and 1x -RT.
#'
#' @param ... Vectors of length 4 describing well contents, e.g. sample_id or
#'   target_id
#' @return tibble (data frame) with 16 rows, and variables well_row, prep_type,
#'   tech_rep, and supplied values.
#' @examples
#' create_rowkey_4_in_16(sample_id=c("sheep","goat","cow","chicken"))
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr %>%
#'
create_rowkey_4_in_16 <- function(...) {
    rowkey <- tibble(well_row = factor(LETTERS[1:16]),
                     prep_type = factor(c(rep("+RT", 12), rep("-RT", 4)),
                                   levels = c("+RT", "-RT")),
                     tech_rep = factor(rep(c(1, 2, 3, 1), each = 4),
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
#' @return tibble (data frame) with 16 rows, and variables well_col, and
#'   supplied values.
#' @examples
#' create_rowkey_8_in_16_plain(sample_id=c("me","you","them","him",
#'                                    "her","dog","cat","monkey"))
#' @family plate creation functions
#'
#' @export
#' @importFrom tibble tibble
#' @importFrom tidyr %>%
#'
create_rowkey_8_in_16_plain <- function(...) {
    rowkey <- tibble(well_row = factor(LETTERS[1:16]))
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
#' @param plate tibble (data frame) with variables well_row, well_col, well.
#'   This would usually be produced by create_blank_plate(). It is possible to
#'   include other information in additional variables.
#' @param rowkey tibble (data frame) describing plate rows, with variables
#'   well_row and others.
#' @param colkey tibble (data frame) describing plate columns, with variables
#'   well_col and others.
#' @param coercefactors if TRUE, coerce well_row in rowkey and well_col in
#'   colkey to factors
#'
#' @return tibble (data frame) with variables well_row, well_col, well, and
#'   others.
#'
#'   This tibble contains all combinations of well_row and well_col found in the
#'   input plate, and all information supplied in rowkey and colkey distributed
#'   across every well of the plate. Return plate is ordered by row well_row
#'   then column well_col.
#'
#'   Note this ordering may cause a problem if well_col is supplied as a
#'   character (1,10,11,...), instead of a factor or integer (1,2,3,...). For
#'   this reason, the function by default converts well_row in `rowkey`, and
#'   well_col in `colkey`, to factors, taking factor levels from `plate`, and
#'   warns the user.
#'
#'   Other tidyqpcr functions require plate plans to contain variables
#'   sample_id, target_id, and prep_type, so `label_plate_rowcol` will warn if
#'   any of these are missing. This is a warning, not an error, because these
#'   variables can be added by users later.
#'
#' @examples
#' label_plate_rowcol(plate = create_blank_plate()) # returns blank plate
#' @family plate creation functions
#'
#' @export
#'
#' @importFrom rlang .data
#'
label_plate_rowcol <- function(plate,
                               rowkey = NULL,
                               colkey = NULL,
                               coercefactors = TRUE) {
    if (!is.null(colkey)) {
        assertthat::assert_that(assertthat::has_name(colkey, "well_col"))
        # Note: should this if clause be a freestanding function?
        # coerce_column_to_factor(df, col, warn=TRUE)?
        if (!is.factor(colkey$well_col) & coercefactors) {
            warning("coercing well_col to a factor with levels from plate$well_col")
            colkey <- dplyr::mutate(
                colkey,
                well_col = factor(.data$well_col,
                                  levels = levels(plate$well_col))
            )
        }
        plate <- dplyr::left_join(plate, colkey, by = "well_col")
    }
    if (!is.null(rowkey)) {
        assertthat::assert_that(assertthat::has_name(rowkey, "well_row"))
        if (!is.factor(rowkey$well_row) & coercefactors) {
            warning("coercing well_row to a factor with levels from plate$well_row")
            rowkey <- dplyr::mutate(
                rowkey,
                well_row = factor(.data$well_row,
                                  levels = levels(plate$well_row))
            )
        }
        plate <- dplyr::left_join(plate, rowkey, by = "well_row")
    }
    # check that plate contains sample_id, target_id, prep_type, warn if not
    if (! "sample_id" %in% names(plate)) {
        warning("plate does not contain variable sample_id")
    }
    if (! "target_id" %in% names(plate)) {
        warning("plate does not have variable target_id")
    }
    if (! "prep_type" %in% names(plate)) {
        warning("plate does not have variable prep_type")
    }
    return(dplyr::arrange(plate, .data$well_row, .data$well_col))
}


#' Display plate plan with sample_id, target_id, prep_type per well
#'
#' @param plate tibble with variables well_col, well_row, sample_id, target_id,
#'   prep_type. Output from label_plate_rowcol.
#'
#' @return ggplot object; major output is to plot it
#'
#' @examples # !Needs a labeled example from label_plate_rowcol...
#' @family plate creation functions
#'
#' @export
#' @importFrom forcats as_factor
#' @importFrom rlang .data
#'
display_plate <- function(plate) {
    rowlevels <- 
        dplyr::pull(plate, .data$well_row) %>%
        as_factor() %>%
        levels()
    #
    ggplot2::ggplot(data = plate,
                    ggplot2::aes(x = as_factor(.data$well_col),
                        y = as_factor(.data$well_row))) +
        ggplot2::geom_tile(ggplot2::aes(fill = .data$target_id), 
                           alpha = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = 
                                   paste(.data$target_id,
                                         .data$sample_id,
                                         .data$prep_type,
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
