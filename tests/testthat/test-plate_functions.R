# Units test for displaying a blank 1536 well

### Test functions give expected output ###

test_that("functions for displaying a blank 1536 well works", {
    blank_plate_schematic <- display_plate(create_blank_plate_1536well())
    
    expect_true(ggplot2::is.ggplot(blank_plate_schematic))
})

# Units test for creating, labeling and displaying a 384 well,
# primer calibration plate with 4 * dilution

### Create dataset of expected function outputs ###

simulated_colkey <- tibble(well_col = factor(1:24),
                      dilution = rep(c(1, 0.2, 0.04,
                                       0.008, 1, 1),
                                     times = 4),
                      dilution_nice = rep(c("1x", "5x", "25x",
                                            "125x", "-RT", "NT"),
                                          times = 4),
                      prep_type = factor(rep(c("+RT", "+RT",
                                               "+RT", "+RT",
                                               "-RT", "NT"),
                                             4),
                                         levels = c("+RT", "-RT", "NT")),
                      biol_rep = factor(rep(c("A", "B"),
                                            each = 12),
                                        levels = c("A", "B")),
                      tech_rep = factor(rep(rep(1:2, each = 6),
                                            times = 2))) %>%
    dplyr::mutate(sample_id = paste(biol_rep, dilution_nice, sep = "_"))

simulated_rowkey <- tibble(well_row = factor(LETTERS[1:16]),
                           target_id = c("Target_1", "Target_2", "Target_3", "Target_4",
                                         "Target_5", "Target_6", "Target_7", "Target_8",
                                         "Target_1", "Target_2", "Target_3", "Target_4",
                                         "Target_5", "Target_6", "Target_7", "Target_8"))

simulated_label_plate <- tibble(well = paste0(rep(LETTERS[1:16], each = 24),
                                               rep(1:24, times = 16)),
                           well_row = factor(rep(LETTERS[1:16],
                                                 each = 24)),
                           well_col = factor(rep(1:24, times = 16)),
                           dilution = rep(c(1, 0.2, 0.04,
                                            0.008, 1, 1),
                                          times = 64),
                           dilution_nice = rep(c("1x", "5x", "25x",
                                                 "125x", "-RT", "NT"),
                                               times = 64),
                           prep_type = factor(rep(c("+RT", "+RT", "+RT",
                                                    "+RT", "-RT", "NT"),
                                                  times = 64),
                                              levels = c("+RT", "-RT", "NT")),
                           biol_rep = factor(rep(rep(c("A", "B"),
                                                     each = 12),
                                                 times = 16),
                                             levels = c("A", "B")),
                           tech_rep = factor(rep(rep(1:2, each = 6),
                                                 times = 32)),
                           sample_id = paste(biol_rep, dilution_nice, sep = "_"),
                           target_id = rep(c("Target_1", "Target_2", "Target_3", "Target_4",
                                             "Target_5", "Target_6", "Target_7", "Target_8",
                                             "Target_1", "Target_2", "Target_3", "Target_4",
                                             "Target_5", "Target_6", "Target_7", "Target_8"),
                                           each = 24))

### Test functions give expected output ###

test_that("functions for creating a 384 well, primer calibration plate work with 4 * dilution", {
    calculated_blank_384_plate <- create_blank_plate()

    calculated_colkey <- create_colkey_4diln_2ctrl_in_24() %>%
        dplyr::mutate(sample_id = paste(biol_rep, dilution_nice, sep = "_"))

    calculated_rowkey <- create_rowkey_8_in_16_plain(target_id = c("Target_1", "Target_2",
                                                                   "Target_3", "Target_4",
                                                                   "Target_5", "Target_6",
                                                                   "Target_7", "Target_8"))

    calculated_labeled_384_plate <- label_plate_rowcol(calculated_blank_384_plate,
                                                       calculated_rowkey,
                                                       calculated_colkey)
    
    calculated_plate_empty <- display_plate(calculated_labeled_384_plate)
    
    calculated_plate_schematic <- display_plate_qpcr(calculated_labeled_384_plate)
    
    calculated_plate_value <- display_plate_value(calculated_labeled_384_plate, 
                                                  value = "dilution")

    expect_equal(calculated_colkey,
                 simulated_colkey)

    expect_equal(calculated_rowkey,
                 simulated_rowkey)

    expect_equal(calculated_labeled_384_plate,
                 simulated_label_plate)

    expect_true(ggplot2::is.ggplot(calculated_plate_empty))
    expect_true(ggplot2::is.ggplot(calculated_plate_schematic))
    expect_true(ggplot2::is.ggplot(calculated_plate_value))
})

# Unit test for creating, labeling and displaying a 384 well,
# primer calibration plate with 6 * dilution

### Create dataset of expected function outputs ###

simulated_colkey <- tibble(well_col = factor(1:24),
                           dilution = rep(c(1, 0.2, 0.04, 0.008,
                                            0.0016, 0.00032, 1, 1),
                                          times = 3),
                           dilution_nice = rep(c("1x", "5x", "25x", "125x",
                                                 "625x", "3125x", "-RT", "NT"),
                                               times = 3),
                           prep_type = factor(rep(c("+RT", "+RT", "+RT", "+RT",
                                                    "+RT", "+RT", "-RT", "NT"),
                                                  times = 3),
                                              levels = c("+RT", "-RT", "NT")),
                           tech_rep = factor(rep(rep(1:3, each = 8)))) %>%
    dplyr::mutate(sample_id = paste(dilution_nice, tech_rep, sep = "_"))

simulated_rowkey <- tibble(well_row = factor(LETTERS[1:16]),
                           target_id = c("Target_1", "Target_2", "Target_3", "Target_4",
                                         "Target_5", "Target_6", "Target_7", "Target_8",
                                         "Target_1", "Target_2", "Target_3", "Target_4",
                                         "Target_5", "Target_6", "Target_7", "Target_8"))

simulated_label_plate <- tibble(well = paste0(rep(LETTERS[1:16], each = 24),
                                               rep(1:24, times = 16)),
                                well_row = factor(rep(LETTERS[1:16], each = 24)),
                                well_col = factor(rep(1:24, times = 16)),
                                dilution = rep(c(1, 0.2, 0.04, 0.008,
                                                 0.0016, 0.00032, 1, 1),
                                               times = 48),
                                dilution_nice = rep(c("1x", "5x", "25x", "125x",
                                                      "625x", "3125x", "-RT", "NT"),
                                                    times = 48),
                                prep_type = factor(rep(c("+RT", "+RT", "+RT", "+RT",
                                                         "+RT", "+RT", "-RT", "NT"),
                                                       times = 48),
                                                   levels = c("+RT", "-RT", "NT")),
                                tech_rep = factor(rep(rep(1:3, each = 8),
                                                      times = 16)),
                                sample_id = paste(dilution_nice, tech_rep, sep = "_"),
                                target_id = rep(c("Target_1", "Target_2", "Target_3", "Target_4",
                                                  "Target_5", "Target_6", "Target_7", "Target_8",
                                                  "Target_1", "Target_2", "Target_3", "Target_4",
                                                  "Target_5", "Target_6", "Target_7", "Target_8"),
                                                each = 24))

### Test functions give expected output ###

test_that("functions for creating a 384 well, primer calibration plate work with 6 * dilution", {
    calculated_blank_384_plate <- create_blank_plate()

    calculated_colkey <- create_colkey_6diln_2ctrl_in_24() %>%
        dplyr::mutate(sample_id = paste(dilution_nice, tech_rep, sep = "_"))

    calculated_rowkey <- create_rowkey_8_in_16_plain(target_id = c("Target_1", "Target_2",
                                                                   "Target_3", "Target_4",
                                                                   "Target_5", "Target_6",
                                                                   "Target_7", "Target_8"))

    calculated_labeled_384_plate <- label_plate_rowcol(calculated_blank_384_plate,
                                                       calculated_rowkey,
                                                       calculated_colkey)

    expect_equal(calculated_colkey,
                 simulated_colkey)

    expect_equal(calculated_rowkey,
                 simulated_rowkey)

    expect_equal(calculated_labeled_384_plate,
                 simulated_label_plate)
})

# Unit test for creating standard 24 by 16 384 well plate
# with tech reps across rows

### Create dataset of expected function outputs ###

# colkey has well_col a number not a factor, testing coercion
simulated_colkey_6_in_24 <- 
    tibble(well_col = 1:24,
           biol_rep = factor(rep(c("A", "B",
                                   "C", "D"),
                                 each = 6),
                             levels = c("A", "B",
                                        "C", "D")),
           sample_id = rep(c("Sample_1", "Sample_2",
                             "Sample_3", "Sample_4",
                             "Sample_5", "Sample_6"),
                           times = 4))


simulated_rowkey_4_in_16 <- 
    tibble(well_row = factor(LETTERS[1:16]),
           prep_type = factor(rep(c("+RT", "+RT",
                                    "+RT", "-RT"),
                                  each = 4),
                              levels = c("+RT", "-RT")),
           tech_rep = factor(rep(c(1, 2, 3, 1),
                                 each = 4),
                             levels = c(1:3)),
           target_id = rep(c("Target_1", "Target_2",
                             "Target_3", "Target_4"),
                           times = 4))

simulated_labeled_384_plate_row4_in_16 <- 
    tibble(well = paste0(rep(LETTERS[1:16], each = 24),
                         rep(1:24, times = 16)),
           well_row = factor(rep(LETTERS[1:16], each = 24)),
           well_col = factor(rep(1:24, times = 16)),
           biol_rep = factor(rep(rep(c("A", "B",
                                       "C", "D"),
                                     each = 6),
                                 times = 16),
                             levels = c("A", "B",
                                        "C", "D")),
           sample_id = rep(c("Sample_1", "Sample_2",
                             "Sample_3", "Sample_4",
                             "Sample_5", "Sample_6"),
                           times = 64),
           prep_type = factor(rep(c("+RT", "+RT",
                                    "+RT", "-RT"),
                                  each = 96),
                              levels = c("+RT", "-RT")),
           tech_rep = factor(c(rep(1:3, each = 96),
                               rep(1, times = 96))),
           target_id = rep(rep(c("Target_1", "Target_2",
                                 "Target_3", "Target_4"),
                               each = 24),
                           times = 4))

### Test functions give expected output ###

test_that("functions for creating a standard 384 well with tech reps across rows", {
    target_id_levels <- c("Target_1", "Target_2",
                          "Target_3", "Target_4")
    calculated_rowkey_4_in_16 <- 
        create_rowkey_4_in_16(target_id = target_id_levels)

    calculated_labeled_384_plate_row4_in_16 <- 
        label_plate_rowcol(create_blank_plate(),
                           calculated_rowkey_4_in_16,
                           simulated_colkey_6_in_24)

    expect_equal(calculated_rowkey_4_in_16,
                 simulated_rowkey_4_in_16)

    expect_equal(calculated_labeled_384_plate_row4_in_16,
                 simulated_labeled_384_plate_row4_in_16)

})

# Unit test for creating standard 24 by 16 384 well plate
# with tech reps across columns

### Create dataset of expected function outputs ###

simulated_colkey_6t_in_24 <- 
    tibble(well_col = factor(1:24),
           prep_type = factor(rep(c("+RT", "+RT",
                                    "+RT", "-RT"),
                                  each = 6),
                              levels = c("+RT", "-RT")),
           tech_rep = factor(rep(c(1, 2, 3, 1),
                                 each = 6),
                             levels = c(1:3)),
           sample_id = rep(c("Sample_1", "Sample_2",
                             "Sample_3", "Sample_4",
                             "Sample_5", "Sample_6"),
                           times = 4))

# rowkey has well_row a character not a factor, testing coercion
simulated_rowkey_T16 <- 
    tibble(well_row = LETTERS[1:16],
           target_id = c("Target_1", "Target_2", "Target_3", "Target_4",
                         "Target_5", "Target_6", "Target_7", "Target_8",
                         "Target_1", "Target_2", "Target_3", "Target_4",
                         "Target_5", "Target_6", "Target_7", "Target_8"))

simulated_labeled_384_plate_col6t_in_24 <- 
    tibble(well = paste0(rep(LETTERS[1:16],
                             each = 24),
                         rep(1:24, times = 16)),
           well_row = factor(rep(LETTERS[1:16], each = 24)),
           well_col = factor(rep(1:24, times = 16)),
           prep_type = factor(rep(rep(c("+RT", "+RT",
                                        "+RT", "-RT"),
                                      each = 6),
                                  times = 16),
                              levels = c("+RT", "-RT")),
           tech_rep = factor(rep(rep(c(1:3, 1),
                                     each = 6),
                                 times = 16),
                             levels = c(1, 2, 3)),
           sample_id = rep(c("Sample_1", "Sample_2", "Sample_3",
                             "Sample_4", "Sample_5", "Sample_6"),
                           times = 64),
           target_id = rep(c("Target_1", "Target_2", "Target_3", "Target_4",
                             "Target_5", "Target_6", "Target_7", "Target_8",
                             "Target_1", "Target_2", "Target_3", "Target_4",
                             "Target_5", "Target_6", "Target_7", "Target_8"),
                           each = 24))

### Test functions give expected output ###

test_that("functions for creating a standard 384 well with tech reps across columns", {
    sample_id_levels <- c("Sample_1", "Sample_2", "Sample_3",
                          "Sample_4", "Sample_5", "Sample_6")

    calculated_colkey_6t_in_24 <- 
        create_colkey_6_in_24(sample_id = sample_id_levels)

    calculated_labeled_384_plate_col6t_in_24 <- 
        label_plate_rowcol(create_blank_plate(),
                           simulated_rowkey_T16,
                           calculated_colkey_6t_in_24)

    expect_equal(calculated_colkey_6t_in_24,
                 simulated_colkey_6t_in_24)

    expect_equal(calculated_labeled_384_plate_col6t_in_24,
                 simulated_labeled_384_plate_col6t_in_24)
})
