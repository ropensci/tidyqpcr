# Units test for creating, labeling and displaying a 384 well, primer calibration plate 

test_colkey <- tibble(well_col = factor(1:24),
                      dilution = rep(c(1,0.2,0.04,0.008,1,1), times = 4),
                      dilution_nice = rep(c("1x", "5x", "25x", "125x", "-RT", "NT"), times = 4),
                      prep_type = factor(rep(c("+RT", 
                                               "+RT", 
                                               "+RT", 
                                               "+RT",
                                               "-RT",
                                               "NT"), 4),
                                         levels = c("+RT", "-RT", "NT")),
                      biol_rep = factor(rep(c("A", "B"), each = 12), levels = c("A", "B")),
                      tech_rep = factor(rep(rep(1:2, each = 6), times = 2)))

test_rowkey <- tibble(well_row = LETTERS[1:16], target_id = c("RPS3_ORF",
                                                     "mCherry_ORF",
                                                     "mTurquoies_ORF",
                                                     "PGK1_ORF",
                                                     "SRO9_ORF",
                                                     "TSA1_ORF",
                                                     "PIR1_ORF",
                                                     "HSP30_ORF",
                                                     "RPS3_ORF",
                                                     "mCherry_ORF",
                                                     "mTurquoies_ORF",
                                                     "PGK1_ORF",
                                                     "SRO9_ORF",
                                                     "TSA1_ORF",
                                                     "PIR1_ORF",
                                                     "HSP30_ORF"))

test_label_plate <- tibble(well = paste0( rep(LETTERS[1:16], each = 24), rep(1:24, times = 16)),
                           well_row = factor(rep(LETTERS[1:16], each = 24)),
                           well_col = factor(rep(1:24, times = 16)),
                           dilution = rep(c(1,0.2,0.04,0.008,1,1), times = 64),
                           dilution_nice = rep(c("1x", "5x", "25x", "125x", "-RT", "NT"), times = 64),
                           prep_type = factor(rep(c("+RT", 
                                                    "+RT", 
                                                    "+RT", 
                                                    "+RT",
                                                    "-RT",
                                                    "NT"), times = 64),
                                              levels = c("+RT", "-RT", "NT")),
                           biol_rep = factor(rep(rep(c("A", "B"), each = 12), times = 16), levels = c("A", "B")),
                           tech_rep = factor(rep(rep(1:2, each = 6), times =32)),
                           target_id = rep(c("RPS3_ORF",
                                             "mCherry_ORF",
                                             "mTurquoies_ORF",
                                             "PGK1_ORF",
                                             "SRO9_ORF",
                                             "TSA1_ORF",
                                             "PIR1_ORF",
                                             "HSP30_ORF",
                                             "RPS3_ORF",
                                             "mCherry_ORF",
                                             "mTurquoies_ORF",
                                             "PGK1_ORF",
                                             "SRO9_ORF",
                                             "TSA1_ORF",
                                             "PIR1_ORF",
                                             "HSP30_ORF"), each = 24))

test_that("functions for creating a 384 well, primer calibration plate work", {
    blank_384_plate <- create_blank_plate()
    
    colkey <- create_colkey_4diln_2ctrl_in_24()
    
    rowkey <- create_rowkey_8_in_16_plain(target_id = c("RPS3_ORF",
                                                        "mCherry_ORF",
                                                        "mTurquoies_ORF",
                                                        "PGK1_ORF",
                                                        "SRO9_ORF",
                                                        "TSA1_ORF",
                                                        "PIR1_ORF",
                                                        "HSP30_ORF"))
    
    labeled_384_plate <- label_plate_rowcol(blank_384_plate, rowkey, colkey)
    
    plate_schematic <- display_plate(labeled_384_plate)
    
    expect_equal(colkey, test_colkey)
    
    expect_equal(rowkey, test_rowkey)
    
    expect_equal(labeled_384_plate, test_label_plate)
    
    expect_true(is.ggplot(plate_schematic))
})
