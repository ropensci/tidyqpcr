# Unit test for the deltacq and deltadeltacq calculation functions

### Create dataset of expected function outputs ###

simulated_48_well_plate_plan <- create_blank_plate_96well() %>%
    dplyr::filter(well_row %in% c("A", "B", "C", "D")) %>%
    dplyr::mutate(target_id = rep(c("Target_1", "Target_2",
                                    "Target_3", "Target_4"),
                                  each = 12),
                  sample_id = rep(rep(c("Sample_1", "Sample_2", "Sample_3"),
                                      each = 4),
                                  times = 4),
                  tech_rep = rep(c(1, 2, 3, 1),
                                 times = 12),
                  prep_type = rep(c("+RT", "+RT",
                                    "+RT", "-RT"),
                                  times = 12))

simulated_48_well_plate_with_cq <- simulated_48_well_plate_plan %>%
    dplyr::mutate(cq = c(8, 8, 9, 20, 10, 11, 11, 20, 14, 14, 14, 20,
                         6, 6, 6, 19, 12, 11, 13, 20, 15, 15, 17, 20,
                         5, 6, 7, 20, 10, 10, 10, 19, 10, 9, 8, 20,
                         7, 7, 11, 20, 9, 10, 9, 20, 9, 5, 9, 20))

simulated_48_well_plate_with_deltacq <- simulated_48_well_plate_with_cq %>%
    dplyr::filter(prep_type == "+RT") %>%
    dplyr::mutate(ref_cq = rep(c(6, 6, 6, 10, 10, 10, 9, 9, 9), times = 4),
                  delta_cq = c(2, 2, 3, 0, 1, 1, 5, 5, 5,
                               0, 0, 0, 2, 1, 3, 6, 6, 8,
                               -1, 0, 1, 0, 0, 0, 1, 0, -1,
                               1, 1, 5, -1, 0, -1, 0, -4, 0),
                  rel_abund = 2^-delta_cq)

simulated_48_well_plate_with_deltadeltacq_positive <- simulated_48_well_plate_with_deltacq %>%
    dplyr::mutate(ref_delta_cq = rep(c(2, 2, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1), each = 3),
                  deltadelta_cq = c(0, 0, -1, 2, 1, 1, -3, -3, -3,
                                    0, 0, 0, -2, -1, -3, -6, -6, -8,
                                    1, 0, -1, 0, 0, 0, -1, 0, 1,
                                    0, 0, -4, 2, 1, 2, 1, 5, 1),
                  fold_change = 2^deltadelta_cq)

simulated_48_well_plate_with_deltadeltacq_negative <- simulated_48_well_plate_with_deltacq %>%
    dplyr::mutate(ref_delta_cq = rep(c(2, 2, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1), each = 3),
                  deltadelta_cq = c(0, 0, 1, -2, -1, -1, 3, 3, 3,
                                    0, 0, 0, 2, 1, 3, 6, 6, 8,
                                    -1, 0, 1, 0, 0, 0, 1, 0, -1,
                                    0, 0, 4, -2, -1, -2, -1, -5, -1),
                  fold_change = 2^deltadelta_cq)

### Test functions give expected output ###

test_that("Unit test for the calculate_deltacq and calculate_deltadeltacq functions", {
    calculated_48_well_plate_with_deltacq <- 
        calculate_deltacq_bysampleid(simulated_48_well_plate_with_cq %>%
                                         dplyr::filter(prep_type == "+RT"), 
                                     ref_target_ids = "Target_3") %>%
        dplyr::arrange(well_row, well_col)
    
    calculated_48_well_plate_with_deltadeltacq_positive <- 
        calculate_deltadeltacq_bytargetid(calculated_48_well_plate_with_deltacq,
                                          ref_sample_ids = "Sample_1") %>%
        dplyr::arrange(well_row, well_col)
    
    calculated_48_well_plate_with_deltadeltacq_negative <- 
        calculate_deltadeltacq_bytargetid(calculated_48_well_plate_with_deltacq,
                                          ref_sample_ids = "Sample_1", 
                                          ddcq_positive = FALSE) %>%
        dplyr::arrange(well_row, well_col)

    expect_equal(calculated_48_well_plate_with_deltacq,
                 simulated_48_well_plate_with_deltacq)
    expect_equal(calculated_48_well_plate_with_deltadeltacq_negative,
                 simulated_48_well_plate_with_deltadeltacq_negative)
    expect_equal(calculated_48_well_plate_with_deltadeltacq_positive,
                 simulated_48_well_plate_with_deltadeltacq_positive)
})
