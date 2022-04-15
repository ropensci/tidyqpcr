# Unit test for the debaseline function on amplification curve data

### Create dataset of expected function outputs ###

plate_wells <- expand.grid(well_row = c("A", "B", "C", "D"), well_col = 1:4) %>%
    tidyr::unite(well, c(well_row, well_col), sep = "") %>%
    dplyr::pull(well)

simulated_amplification_curve_data <- 
    tibble(well = rep(plate_wells, times = 15),
           cycle = rep(1:15, each = 16),
           program_no = 2,
           fluor_raw = rep(c(10, 15, 13, 20, 20, 50, 60, 70,
                             90, 130, 170, 280, 440, 670, 1020),
                           each = 16),
           fluor_base = rep(35, 240))

simulated_debased_amplification_curve_data <- simulated_amplification_curve_data %>%
    dplyr::mutate(fluor_signal = rep(c(-25, -20, -22, -15, -15, 15, 25, 35,
                                       55, 95, 135, 245, 405, 635, 985),
                                     each = 16))

### Test functions give expected output ###

test_that("debaseline functions correctly", {
    calculated_debased_amplification_curve_data <- debaseline(simulated_amplification_curve_data)
    
    expect_equal(calculated_debased_amplification_curve_data, simulated_debased_amplification_curve_data)
})

# Unit test for the calculate drdt function for melt curve data

### Create dataset of expected function outputs ###

plate_wells <- expand.grid(well_row = c("A", "B", "C", "D"), well_col = 1:4) %>%
    tidyr::unite(well, c(well_row, well_col), sep = "") %>%
    dplyr::pull(well)

simulated_melt_curve_data_diff <- 
    tibble(well = rep(plate_wells, times = 15),
           temperature = rep(1:15, each = 16),
           fluor_raw = rep(c(10, 15, 13, 20, 20, 50, 60, 70,
                             90, 130, 170, 280, 440, 670, 1020),
                           each = 16))

simulated_melt_curve_data_spline <- tibble(well = rep(plate_wells, times = 15),
                                         temperature = rep(1:15, each = 16),
                                         fluor_raw = rep(1:15, each = 16))

simulated_dRdT_melt_curve_data_diff <- simulated_melt_curve_data_diff %>%
    dplyr::mutate(dRdT = rep(c(-5, 2, -7, 0, -30, -10, -10, -20, -40,
                               -40, -110, -160, -230, -350, NA),
                             each = 16)) %>%
    dplyr::arrange(well, temperature)

simulated_dRdT_melt_curve_data_spline <- simulated_melt_curve_data_spline %>%
    dplyr::mutate(dRdT = rep(-1, each = 240)) %>%
    dplyr::arrange(well, temperature)

### Test functions give expected output ###

test_that("debaseline functions correctly", {
    calculated_dRdT_melt_curve_data_diff <- calculate_drdt_plate(simulated_melt_curve_data_diff, method = "diff")
    
    calculated_dRdT_melt_curve_data_spline <- calculate_drdt_plate(simulated_melt_curve_data_spline, method = "spline")
    
    expect_equal(calculated_dRdT_melt_curve_data_diff, simulated_dRdT_melt_curve_data_diff)
    
    expect_equal(calculated_dRdT_melt_curve_data_spline, simulated_dRdT_melt_curve_data_spline)
})
