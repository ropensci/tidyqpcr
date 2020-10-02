test_that("Check debaseline subtracts median of first 10 cycles from all cycles for each well", {
    
  raw_fluor_tibble <- tibble(well =rep(c("A1","B2","C5"), each = 20), 
             program_no = 2, 
             cycle = rep(seq(1,20), times = 3), 
             fluor_raw = 2^rep(seq(0.5,10,0.5), times = 3))
  
  debased_fluor_tibble <- tibble(well =rep(c("A1","B2","C5"), each = 20), 
                                 program_no = 2, 
                                 cycle = rep(seq(1,20), times = 3), 
                                 fluor_raw = 2^rep(seq(0.5,10,0.5),times = 3),
                                 fluor_base = median(2^rep(seq(0.5,5,0.5))),
                                 fluor_signal = fluor_raw - fluor_base)
  
  expect_equal(debaseline(raw_fluor_tibble),debased_fluor_tibble)
})

test_that("Check calculate_dydx_1 only accepts two numerical arguments of equal length", {
    expect_error(calculate_dydx_1("x",1), "x is not a numeric or integer vector")
    expect_error(calculate_dydx_1(1,"y"), "y is not a numeric or integer vector")
    expect_error(calculate_dydx_1(c(1,2),1))
})

test_that("Check calculate_dydx_1 uses correct method or prints an error", {
    expect_equal(calculate_dydx_1(c(1,5,8,10),c(4,12,11,9), method = "diff"), c(-2,1/3,1,NA))
    
    
    expect_equal(calculate_dydx_1(c(1,5,8,10),c(4,12,11,9), method = "spline"), 
                 -1 * stats::predict(object = stats::smooth.spline(c(1,5,8,10),c(4,12,11,9)), x = c(1,5,8,10), deriv = 1)$y)
    
    expect_error(calculate_dydx_1(c(1,5,8,10),c(4,12,11,9), method = "line"))
})

calculate_drdt_plate("Check calculate_drdt_plate runs across all wells", {
    
    plate_fluor_data <- tibble(well =rep(c("A1","B2","C5"), each = 20), 
           program_no = 2, 
           cycle = rep(seq(1,20), times = 3), 
           fluor_raw = 2^rep(seq(0.5,10,0.5), times = 3),
           temperature = rep(seq(21,40), times = 3))
    
    plate_dRdT_data <- tibble(well =rep(c("A1","B2","C5"), each = 20), 
                               program_no = 2, 
                               cycle = rep(seq(1,20), times = 3), 
                               fluor_raw = 2^rep(seq(0.5,10,0.5), times = 3),
                               temperature = rep(seq(21,40), times = 3),
                               dRdT = rep(-1 * stats::predict(object = stats::smooth.spline(seq(21,40),2^seq(0.5,10,0.5)), x = seq(21,40), deriv = 1)$y, times = 3))
    
    expect_equal(calculate_drdt_plate(plate_fluor_data), plate_dRdT_data)
})

calculate_drdt_plate("Check getdRdTall calls calculate_drdt_plate correctly", {
    
    plate_fluor_data <- tibble(well =rep(c("A1","B2","C5"), each = 20), 
                               program_no = 2, 
                               cycle = rep(seq(1,20), times = 3), 
                               fluor_raw = 2^rep(seq(0.5,10,0.5), times = 3),
                               temperature = rep(seq(21,40), times = 3))
    
    plate_dRdT_data <- tibble(well =rep(c("A1","B2","C5"), each = 20), 
                              program_no = 2, 
                              cycle = rep(seq(1,20), times = 3), 
                              fluor_raw = 2^rep(seq(0.5,10,0.5), times = 3),
                              temperature = rep(seq(21,40), times = 3),
                              dRdT = rep(-1 * stats::predict(object = stats::smooth.spline(seq(21,40),2^seq(0.5,10,0.5)), x = seq(21,40), deriv = 1)$y, times = 3))
    
    expect_equal(getdRdTall(plate_fluor_data),plate_dRdT_data)
})

