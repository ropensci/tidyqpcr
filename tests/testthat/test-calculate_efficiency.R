test_that("calculate_efficiency outputs warning to use calculate_efficiency_plate if more than one target_id found", {
  warning_tibble <- tibble(target_id = rep(c("PGK1_ORF", "RPPS3_ORF"), each = 3),
                           cq = c(8,10,12,7,8,9),
                           biol_rep = 1,
                           dilution = c(1,0.5,0.25,1,0.5,0.25))
  
  expect_warning(calculate_efficiency(warning_tibble), "multiple target_ids, did you mean calculate_efficiency_plate?")
})

test_that("calculate_efficiency calculates correct values for a given target_id", {
    test_dilution_tibble <- tibble(target_id = rep("PGK1_ORF", times = 9),
                             cq = c(8,7.7,7.5,10,10.1,10.2,11.5,12.5,12),
                             biol_rep = rep(c(1,2,3), times = 3),
                             dilution = rep(c(1,0.5,0.25), each = 3))
    
    test_efficiency_tibble <- tibble(efficiency = c(`log2(dilution)` = 2.133333), # Output of lm$coefficient has underlying structure
               efficiency.sd = 0.1442306,
               r.squared = 0.9733133)
    
    expect_equal(calculate_efficiency(test_dilution_tibble), test_efficiency_tibble, tolerance = 0.002)
})

test_that("calculate_efficiency_bytargetid calculates correct values for a given target_id", {
    test_multi_dilution_tibble <- tibble(target_id = rep(c("PGK1_ORF", "RPS3_ORF"), each = 9),
                                   cq = c(8,7.7,7.5,10,10.1,10.2,11.5,12.5,12,6,6.7,6.9,8.3,9.3,7.6,10.9,12.5,12),
                                   biol_rep = rep(c(1,2,3), times = 6),
                                   dilution = rep(rep(c(1,0.5,0.25), each = 3), times =2))
    
    test_multi_efficiency_tibble <- tibble(target_id = c("PGK1_ORF", "RPS3_ORF"),
                                     efficiency = c(`log2(dilution)` = 2.133333, `log2(dilution)` = 2.6333), # Output of lm$coefficient has underlying named list structure
                                     efficiency.sd =c(0.1442306, 0.3391),
                                     r.squared = c(0.9733133,0.91))
    
    expect_equal(calculate_efficiency_bytargetid(test_multi_dilution_tibble, use_prep_types=NA), test_multi_efficiency_tibble, tolerance = 0.002)
})
