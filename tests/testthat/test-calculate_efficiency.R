# Unit test for the calculate primer efficiency function

### Create dataset of expected function outputs ###

simulated_dilution_cq_dataset <- tibble(dilution = rep(c(1, 0.5, 0.25,
                                                         0.125, 1, 1),
                                                       times = 4),
                                        biol_rep = rep(c("A", "B",
                                                         "A", "B"),
                                                       each = 6),
                                        prep_type = rep(c("+RT", "+RT", "+RT",
                                                          "+RT", "-RT", "NT"),
                                                        times = 4),
                                        target_id = rep(c("Target_1", "Target_1",
                                                          "Target_2", "Target_2"),
                                                        each = 6),
                                        cq = rep(c(3, 6, 9,
                                                   12, 3, 3),
                                                 times = 4))

simulated_dilution_efficiency_dataset <- tibble(target_id = c("Target_1", "Target_2"),
                                                efficiency = c("log2(dilution)" = 3,
                                                               "log2(dilution)" = 3),
                                                efficiency.sd = 0,
                                                r.squared = 1)

### Test functions give expected output ###
test_that("functions for calculating efficiency.", {
    calculated_dilution_efficiency_dataset <- calculate_efficiency_bytargetid(simulated_dilution_cq_dataset)

    expect_equal(calculated_dilution_efficiency_dataset,
                 simulated_dilution_efficiency_dataset)
})
