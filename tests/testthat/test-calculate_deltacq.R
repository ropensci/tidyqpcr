library(tidyverse)

test_that("calculate_normvalue applies norm function to correct groups", {
    test_cq_tibble <- tibble(target_id = rep(c("PGK1_ORF", "RPS3_ORF", "URA3_ORF"), each =6),
                                  sample_id = rep(rep(c("wildtype", "delta_SSD1"), each=3), times = 3),
                                  tech_rep = rep(c(1,2,3), times= 6),
                               cq = c(8.4,8.7,9,11.3,10.7,11.2,9.8,9.5,9.3,12.1,12.2,12.6,10.3,10.4,10.2,10.5,10.1,10.6))
    
    test_norm_tibble <- test_cq_tibble %>%
        mutate(value_to_norm_by = rep(10.35,times=18))
    expect_equal(calculate_normvalue(test_cq_tibble, "URA3_ORF", value_name = "cq", id_name = "target_id"), test_norm_tibble)
})

test_that("calculate_deltacq_bysampleid groups genes correctly for calculation", {
    test_cq_tibble <- tibble(target_id = rep(c("PGK1_ORF", "RPS3_ORF", "URA3_ORF"), each =6),
                             sample_id = rep(rep(c("wildtype", "delta_SSD1"), each=3), times = 3),
                             tech_rep = rep(c(1,2,3), times= 6),
                             cq = c(8.4,8.7,9,11.3,10.7,11.2,9.8,9.5,9.3,12.1,12.2,12.6,10.3,10.4,10.2,10.5,10.1,10.6))
    
    test_deltacq_tibble <- test_cq_tibble %>%
        mutate(ref_cq = rep(rep(c(10.3,10.5), each=3),times=3),
               delta_cq = c(-1.9,-1.6,-1.3,0.8,0.2,0.7,-0.5,-0.8,-1,1.6,1.7,2.1,0,0.1,-0.1,0,-0.4,0.1),
               rel_abund = 2^-delta_cq) %>% 
        arrange(sample_id)
        
    expect_equal(calculate_deltacq_bysampleid(test_cq_tibble, "URA3_ORF"), test_deltacq_tibble)
})

test_that("calculate_deltadeltacq_bytargetid groups genes correctly for calculation", {
    
    test_deltacq_tibble <- tibble(target_id = rep(c("PGK1_ORF", "RPS3_ORF", "URA3_ORF"), each =4),
                                  sample_id = rep(c("wildtype", "delta_SSD1", "delta_HSP26", "delta_HSP104"), times=3),
                                  delta_cq = c(6,1,3.7,5.5,4,3.9,3.8,4.5,2,1.9,0.3,2))
    
    test_deltadeltacq_tibble <- test_deltacq_tibble %>%
        mutate(ref_delta_cq = c(6,6,6,6,4,4,4,4,2,2,2,2), 
               deltadelta_cq = c(0,-5,-2.3,-0.5,0,-0.1,-0.2,0.5,0,-0.1,-1.7,0),
               fold_change = 2^-deltadelta_cq)
    
    expect_equal(calculate_deltadeltacq_bytargetid(test_deltacq_tibble, "wildtype"), test_deltadeltacq_tibble)
})
