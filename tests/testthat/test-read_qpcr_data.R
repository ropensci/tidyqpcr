# Unit tests for reading qpcr data from files

### Test functions run and load data frames ###

cq_data_names <- c("include", "color", "well", "sample_info", 
                   "cq", "concentration", "standard", "status") 

fluor_data_names <- c("well", "sample_info", "program_no", "segment_no",
                      "cycle", "time", "temperature", "fluor_raw")

test_that("functions for loading data from Roche Lightcycler software", {
    cq_data_96well <- 
        read_lightcycler_1colour_cq(
            system.file("extdata/Stuart_dAgr_glyS_spoVG_5S_individualWells_Cq.txt.gz",
                        package = "tidyqpcr"))
    
    expect_s3_class(cq_data_96well, "data.frame")
    
    expect_named(cq_data_96well, cq_data_names)
    
    fluor_data_384well <- 
        read_lightcycler_1colour_raw(
            system.file("extdata/Edward_qPCR_Nrd1_calibration_2019-02-02.txt.gz",
                        package = "tidyqpcr")
        )
    
    expect_s3_class(fluor_data_384well, "data.frame")
    
    expect_named(fluor_data_384well, fluor_data_names)
})
