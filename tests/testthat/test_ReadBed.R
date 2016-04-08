test_that("ReadBed of non-existing file produces error", {
    expect_error(ReadBED("ABC"));
    
});

bedFile <- system.file("extdata",
                       "miCLIP_m6A_Linder2015_hg38.bed",
                       package = "RNAModR");
test_that("ReadBed of exisiting BED file produces GRanges", {
    expect_is(ReadBED(bedFile), "GRanges"); 
    expect_equal(length(ReadBED(bedFile)), 15167);
});
