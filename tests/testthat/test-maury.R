context("maury")

test_that("check no read-count sites", {
    f1 <- system.file("extdata", 'ex1.bam',
        package = "maury", mustWork = TRUE)
    vcf <- system.file("extdata", 'ex1.vcf',
        package = "maury", mustWork = TRUE)
    expect_equal(maury(vcf, f1, f1), c(1))
})

