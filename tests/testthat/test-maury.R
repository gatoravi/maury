context("maury")

test_that("check no genotypes case", {
    bam <- system.file("extdata", 'ex1.bam',
        package = "maury", mustWork = TRUE)
    vcf <- system.file("extdata", 'ex1.vcf',
        package = "maury", mustWork = TRUE)
    expect_equal(maury("test_sample", vcf, bam, bam,
                        min_rd = 0, min_mq = 0, min_bq = 0),
                 list(sample = "test_sample", matches = 1,
                    total = 1, score = 1))
})

