#Maury

##Introduction
Maury is an R package that can be used to quantify relatedness of
multiple samples from Next Generation Sequencing data(NGS).

Maury works by looking at genotypes at sites of
common polymorphisms. Maury outputs the proportion of
genotypes that are concordant between a pair of samples. It uses
Maximum Likelihood Estimation(MLE) to infer the genotype of a sample at a particular
locus by looking at the read depths of the reference and alternate alleles.

A sample VCF with a recommended list of common polymorphisms
using information from Pengelly et al. (A SNP profiling panel for sample tracking in whole-exome sequencing studies.
http://www.ncbi.nlm.nih.gov/pubmed/24070238) is distributed with this package under `inst/extdata/`.

Maury can be used in a multitude
of applications, for example detecting tumor/normal
sample-swaps.


##Install
```r
devtools::install_github("gatoravi/maury")
```

##Example Usage
```r
library(maury)
bam1 <- system.file("extdata", 'ex1.bam',
    package = "maury", mustWork = TRUE)
bam2 <- system.file("extdata", 'ex1.bam',
    package = "maury", mustWork = TRUE)
vcf <- system.file("extdata", 'ex1.vcf',
    package = "maury", mustWork = TRUE)
maury("test_sample", vcf, bam1, bam2,
    min_rd = 0, min_mq = 0, min_bq = 0)
```

##Example output
For an example analysis using maury, please refer to
https://github.com/gatoravi/maury/wiki

