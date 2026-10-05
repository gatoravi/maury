# Maury

Maury is an R package that measures how genetically related two sequencing
samples are, directly from their BAM files. Its main use is catching **sample
swaps and mix-ups** in next-generation sequencing (NGS) studies.

## Why

In large sequencing projects, samples get mislabeled. A tumor gets paired with
the wrong normal, two "unrelated" individuals turn out to be the same person, or
an RNA-seq library is attached to the wrong DNA sample. These errors silently
corrupt downstream results such as somatic mutation calls, eQTL analyses and
association studies.

Maury checks sample identity by comparing genotypes at a small panel of common
polymorphisms. Samples from the same individual share nearly all genotypes,
first-degree relatives share fewer, and unrelated individuals share fewer
still.

## How it works

1. **Pileup.** For each site in a VCF of common SNPs, Maury counts reference
   and alternate allele reads in both BAMs using
   [Rsamtools](https://bioconductor.org/packages/Rsamtools/), with mapping and
   base quality filters.
2. **Genotyping.** At sites with enough coverage, it picks the most likely
   genotype (hom-ref, het or hom-alt) by maximum likelihood under a binomial
   model with a small sequencing error rate.
3. **Concordance.** It reports the proportion of sites where the two samples
   have the same genotype.

## SNP panels

Maury ships with SNP panels in `inst/extdata/` based on the sample-tracking
panel from Pengelly *et al.*,
[*A SNP profiling panel for sample tracking in whole-exome sequencing studies*](http://www.ncbi.nlm.nih.gov/pubmed/24070238),
Genome Medicine (2013):

| File | Description |
| --- | --- |
| `optimized_grch37.vcf` | Recommended panel, GRCh37 coordinates |
| `optimized_grch38.vcf` | Recommended panel, GRCh38 coordinates |
| `all_grch37.vcf` | Larger panel, GRCh37 coordinates |

Use the panel that matches the reference genome your BAMs were aligned to.

## Install

Maury depends on Bioconductor packages, so install those first:

```r
install.packages(c("BiocManager", "devtools"))
BiocManager::install(c("Rsamtools", "VariantAnnotation", "GenomicRanges"))
devtools::install_github("gatoravi/maury")
```

## Usage

```r
library(maury)

vcf  <- system.file("extdata", "optimized_grch37.vcf", package = "maury")
maury(sample = "patient1_tumor_vs_normal",
      vcf  = vcf,
      bam1 = "patient1_tumor.bam",
      bam2 = "patient1_normal.bam",
      op_f = "patient1_concordance.tsv")
```

A runnable example with the bundled test data:

```r
library(maury)
bam <- system.file("extdata", "ex1.bam", package = "maury", mustWork = TRUE)
vcf <- system.file("extdata", "ex1.vcf", package = "maury", mustWork = TRUE)
maury("test_sample", vcf, bam, bam, min_rd = 0, min_mq = 0, min_bq = 0)
```

### Parameters

| Argument | Default | Meaning |
| --- | --- | --- |
| `sample` |  | Label for this comparison, written to the output |
| `vcf` |  | VCF of sites to genotype |
| `bam1`, `bam2` |  | Indexed BAM files to compare |
| `op_f` | `""` | Optional path for a tab-separated output file |
| `min_rd` | 20 | Minimum read depth to genotype a site |
| `seq_error` | 0.001 | Sequencing error rate used in the likelihood |
| `min_mq` | 20 | Minimum mapping quality |
| `min_bq` | 20 | Minimum base quality |

### Output

A data frame with one row per comparison: the sample label, the number of
matching genotypes, the number of sites genotyped in both samples, and the
proportion that match. A proportion close to 1 means the two BAMs very likely
come from the same individual.

For a sample analysis, see the [wiki](https://github.com/gatoravi/maury/wiki).

## License

MIT
