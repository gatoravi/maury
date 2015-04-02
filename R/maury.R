#' Compare two samples' genotypes across a list of given positions
#'
#' @export
#' @param vcf VCF with list of sites to compare genotypes against.
#' @param bam1 BAM of first sample.
#' @param bam2 BAM of second sample.
#' @return p-value of samples being identical.
#' @examples
#' \donttest{
#' maury(vcf = "optimized1.vcf", bam1 = "sample1.bam", bam2 = "sample2.bam")
#' }
maury <- function(vcf, bam1, bam2, genome = "hg19") {
    vcf_data <- VariantAnnotation::readVcf(
        vcf, "hg19", collapsed = TRUE)
    pfiles <- Rsamtools::PileupFiles(c(bam1, bam2))
    bams <- names(pfiles)
    pileup <- get_pileup(pfiles, vcf_data)
    genotypes <- sapply(pileup, get_genotypes, bams)
    gt1 <- unlist(genotypes["gt1", ])
    gt2 <- unlist(genotypes["gt2", ])
    gt1 <- gt1[!is.na(gt1)]
    gt2 <- gt2[!is.na(gt2)]
    score <- get_score(gt1, gt2)
    score
}

#' calculate proportion of genotypes that match
#'
#' @param gt1 vector of genotypes
#' @param gt2 vector of genotypes
#' @return proportion of gts which are concordant
get_score <- function(gt1, gt2) {
    if(length(gt1) & length(gt2)) {
        #TODO:add assert to check lengths are equal
        return(sum(mapply(identical, gt1, gt2))/length(gt1))
    } else {
        print("No genotypes found")
        return(NA)
    }
}

#' Return best genotype using MLE
#'
#' @param ref_rd readcount for the reference allele
#' @param alt_rd readcount for the alternate allele
#' @param seq_error sequencing error-rate[0.001]
#' @return genotype 0=R/R, 1=R/A, 2=A/A
get_best_genotype <- function(ref_rd, alt_rd, seq_error = 0.001) {
    max_genotype <- which.max(c(
        ref_rd * log(1 - seq_error) + alt_rd * log(seq_error),
        ref_rd * log(0.5) + alt_rd * log(0.5),
        ref_rd * log(seq_error) + alt_rd * log(1 - seq_error))) - 1
}

#' callback for applyPileups
#'
#' @param x list with elements describing the current pile-up.
#' @return list with seqnames, positions and depths for alleles
calc_depths <- function(x) {
    depths <- apply(x[["seq"]], 2, function(y) {
        if(ncol(y)) {
            y <- y[c("A", "C", "G", "T", "N"), ,
                drop = FALSE]
        #no ReadDepth case.
        } else {
            y <- rep(0, 5)
        }
        names(y) <- c("A", "C", "G", "T", "N")
        y
    })
    list(seqnames = x[["seqnames"]], pos=x[["pos"]], depths = depths)
}

#' Compare genotypes of two samples
#'
#' @param apileup1 list with ref,alt alleles and read-depths
#' @param bams vector containing BAMs from both samples
#' @param min_rd minimum read-depth to calculate genotypes[20]
#' @return list of both genotypes
get_genotypes <- function(apileup1, bams, min_rd = 20) {
    ref <- apileup1$ref
    alt <- apileup1$alt
    depths <- apileup1$depths
    bam1 <- bams[1]
    bam2 <- bams[2]
    refc1 <- depths[ref, bam1]
    altc1 <- depths[alt, bam1]
    refc2 <- depths[ref, bam2]
    altc2 <- depths[alt, bam2]
    #ignore no coverage sites.
    if((refc1 + altc1) < min_rd |
        (refc2 + altc2) < min_rd) {
        return(list(gt1 = NA, gt2 = NA))
    } else {
        gt1 = get_best_genotype(refc1, altc1)
        gt2 = get_best_genotype(refc2, altc2)
        return(list(gt1 = gt1, gt2 = gt2))
    }
}

#' Append ref, alt alleles to Pileup list
#'
#' @param apileup1 list of lists, each sub-list is for a different chr-pos
#' @param refs vector of reference-alleles at different chr-pos
#' @param alts vector of alternate-alleles at different chr-pos
#' @return appended list
append_ref_alt <- function(apileup, refs, alts) {
    for(i in 1:length(apileup)) {
        apileup[[i]]$ref = refs[i]
        apileup[[i]]$alt = alts[i]
    }
    apileup
}

#' get pileup at supplied list of sites.
#'
#' @param pfiles reference to BAM files using Rsamtools::PileupFiles
#' @param vcf_data VCF object returned by VariantAnnotation::readVcf
#' @return appended list
get_pileup <- function(pfiles, vcf_data) {
    evcf <- VariantAnnotation::expand(vcf_data)
    ranges <- GenomicRanges::rowData(evcf)
    ap_param <- Rsamtools::ApplyPileupsParam(
        minMapQuality = 20,
        minBaseQuality = 20,
        what = "seq",
        which = ranges)
    apileup <- Rsamtools::applyPileups(
        pfiles, calc_depths, param = ap_param)
    refs <- as.character(VariantAnnotation::ref(evcf))
    alts <- as.character(VariantAnnotation::alt(evcf))
    stopifnot(length(apileup) == length(refs),
              length(apileup) == length(alts))
    apileup <- append_ref_alt(apileup, refs, alts)
}
