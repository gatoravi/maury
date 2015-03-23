#' Compare two samples' genotypes across a list of given positions
#'
#' @export
#' @param sites VCF with list of sites to compare genotypes against.
#' @examples
#' \donttest{
#' maury(vcf = "optimized1.vcf", bam1 = "sample1.bam", bam2 = "sample2.bam")
#' }
maury <- function(vcf_f, bam1, bam2, genome = "hg19") {
    vcf_data <- VariantAnnotation::readVcf(
        vcf_f, "hg19", collapsed = TRUE)
    pfiles <- Rsamtools::PileupFiles(c(bam1, bam2))
    bams <- names(pfiles)
    pileup <- get_pileupfreq(pfiles, vcf_data)
    pvalues <- sapply(pileup, apply_binomial, bams)
    apply_fishers_method(pvalues[!is.na(pvalues)])
}

apply_fishers_method <- function(pvalues) {
    1 - pchisq(-2 * sum(log(pvalues)), df = 2 * length(pvalues))
}

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
    #no ReadDepth case.
    if(length(x[["pos"]]) == 0) {
        return(list(seqnames = "empty", pos = "empty", depths = depths))
    }
    list(seqnames = x[["seqnames"]], pos=x[["pos"]], depths = depths)
}

apply_binomial <- function(apileup1, bams) {
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
    if((refc1 + altc1) == 0 |
        (refc2 + altc2) == 0) {
            return(NA)
    } else {
        ft = fisher.test(matrix(c(refc1, altc1,
            refc2, altc2), nrow=2))
        return(ft$p.value)
    }
}

append_ref_alt <- function(apileup, refs, alts) {
    for(i in 1:length(apileup)) {
        apileup[[i]]$ref = refs[i]
        apileup[[i]]$alt = alts[i]
    }
    apileup
}

get_pileupfreq <- function(pfiles, vcf_data) {
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
