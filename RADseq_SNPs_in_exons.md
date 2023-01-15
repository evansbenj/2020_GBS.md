This script will identify RADseq SNPs in exons:
```R

################################################
#
#
#  Deterine which RADseq SNPs are in exons
#
#
################################################


#BiocManager::install("genomation")
#library(genomation)
library("GenomicRanges")
library(genomation)
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/plink_somemappedtoXLv9.2_othersmappedtoAustinXB")
# read in the data from plink output
dat <-read.table("pygmaeus_filtered_removed_allchrs.vcf.gz_plink_noNAs.assoc", header=TRUE)
# make a -logP column
dat$minuslogP <- -log(dat$P)
# reverse sort the df by -logP
head(dat[order(-dat$minuslogP),],n=20)


gff <- gffToGRanges("XENLA_9.2_Xenbase.gtf", filter = NULL, zero.based = FALSE, ensembl = FALSE)


# get the RADseq SNPs that have a strong bias
subset<-dat[(dat$CHR == 'chr8L')&(dat$minuslogP >= 20),];subset
dim(subset)
# subset gtf and keep only the ranges exons in chr8L that are exons
chr8L_exons <- gff %>% subset(seqnames == 'chr8L') %>% subset(type == 'exon')
# make a GenomicRanges object from the RADseq SNPs
SL_sites <- makeGRangesFromDataFrame(subset,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="BP",
                         end.field=c("BP"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)

o = findOverlaps(SL_sites,chr8L_exons);o
SL_sites = split(SL_sites[queryHits(o)], 1:length(o)) # You can't mendoapply on a GRanges object
chr8L_exons = split(chr8L_exons[subjectHits(o)], 1:length(o))
foo = function(x, y) {
    rv = x
    start(rv) = max(start(x), start(y))
    end(rv) = min(end(x), end(y))
    return(rv)
}
unlist(mendoapply(foo, SL_sites, y=chr8L_exons))
```
