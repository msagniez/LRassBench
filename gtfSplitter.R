library(GenomicRanges)
library(rtracklayer)
gtf.gr <- import("/Users/anshul/Documents/PhD_UDEM/work-year2/Benchmarking/rnasequin_annotation_2.4.gtf", format="gtf") ##OG gtf
gtf.gr <- sortSeqlevels(gtf.gr)
gtf.gr@ranges ##1411 rows; starts 144926; ends 10233186
### length of fasta: 10567884 / 3 = 3522628
q1=GRanges(seqnames="chrIS",
          ranges=IRanges(start = 1, end = 3522628)) 
q2=GRanges(seqnames="chrIS",
           ranges=IRanges(start = 3522628, end = 7045256)) 
q3=GRanges(seqnames="chrIS",
           ranges=IRanges(start = 7045256, end = 10567884)) 

gtf.gr1 <- subsetByOverlaps(gtf.gr, q1) ##286 ##14 genes
gtf.gr2 <- subsetByOverlaps(gtf.gr, q2) ##664 ##35 genes
gtf.gr3 <- subsetByOverlaps(gtf.gr, q3) ##461 ##27 genes
  
# export(gtf.gr1, "rnasequin_annotation_2.4_g1.gtf")
# export(gtf.gr2, "rnasequin_annotation_2.4_g2.gtf")
# export(gtf.gr3, "rnasequin_annotation_2.4_g3.gtf")

##get gene IDs and subset by that

data1_genes <- unique(gtf.gr1$gene_id)[1:round(length(unique(gtf.gr1$gene_id))/2)]
data2_genes <- unique(gtf.gr2$gene_id)[1:round(length(unique(gtf.gr2$gene_id))/2)]
data3_genes <- unique(gtf.gr3$gene_id)[1:round(length(unique(gtf.gr3$gene_id))/2)]

gtf.gr1_half <- gtf.gr1[which(gtf.gr1$gene_id %in% data1_genes),] ##135 ##7 genes
gtf.gr2_half <- gtf.gr2[which(gtf.gr2$gene_id %in% data2_genes),] ##367 ##18 genes
gtf.gr3_half <- gtf.gr3[which(gtf.gr3$gene_id %in% data3_genes),] ##240 ##14 genes

# export(gtf.gr1_half, "rnasequin_annotation_2.4_g1_half.gtf")
# export(gtf.gr2_half, "rnasequin_annotation_2.4_g2_half.gtf")
# export(gtf.gr3_half, "rnasequin_annotation_2.4_g3_half.gtf")

full_half_null<-c(gtf.gr1,gtf.gr2_half)
full_null_half<-c(gtf.gr1,gtf.gr3_half)
half_full_null<-c(gtf.gr1_half,gtf.gr2)
half_null_full<-c(gtf.gr1_half,gtf.gr3)
null_full_half<-c(gtf.gr2,gtf.gr3_half)
null_half_full<-c(gtf.gr2_half,gtf.gr3)

export(full_half_null, "rnasequin_annotation_2.4_full_half_null.gtf")
export(full_null_half, "rnasequin_annotation_2.4_full_null_half.gtf")
export(half_full_null, "rnasequin_annotation_2.4_half_full_null.gtf")
export(half_null_full, "rnasequin_annotation_2.4_half_null_full.gtf")
export(null_full_half, "rnasequin_annotation_2.4_null_full_half.gtf")
export(null_half_full, "rnasequin_annotation_2.4_null_half_full.gtf")



