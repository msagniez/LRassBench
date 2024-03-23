#Prepare table with bash
#FP = expected abundance = 0
#GFF
#echo -e "Sequin\tMixA\tType\tTranscript\tClasscode\tTPM\tnbReads" > headerGFF.tsv
#for file in *.sf ; do join -1 1 -2 5 <( sort -k1,1 $file ) <( sort -k5,5 ${file%*_salmon.sf}.tmap ) | sed 's/ /\t/g' | awk '$8=="=" || $8=="c" {print "TP\t"$1"\t"$7"\t"$8"\t"$4"\t"$5}' - > GFF-${file%*_salmon.sf}_quantif-TP.tsv; join -t $'\t' -1 1 -2 3 <( sort -k1,1 rnasequin_isoforms_2.5.tsv ) <( sort -k3,3 GFF-${file%*_salmon.sf}_quantif-TP.tsv ) | cut -d$'\t' -f5,6,1,7,8,9,3 > GFF-${file%*_salmon.sf}_quantif-TP2.tsv; join -1 1 -2 5 <( sort -k1,1 $file ) <( sort -k5,5 ${file%*_salmon.sf}.tmap ) | sed 's/ /\t/g' | awk '$8!="=" && $8!="c" && $8!="u" {print $7"\t0\tFP\t"$1"\t"$8"\t"$4"\t"$5"\t"}' - > GFF-${file%*_salmon.sf}_quantif-FP.tsv; cat headerGFF.tsv GFF-${file%*_salmon.sf}_quantif-TP2.tsv GFF-${file%*_salmon.sf}_quantif-FP.tsv > GFF-${file%*_salmon.sf}_quantif.tsv; rm GFF-${file%*_salmon.sf}_quantif-FP.tsv GFF-${file%*_salmon.sf}_quantif-#TP.tsv GFF-${file%*_salmon.sf}_quantif-TP2.tsv; done
#rm headerGFF.tsv


#SQ3
#echo -e "Sequin\tMixA\tType\tTranscript\tClasscode1\tClasscode2\tTPM\tnbReads" > headerSQ3.tsv
#for file in *.sf ; do join -1 1 -2 1 <( sort -k1,1 $file ) <( sort -k1,1 ${file%*_salmon.sf}_classification.txt | sed 's/%7C/|/g' - ) | sed 's/ /\t/g' | awk '$10=="full-splice_match" || ( $10=="incomplete-splice_match" && $19!="intron_retention" ) {print "TP\t"$1"\t"$12"\t"$10"\t"$19"\t"$4"\t"$5}' - > SQ3-${file%*_salmon.sf}_quantif-TP.tsv; join -t $'\t' -1 1 -2 3 <( sort -k1,1 rnasequin_isoforms_2.5.tsv ) <( sort -k3,3 SQ3-${file%*_salmon.sf}_quantif-TP.tsv ) | cut -d$'\t' -f5,6,1,7,8,9,10,3 > SQ3-${file%*_salmon.sf}_quantif-TP2.tsv; join -1 1 -2 1 <( sort -k1,1 $file ) <( sort -k1,1 ${file%*_salmon.sf}_classification.txt | sed 's/%7C/|/g' - ) | sed 's/ /\t/g' | awk '$10!="full-splice_match" && ( $10!="incomplete-splice_match" || $19=="intron_retention" ) && $10!="intergenic" {print $12"\t0\tFP\t"$1"\t"$10"\t"$19"\t"$4"\t"$5}' - > SQ3-${file%*_salmon.sf}_quantif-FP.tsv; cat headerSQ3.tsv SQ3-${file%*_salmon.sf}_quantif-TP2.tsv SQ3-${file%*_salmon.sf}_quantif-FP.tsv > SQ3-${file%*_salmon.sf}_quantif.tsv; rm SQ3-${file%*_salmon.sf}_quantif-FP.tsv SQ3-${file%*_salmon.sf}_quantif-TP.tsv SQ3-${file%*_salmon.sf}_quantif-TP2.tsv; done
#rm headerSQ3.tsv

#Plot figure with R
#RStudio
library(ggplot2)
library(ggpubr)
library(cowplot)

setwd("D:/Assembly-Benchmarking/Fig5-Quantification/Upt2")

#All sequins abundances sum = 15196.5066 for TPM comversion

GFF_LSK114_Bambu <- read.csv("GFF-LSK114_Bambu_quantif.tsv", header=T, sep="\t")
GFF_LSK114_Bambu$nbreadslog <- log10(GFF_LSK114_Bambu$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_LSK114_Bambu$MixAlog <- log10(((GFF_LSK114_Bambu$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_LSK114_Bambu_x <- subset(subset(GFF_LSK114_Bambu, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_LSK114_Bambu_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_LSK114_Bambu_plot <- ggplot(GFF_LSK114_Bambu,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_Bambu_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_Bambu_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_Bambu_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

GFF_LSK114_isoQuant <- read.csv("GFF-LSK114_isoQuant_quantif.tsv", header=T, sep="\t")
GFF_LSK114_isoQuant$nbreadslog <- log10(GFF_LSK114_isoQuant$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_LSK114_isoQuant$MixAlog <- log10(((GFF_LSK114_isoQuant$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_LSK114_isoQuant_x <- subset(subset(GFF_LSK114_isoQuant, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_LSK114_isoQuant_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_LSK114_isoQuant_plot <- ggplot(GFF_LSK114_isoQuant,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_isoQuant_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_isoQuant_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_isoQuant_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

GFF_LSK114_Flair <- read.csv("GFF-LSK114_FLAIR-noRef_quantif.tsv", header=T, sep="\t")
GFF_LSK114_Flair$nbreadslog <- log10(GFF_LSK114_Flair$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_LSK114_Flair$MixAlog <- log10(((GFF_LSK114_Flair$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_LSK114_Flair_x <- subset(subset(GFF_LSK114_Flair, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_LSK114_Flair_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_LSK114_Flair_plot <- ggplot(GFF_LSK114_Flair,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_Flair_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_Flair_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_Flair_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

GFF_LSK114_Stringtie2 <- read.csv("GFF-LSK114_Stringtie2-noRef_quantif.tsv", header=T, sep="\t")
GFF_LSK114_Stringtie2$nbreadslog <- log10(GFF_LSK114_Stringtie2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_LSK114_Stringtie2$MixAlog <- log10(((GFF_LSK114_Stringtie2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_LSK114_Stringtie2_x <- subset(subset(GFF_LSK114_Stringtie2, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_LSK114_Stringtie2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_LSK114_Stringtie2_plot <- ggplot(GFF_LSK114_Stringtie2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_Stringtie2_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_Stringtie2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_Stringtie2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

SQ3_LSK114_Bambu <- read.csv("SQ3-LSK114_Bambu_quantif.tsv", header=T, sep="\t")
SQ3_LSK114_Bambu$nbreadslog <- log10(SQ3_LSK114_Bambu$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_LSK114_Bambu$MixAlog <- log10(((SQ3_LSK114_Bambu$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_LSK114_Bambu_x <- subset(subset(SQ3_LSK114_Bambu, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_LSK114_Bambu_x[,9:10]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_LSK114_Bambu_plot <- ggplot(SQ3_LSK114_Bambu,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_Bambu_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_Bambu_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_Bambu_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

SQ3_LSK114_isoQuant <- read.csv("SQ3-LSK114_isoQuant_quantif.tsv", header=T, sep="\t")
SQ3_LSK114_isoQuant$nbreadslog <- log10(SQ3_LSK114_isoQuant$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_LSK114_isoQuant$MixAlog <- log10(((SQ3_LSK114_isoQuant$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_LSK114_isoQuant_x <- subset(subset(SQ3_LSK114_isoQuant, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_LSK114_isoQuant_x[,9:10]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_LSK114_isoQuant_plot <- ggplot(SQ3_LSK114_isoQuant,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_isoQuant_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_isoQuant_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_isoQuant_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

SQ3_LSK114_Flair <- read.csv("SQ3-LSK114_FLAIR-noRef_quantif.tsv", header=T, sep="\t")
SQ3_LSK114_Flair$nbreadslog <- log10(SQ3_LSK114_Flair$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_LSK114_Flair$MixAlog <- log10(((SQ3_LSK114_Flair$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_LSK114_Flair_x <- subset(subset(SQ3_LSK114_Flair, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_LSK114_Flair_x[,9:10]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_LSK114_Flair_plot <- ggplot(SQ3_LSK114_Flair,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_Flair_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_Flair_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_Flair_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

SQ3_LSK114_Stringtie2 <- read.csv("SQ3-LSK114_Stringtie2-noRef_quantif.tsv", header=T, sep="\t")
SQ3_LSK114_Stringtie2$nbreadslog <- log10(SQ3_LSK114_Stringtie2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_LSK114_Stringtie2$MixAlog <- log10(((SQ3_LSK114_Stringtie2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_LSK114_Stringtie2_x <- subset(subset(SQ3_LSK114_Stringtie2, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_LSK114_Stringtie2_x[,9:10]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_LSK114_Stringtie2_plot <- ggplot(SQ3_LSK114_Stringtie2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_Stringtie2_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_Stringtie2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_Stringtie2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

GFF_RNA004_Bambu <- read.csv("GFF-RNA004_Bambu_quantif.tsv", header=T, sep="\t")
GFF_RNA004_Bambu$nbreadslog <- log10(GFF_RNA004_Bambu$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_RNA004_Bambu$MixAlog <- log10(((GFF_RNA004_Bambu$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_RNA004_Bambu_x <- subset(subset(GFF_RNA004_Bambu, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_RNA004_Bambu_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_RNA004_Bambu_plot <- ggplot(GFF_RNA004_Bambu,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_Bambu_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_Bambu_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_Bambu_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

GFF_RNA004_isoQuant <- read.csv("GFF-RNA004_isoQuant_quantif.tsv", header=T, sep="\t")
GFF_RNA004_isoQuant$nbreadslog <- log10(GFF_RNA004_isoQuant$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_RNA004_isoQuant$MixAlog <- log10(((GFF_RNA004_isoQuant$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_RNA004_isoQuant_x <- subset(subset(GFF_RNA004_isoQuant, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_RNA004_isoQuant_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_RNA004_isoQuant_plot <- ggplot(GFF_RNA004_isoQuant,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_isoQuant_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_isoQuant_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_isoQuant_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

GFF_RNA004_Flair <- read.csv("GFF-RNA004_FLAIR-noRef_quantif.tsv", header=T, sep="\t")
GFF_RNA004_Flair$nbreadslog <- log10(GFF_RNA004_Flair$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_RNA004_Flair$MixAlog <- log10(((GFF_RNA004_Flair$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_RNA004_Flair_x <- subset(subset(GFF_RNA004_Flair, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_RNA004_Flair_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_RNA004_Flair_plot <- ggplot(GFF_RNA004_Flair,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_Flair_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_Flair_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_Flair_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

GFF_RNA004_Stringtie2 <- read.csv("GFF-RNA004_Stringtie2-noRef_quantif.tsv", header=T, sep="\t")
GFF_RNA004_Stringtie2$nbreadslog <- log10(GFF_RNA004_Stringtie2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_RNA004_Stringtie2$MixAlog <- log10(((GFF_RNA004_Stringtie2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_RNA004_Stringtie2_x <- subset(subset(GFF_RNA004_Stringtie2, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_RNA004_Stringtie2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_RNA004_Stringtie2_plot <- ggplot(GFF_RNA004_Stringtie2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_Stringtie2_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_Stringtie2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_Stringtie2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

SQ3_RNA004_Bambu <- read.csv("SQ3-RNA004_Bambu_quantif.tsv", header=T, sep="\t")
SQ3_RNA004_Bambu$nbreadslog <- log10(SQ3_RNA004_Bambu$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_RNA004_Bambu$MixAlog <- log10(((SQ3_RNA004_Bambu$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_RNA004_Bambu_x <- subset(subset(SQ3_RNA004_Bambu, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_RNA004_Bambu_x[,9:10]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_RNA004_Bambu_plot <- ggplot(SQ3_RNA004_Bambu,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_Bambu_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_Bambu_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_Bambu_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

SQ3_RNA004_isoQuant <- read.csv("SQ3-RNA004_isoQuant_quantif.tsv", header=T, sep="\t")
SQ3_RNA004_isoQuant$nbreadslog <- log10(SQ3_RNA004_isoQuant$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_RNA004_isoQuant$MixAlog <- log10(((SQ3_RNA004_isoQuant$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_RNA004_isoQuant_x <- subset(subset(SQ3_RNA004_isoQuant, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_RNA004_isoQuant_x[,9:10]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_RNA004_isoQuant_plot <- ggplot(SQ3_RNA004_isoQuant,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_isoQuant_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_isoQuant_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_isoQuant_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

SQ3_RNA004_Flair <- read.csv("SQ3-RNA004_FLAIR-noRef_quantif.tsv", header=T, sep="\t")
SQ3_RNA004_Flair$nbreadslog <- log10(SQ3_RNA004_Flair$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_RNA004_Flair$MixAlog <- log10(((SQ3_RNA004_Flair$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_RNA004_Flair_x <- subset(subset(SQ3_RNA004_Flair, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_RNA004_Flair_x[,9:10]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_RNA004_Flair_plot <- ggplot(SQ3_RNA004_Flair,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_Flair_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_Flair_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_Flair_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))

SQ3_RNA004_Stringtie2 <- read.csv("SQ3-RNA004_Stringtie2-noRef_quantif.tsv", header=T, sep="\t")
SQ3_RNA004_Stringtie2$nbreadslog <- log10(SQ3_RNA004_Stringtie2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_RNA004_Stringtie2$MixAlog <- log10(((SQ3_RNA004_Stringtie2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_RNA004_Stringtie2_x <- subset(subset(SQ3_RNA004_Stringtie2, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_RNA004_Stringtie2_x[,9:10]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_RNA004_Stringtie2_plot <- ggplot(SQ3_RNA004_Stringtie2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_Stringtie2_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_Stringtie2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_Stringtie2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred"))



plot_grid(GFF_LSK114_Bambu_plot,GFF_LSK114_isoQuant_plot,GFF_LSK114_Stringtie2_plot,GFF_LSK114_Flair_plot,SQ3_LSK114_Bambu_plot,SQ3_LSK114_isoQuant_plot,SQ3_LSK114_Stringtie2_plot,SQ3_LSK114_Flair_plot,GFF_RNA004_Bambu_plot,GFF_RNA004_isoQuant_plot,GFF_RNA004_Stringtie2_plot,GFF_RNA004_Flair_plot,SQ3_RNA004_Bambu_plot,SQ3_RNA004_isoQuant_plot,SQ3_RNA004_Stringtie2_plot,SQ3_RNA004_Flair_plot, nrow=4)


#PDF 10x10 inches
