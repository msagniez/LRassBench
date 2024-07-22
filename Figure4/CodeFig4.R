#Prepare table with bash
#FP = expected abundance = 0
#GFF
#echo -e "Sequin\tMixA\tType\tTranscript\tClasscode\tTPM\tnbReads" > headerGFF.tsv
#for file in *.txt ; do join -1 1 -2 5 <( sort -k1,1 $file ) <( sort -k5,5 ../${file%*.txt}.tmap ) | sed 's/ /\t/g' | awk '$8=="=" || $8=="c" {print "TP\t"$1"\t"$7"\t"$8"\t"$4"\t"$5}' - > ../GFF-${file%*.txt}_quantif-TP.tsv; join -t $'\t' -1 1 -2 3 <( sort -k1,1 /mnt/d/Anaquin-main/data/transcriptome/rnasequin_isoforms_2.5.tsv ) <( sort -k3,3 ../GFF-${file%*.txt}_quantif-TP.tsv ) | cut -d$'\t' -f5,6,1,7,8,9,3 > ../GFF-${file%*.txt}_quantif-TP2.tsv; join -1 1 -2 5 <( sort -k1,1 $file ) <( sort -k5,5 ../${file%*.txt}.tmap ) | sed 's/ /\t/g' | awk '$8!="=" && $8!="c" && $8!="u" {print $7"\t0\tFP\t"$1"\t"$8"\t"$4"\t"$5"\t"}' - > ../GFF-${file%*.txt}_quantif-FP.tsv; cat headerGFF.tsv ../GFF-${file%*.txt}_quantif-TP2.tsv ../GFF-${file%*.txt}_quantif-FP.tsv > ../GFF-${file%*.txt}_quantif.tsv; rm ../GFF-${file%*.txt}_quantif-FP.tsv ../GFF-${file%*.txt}_quantif-TP.tsv ../GFF-${file%*.txt}_quantif-TP2.tsv; done
#rm headerGFF.tsv


#SQ3
#echo -e "Sequin\tMixA\tType\tTranscript\tClasscode1\tClasscode2\tTPM\tnbReads" > headerSQ3.tsv
#for file in *.txt ; do join -1 1 -2 1 <(sort -k1,1 $file ) <( sort -k1,1 ../${file%*.txt}_classification.txt | sed 's/%7C/|/g' - ) | sed 's/ /\t/g' | awk '$10=="full-splice_match" || ( $10=="incomplete-splice_match" && $19!="intron_retention" ) {print "TP\t"$1"\t"$12"\t"$10"\t"$19"\t"$4"\t"$5}' - > ../SQ3-${file%*.txt}_quantif-TP.tsv; join -t $'\t' -1 1 -2 3 <( sort -k1,1 /mnt/d/Anaquin-main/data/transcriptome/rnasequin_isoforms_2.5.tsv ) <( sort -k3,3 ../SQ3-${file%*.txt}_quantif-TP.tsv ) | cut -d$'\t' -f5,6,1,7,8,9,10,3 > ../SQ3-${file%*.txt}_quantif-TP2.tsv; join -1 1 -2 1 <( sort -k1,1 $file ) <( sort -k1,1 ../${file%*.txt}_classification.txt | sed 's/%7C/|/g' - ) | sed 's/ /\t/g' | awk '$10!="full-splice_match" && ( $10!="incomplete-splice_match" || $19=="intron_retention" ) && $10!="intergenic" {print $12"\t0\tFP\t"$1"\t"$10"\t"$19"\t"$4"\t"$5}' - > ../SQ3-${file%*.txt}_quantif-FP.tsv; cat headerSQ3.tsv ../SQ3-${file%*.txt}_quantif-TP2.tsv ../SQ3-${file%*.txt}_quantif-FP.tsv > ../SQ3-${file%*.txt}_quantif.tsv; rm ../SQ3-${file%*.txt}_quantif-FP.tsv ../SQ3-${file%*.txt}_quantif-TP.tsv ../SQ3-${file%*.txt}_quantif-TP2.tsv; done
#rm headerSQ3.tsv

#Plot figure with R
#RStudio
library(ggplot2)
library(ggpubr)
library(cowplot)

setwd("D:/Assembly-Benchmarking/Fig5-Quantification/v2")

#All sequins abundances sum = 15196.5066 for TPM conversion

GFF_LSK114_Bambu <- read.csv("GFF-LSK114_Bambu_quantif-update.csv", header=T, sep=";")
GFF_LSK114_Bambu$nbreadslog <- log10(GFF_LSK114_Bambu$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_LSK114_Bambu$MixAlog <- log10(((GFF_LSK114_Bambu$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_LSK114_Bambu_x <- subset(subset(GFF_LSK114_Bambu, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_LSK114_Bambu_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_LSK114_Bambu_plot <- ggplot(GFF_LSK114_Bambu,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_Bambu_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_Bambu_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_Bambu_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_Bambu_noRef <- read.csv("GFF-LSK114_Bambu-noRef_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_Bambu_noRef$nbreadslog <- log10(GFF_LSK114_Bambu_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_Bambu_noRef$MixAlog <- log10(((GFF_LSK114_Bambu_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_Bambu_noRef_x <- subset(subset(GFF_LSK114_Bambu_noRef, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_Bambu_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_Bambu_noRef_plot <- ggplot(GFF_LSK114_Bambu_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_Bambu_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_Bambu_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_Bambu_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

GFF_LSK114_isoQuant <- read.csv("GFF-LSK114_isoQuant_quantif-update.csv", header=T, sep=";")
GFF_LSK114_isoQuant$nbreadslog <- log10(GFF_LSK114_isoQuant$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_LSK114_isoQuant$MixAlog <- log10(((GFF_LSK114_isoQuant$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_LSK114_isoQuant_x <- subset(subset(GFF_LSK114_isoQuant, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_LSK114_isoQuant_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_LSK114_isoQuant_plot <- ggplot(GFF_LSK114_isoQuant,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_isoQuant_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_isoQuant_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_isoQuant_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_isoQuant_noRef <- read.csv("GFF-LSK114_isoQuant-noRef_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_isoQuant_noRef$nbreadslog <- log10(GFF_LSK114_isoQuant_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_isoQuant_noRef$MixAlog <- log10(((GFF_LSK114_isoQuant_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_isoQuant_noRef_x <- subset(subset(GFF_LSK114_isoQuant_noRef, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_isoQuant_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_isoQuant_noRef_plot <- ggplot(GFF_LSK114_isoQuant_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_isoQuant_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_isoQuant_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_isoQuant_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_Flair <- read.csv("GFF-LSK114_FLAIR_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_Flair$nbreadslog <- log10(GFF_LSK114_Flair$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_Flair$MixAlog <- log10(((GFF_LSK114_Flair$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_Flair_x <- subset(subset(GFF_LSK114_Flair, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_Flair_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_Flair_plot <- ggplot(GFF_LSK114_Flair,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_Flair_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_Flair_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_Flair_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

GFF_LSK114_Flair_noRef <- read.csv("GFF-LSK114_FLAIR-noRef_quantif-update.csv", header=T, sep=";")
GFF_LSK114_Flair_noRef$nbreadslog <- log10(GFF_LSK114_Flair_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_LSK114_Flair_noRef$MixAlog <- log10(((GFF_LSK114_Flair_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_LSK114_Flair_noRef_x <- subset(subset(GFF_LSK114_Flair_noRef, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_LSK114_Flair_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_LSK114_Flair_noRef_plot <- ggplot(GFF_LSK114_Flair_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_Flair_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_Flair_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_Flair_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_FLAMES <- read.csv("GFF-LSK114_FLAMES-filt_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_FLAMES$nbreadslog <- log10(GFF_LSK114_FLAMES$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_FLAMES$MixAlog <- log10(((GFF_LSK114_FLAMES$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_FLAMES_x <- subset(subset(GFF_LSK114_FLAMES, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_FLAMES_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_FLAMES_plot <- ggplot(GFF_LSK114_FLAMES,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_FLAMES_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_FLAMES_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_FLAMES_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_isONclust <- read.csv("GFF-LSK114_isONclust_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_isONclust$nbreadslog <- log10(GFF_LSK114_isONclust$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_isONclust$MixAlog <- log10(((GFF_LSK114_isONclust$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_isONclust_x <- subset(subset(GFF_LSK114_isONclust, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_isONclust_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_isONclust_plot <- ggplot(GFF_LSK114_isONclust,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_isONclust_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_isONclust_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_isONclust_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_isONclust2 <- read.csv("GFF-LSK114_isONclust2_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_isONclust2$nbreadslog <- log10(GFF_LSK114_isONclust2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_isONclust2$MixAlog <- log10(((GFF_LSK114_isONclust2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_isONclust2_x <- subset(subset(GFF_LSK114_isONclust2, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_isONclust2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_isONclust2_plot <- ggplot(GFF_LSK114_isONclust2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_isONclust2_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_isONclust2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_isONclust2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_Mandalorion <- read.csv("GFF-LSK114_Mandalorion_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_Mandalorion$nbreadslog <- log10(GFF_LSK114_Mandalorion$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_Mandalorion$MixAlog <- log10(((GFF_LSK114_Mandalorion$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_Mandalorion_x <- subset(subset(GFF_LSK114_Mandalorion, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_Mandalorion_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_Mandalorion_plot <- ggplot(GFF_LSK114_Mandalorion,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_Mandalorion_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_Mandalorion_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_Mandalorion_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_RATTLE <- read.csv("GFF-LSK114_RATTLE_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_RATTLE$nbreadslog <- log10(GFF_LSK114_RATTLE$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_RATTLE$MixAlog <- log10(((GFF_LSK114_RATTLE$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_RATTLE_x <- subset(subset(GFF_LSK114_RATTLE, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_RATTLE_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_RATTLE_plot <- ggplot(GFF_LSK114_RATTLE,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_RATTLE_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_RATTLE_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_RATTLE_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_RNAbloom <- read.csv("GFF-LSK114_RNAbloom_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_RNAbloom$nbreadslog <- log10(GFF_LSK114_RNAbloom$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_RNAbloom$MixAlog <- log10(((GFF_LSK114_RNAbloom$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_RNAbloom_x <- subset(subset(GFF_LSK114_RNAbloom, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_RNAbloom_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_RNAbloom_plot <- ggplot(GFF_LSK114_RNAbloom,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_RNAbloom_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_RNAbloom_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_RNAbloom_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_RNAbloom2 <- read.csv("GFF-LSK114_RNAbloom2_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_RNAbloom2$nbreadslog <- log10(GFF_LSK114_RNAbloom2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_RNAbloom2$MixAlog <- log10(((GFF_LSK114_RNAbloom2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_RNAbloom2_x <- subset(subset(GFF_LSK114_RNAbloom2, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_RNAbloom2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_RNAbloom2_plot <- ggplot(GFF_LSK114_RNAbloom2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_RNAbloom2_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_RNAbloom2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_RNAbloom2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_Stringtie2 <- read.csv("GFF-LSK114_Stringtie2_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_Stringtie2$nbreadslog <- log10(GFF_LSK114_Stringtie2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_Stringtie2$MixAlog <- log10(((GFF_LSK114_Stringtie2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_Stringtie2_x <- subset(subset(GFF_LSK114_Stringtie2, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_Stringtie2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_Stringtie2_plot <- ggplot(GFF_LSK114_Stringtie2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_Stringtie2_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_Stringtie2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_Stringtie2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

GFF_LSK114_Stringtie2_noRef <- read.csv("GFF-LSK114_Stringtie2-noRef_quantif-update.csv", header=T, sep=";")
GFF_LSK114_Stringtie2_noRef$nbreadslog <- log10(GFF_LSK114_Stringtie2_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_LSK114_Stringtie2_noRef$MixAlog <- log10(((GFF_LSK114_Stringtie2_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_LSK114_Stringtie2_noRef_x <- subset(subset(GFF_LSK114_Stringtie2_noRef, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_LSK114_Stringtie2_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_LSK114_Stringtie2_noRef_plot <- ggplot(GFF_LSK114_Stringtie2_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_Stringtie2_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_Stringtie2_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_Stringtie2_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_TALON <- read.csv("GFF-LSK114_TALON_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_TALON$nbreadslog <- log10(GFF_LSK114_TALON$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_TALON$MixAlog <- log10(((GFF_LSK114_TALON$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_TALON_x <- subset(subset(GFF_LSK114_TALON, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_TALON_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_TALON_plot <- ggplot(GFF_LSK114_TALON,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_TALON_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_TALON_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_TALON_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_LSK114_TALON_reco <- read.csv("GFF-LSK114_TALON_reco_quantif.tsv", header=T, sep="\t")
#GFF_LSK114_TALON_reco$nbreadslog <- log10(GFF_LSK114_TALON_reco$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_LSK114_TALON_reco$MixAlog <- log10(((GFF_LSK114_TALON_reco$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_LSK114_TALON_reco_x <- subset(subset(GFF_LSK114_TALON_reco, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_LSK114_TALON_reco_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_LSK114_TALON_reco_plot <- ggplot(GFF_LSK114_TALON_reco,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_LSK114_TALON_reco_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_LSK114_TALON_reco_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_LSK114_TALON_reco_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

SQ3_LSK114_Bambu <- read.csv("SQ3-LSK114_Bambu_quantif-update.csv", header=T, sep=";")
SQ3_LSK114_Bambu$nbreadslog <- log10(SQ3_LSK114_Bambu$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_LSK114_Bambu$MixAlog <- log10(((SQ3_LSK114_Bambu$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_LSK114_Bambu_x <- subset(subset(SQ3_LSK114_Bambu, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_LSK114_Bambu_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_LSK114_Bambu_plot <- ggplot(SQ3_LSK114_Bambu,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_Bambu_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_Bambu_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_Bambu_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_Bambu_noRef <- read.csv("SQ3-LSK114_Bambu-noRef_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_Bambu_noRef$nbreadslog <- log10(SQ3_LSK114_Bambu_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_Bambu_noRef$MixAlog <- log10(((SQ3_LSK114_Bambu_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_Bambu_noRef_x <- subset(subset(SQ3_LSK114_Bambu_noRef, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_Bambu_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_Bambu_noRef_plot <- ggplot(SQ3_LSK114_Bambu_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_Bambu_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_Bambu_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_Bambu_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

SQ3_LSK114_isoQuant <- read.csv("SQ3-LSK114_isoQuant_quantif-update.csv", header=T, sep=";")
SQ3_LSK114_isoQuant$nbreadslog <- log10(SQ3_LSK114_isoQuant$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_LSK114_isoQuant$MixAlog <- log10(((SQ3_LSK114_isoQuant$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_LSK114_isoQuant_x <- subset(subset(SQ3_LSK114_isoQuant, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_LSK114_isoQuant_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_LSK114_isoQuant_plot <- ggplot(SQ3_LSK114_isoQuant,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_isoQuant_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_isoQuant_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_isoQuant_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_isoQuant_noRef <- read.csv("SQ3-LSK114_isoQuant-noRef_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_isoQuant_noRef$nbreadslog <- log10(SQ3_LSK114_isoQuant_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_isoQuant_noRef$MixAlog <- log10(((SQ3_LSK114_isoQuant_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_isoQuant_noRef_x <- subset(subset(SQ3_LSK114_isoQuant_noRef, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_isoQuant_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_isoQuant_noRef_plot <- ggplot(SQ3_LSK114_isoQuant_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_isoQuant_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_isoQuant_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_isoQuant_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_Flair <- read.csv("SQ3-LSK114_FLAIR_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_Flair$nbreadslog <- log10(SQ3_LSK114_Flair$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_Flair$MixAlog <- log10(((SQ3_LSK114_Flair$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_Flair_x <- subset(subset(SQ3_LSK114_Flair, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_Flair_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_Flair_plot <- ggplot(SQ3_LSK114_Flair,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_Flair_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_Flair_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_Flair_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

SQ3_LSK114_Flair_noRef <- read.csv("SQ3-LSK114_FLAIR-noRef_quantif-update.csv", header=T, sep=";")
SQ3_LSK114_Flair_noRef$nbreadslog <- log10(SQ3_LSK114_Flair_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_LSK114_Flair_noRef$MixAlog <- log10(((SQ3_LSK114_Flair_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_LSK114_Flair_noRef_x <- subset(subset(SQ3_LSK114_Flair_noRef, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_LSK114_Flair_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_LSK114_Flair_noRef_plot <- ggplot(SQ3_LSK114_Flair_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_Flair_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_Flair_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_Flair_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_FLAMES <- read.csv("SQ3-LSK114_FLAMES-filt_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_FLAMES$nbreadslog <- log10(SQ3_LSK114_FLAMES$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_FLAMES$MixAlog <- log10(((SQ3_LSK114_FLAMES$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_FLAMES_x <- subset(subset(SQ3_LSK114_FLAMES, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_FLAMES_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_FLAMES_plot <- ggplot(SQ3_LSK114_FLAMES,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_FLAMES_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_FLAMES_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_FLAMES_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_isONclust <- read.csv("SQ3-LSK114_isONclust_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_isONclust$nbreadslog <- log10(SQ3_LSK114_isONclust$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_isONclust$MixAlog <- log10(((SQ3_LSK114_isONclust$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_isONclust_x <- subset(subset(SQ3_LSK114_isONclust, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_isONclust_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_isONclust_plot <- ggplot(SQ3_LSK114_isONclust,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_isONclust_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_isONclust_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_isONclust_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_isONclust2 <- read.csv("SQ3-LSK114_isONclust2_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_isONclust2$nbreadslog <- log10(SQ3_LSK114_isONclust2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_isONclust2$MixAlog <- log10(((SQ3_LSK114_isONclust2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_isONclust2_x <- subset(subset(SQ3_LSK114_isONclust2, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_isONclust2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_isONclust2_plot <- ggplot(SQ3_LSK114_isONclust2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_isONclust2_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_isONclust2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_isONclust2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_Mandalorion <- read.csv("SQ3-LSK114_Mandalorion_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_Mandalorion$nbreadslog <- log10(SQ3_LSK114_Mandalorion$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_Mandalorion$MixAlog <- log10(((SQ3_LSK114_Mandalorion$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_Mandalorion_x <- subset(subset(SQ3_LSK114_Mandalorion, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_Mandalorion_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_Mandalorion_plot <- ggplot(SQ3_LSK114_Mandalorion,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_Mandalorion_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_Mandalorion_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_Mandalorion_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_RATTLE <- read.csv("SQ3-LSK114_RATTLE_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_RATTLE$nbreadslog <- log10(SQ3_LSK114_RATTLE$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_RATTLE$MixAlog <- log10(((SQ3_LSK114_RATTLE$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_RATTLE_x <- subset(subset(SQ3_LSK114_RATTLE, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_RATTLE_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_RATTLE_plot <- ggplot(SQ3_LSK114_RATTLE,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_RATTLE_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_RATTLE_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_RATTLE_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_RNAbloom <- read.csv("SQ3-LSK114_RNAbloom_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_RNAbloom$nbreadslog <- log10(SQ3_LSK114_RNAbloom$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_RNAbloom$MixAlog <- log10(((SQ3_LSK114_RNAbloom$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_RNAbloom_x <- subset(subset(SQ3_LSK114_RNAbloom, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_RNAbloom_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_RNAbloom_plot <- ggplot(SQ3_LSK114_RNAbloom,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_RNAbloom_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_RNAbloom_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_RNAbloom_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_RNAbloom2 <- read.csv("SQ3-LSK114_RNAbloom2_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_RNAbloom2$nbreadslog <- log10(SQ3_LSK114_RNAbloom2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_RNAbloom2$MixAlog <- log10(((SQ3_LSK114_RNAbloom2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_RNAbloom2_x <- subset(subset(SQ3_LSK114_RNAbloom2, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_RNAbloom2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_RNAbloom2_plot <- ggplot(SQ3_LSK114_RNAbloom2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_RNAbloom2_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_RNAbloom2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_RNAbloom2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_Stringtie2 <- read.csv("SQ3-LSK114_Stringtie2_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_Stringtie2$nbreadslog <- log10(SQ3_LSK114_Stringtie2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_Stringtie2$MixAlog <- log10(((SQ3_LSK114_Stringtie2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_Stringtie2_x <- subset(subset(SQ3_LSK114_Stringtie2, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_Stringtie2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_Stringtie2_plot <- ggplot(SQ3_LSK114_Stringtie2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_Stringtie2_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_Stringtie2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_Stringtie2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

SQ3_LSK114_Stringtie2_noRef <- read.csv("SQ3-LSK114_Stringtie2-noRef_quantif-update.csv", header=T, sep=";")
SQ3_LSK114_Stringtie2_noRef$nbreadslog <- log10(SQ3_LSK114_Stringtie2_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_LSK114_Stringtie2_noRef$MixAlog <- log10(((SQ3_LSK114_Stringtie2_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_LSK114_Stringtie2_noRef_x <- subset(subset(SQ3_LSK114_Stringtie2_noRef, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_LSK114_Stringtie2_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_LSK114_Stringtie2_noRef_plot <- ggplot(SQ3_LSK114_Stringtie2_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_Stringtie2_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_Stringtie2_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_Stringtie2_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_TALON <- read.csv("SQ3-LSK114_TALON_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_TALON$nbreadslog <- log10(SQ3_LSK114_TALON$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_TALON$MixAlog <- log10(((SQ3_LSK114_TALON$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_TALON_x <- subset(subset(SQ3_LSK114_TALON, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_TALON_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_TALON_plot <- ggplot(SQ3_LSK114_TALON,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_TALON_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_TALON_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_TALON_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_LSK114_TALON_reco <- read.csv("SQ3-LSK114_TALON_reco_quantif.tsv", header=T, sep="\t")
#SQ3_LSK114_TALON_reco$nbreadslog <- log10(SQ3_LSK114_TALON_reco$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_LSK114_TALON_reco$MixAlog <- log10(((SQ3_LSK114_TALON_reco$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_LSK114_TALON_reco_x <- subset(subset(SQ3_LSK114_TALON_reco, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_LSK114_TALON_reco_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_LSK114_TALON_reco_plot <- ggplot(SQ3_LSK114_TALON_reco,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_LSK114_TALON_reco_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_LSK114_TALON_reco_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_LSK114_TALON_reco_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

GFF_RNA004_Bambu <- read.csv("GFF-RNA004_Bambu_quantif-update.csv", header=T, sep=";")
GFF_RNA004_Bambu$nbreadslog <- log10(GFF_RNA004_Bambu$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_RNA004_Bambu$MixAlog <- log10(((GFF_RNA004_Bambu$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_RNA004_Bambu_x <- subset(subset(GFF_RNA004_Bambu, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_RNA004_Bambu_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_RNA004_Bambu_plot <- ggplot(GFF_RNA004_Bambu,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_Bambu_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_Bambu_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_Bambu_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_Bambu_noRef <- read.csv("GFF-RNA004_Bambu-noRef_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_Bambu_noRef$nbreadslog <- log10(GFF_RNA004_Bambu_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_Bambu_noRef$MixAlog <- log10(((GFF_RNA004_Bambu_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_Bambu_noRef_x <- subset(subset(GFF_RNA004_Bambu_noRef, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_Bambu_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_Bambu_noRef_plot <- ggplot(GFF_RNA004_Bambu_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_Bambu_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_Bambu_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_Bambu_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

GFF_RNA004_isoQuant <- read.csv("GFF-RNA004_isoQuant_quantif-update.csv", header=T, sep=";")
GFF_RNA004_isoQuant$nbreadslog <- log10(GFF_RNA004_isoQuant$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_RNA004_isoQuant$MixAlog <- log10(((GFF_RNA004_isoQuant$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_RNA004_isoQuant_x <- subset(subset(GFF_RNA004_isoQuant, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_RNA004_isoQuant_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_RNA004_isoQuant_plot <- ggplot(GFF_RNA004_isoQuant,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_isoQuant_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_isoQuant_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_isoQuant_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_isoQuant_noRef <- read.csv("GFF-RNA004_isoQuant-noRef_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_isoQuant_noRef$nbreadslog <- log10(GFF_RNA004_isoQuant_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_isoQuant_noRef$MixAlog <- log10(((GFF_RNA004_isoQuant_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_isoQuant_noRef_x <- subset(subset(GFF_RNA004_isoQuant_noRef, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_isoQuant_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_isoQuant_noRef_plot <- ggplot(GFF_RNA004_isoQuant_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_isoQuant_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_isoQuant_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_isoQuant_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_Flair <- read.csv("GFF-RNA004_FLAIR_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_Flair$nbreadslog <- log10(GFF_RNA004_Flair$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_Flair$MixAlog <- log10(((GFF_RNA004_Flair$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_Flair_x <- subset(subset(GFF_RNA004_Flair, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_Flair_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_Flair_plot <- ggplot(GFF_RNA004_Flair,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_Flair_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_Flair_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_Flair_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

GFF_RNA004_Flair_noRef <- read.csv("GFF-RNA004_FLAIR-noRef_quantif-update.csv", header=T, sep=";")
GFF_RNA004_Flair_noRef$nbreadslog <- log10(GFF_RNA004_Flair_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_RNA004_Flair_noRef$MixAlog <- log10(((GFF_RNA004_Flair_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_RNA004_Flair_noRef_x <- subset(subset(GFF_RNA004_Flair_noRef, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_RNA004_Flair_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_RNA004_Flair_noRef_plot <- ggplot(GFF_RNA004_Flair_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_Flair_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_Flair_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_Flair_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_FLAMES <- read.csv("GFF-RNA004_FLAMES-filt_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_FLAMES$nbreadslog <- log10(GFF_RNA004_FLAMES$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_FLAMES$MixAlog <- log10(((GFF_RNA004_FLAMES$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_FLAMES_x <- subset(subset(GFF_RNA004_FLAMES, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_FLAMES_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_FLAMES_plot <- ggplot(GFF_RNA004_FLAMES,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_FLAMES_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_FLAMES_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_FLAMES_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_isONclust <- read.csv("GFF-RNA004_isONclust_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_isONclust$nbreadslog <- log10(GFF_RNA004_isONclust$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_isONclust$MixAlog <- log10(((GFF_RNA004_isONclust$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_isONclust_x <- subset(subset(GFF_RNA004_isONclust, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_isONclust_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_isONclust_plot <- ggplot(GFF_RNA004_isONclust,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_isONclust_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_isONclust_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_isONclust_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_Mandalorion <- read.csv("GFF-RNA004_Mandalorion_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_Mandalorion$nbreadslog <- log10(GFF_RNA004_Mandalorion$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_Mandalorion$MixAlog <- log10(((GFF_RNA004_Mandalorion$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_Mandalorion_x <- subset(subset(GFF_RNA004_Mandalorion, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_Mandalorion_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_Mandalorion_plot <- ggplot(GFF_RNA004_Mandalorion,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_Mandalorion_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_Mandalorion_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_Mandalorion_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_RATTLE <- read.csv("GFF-RNA004_RATTLE_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_RATTLE$nbreadslog <- log10(GFF_RNA004_RATTLE$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_RATTLE$MixAlog <- log10(((GFF_RNA004_RATTLE$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_RATTLE_x <- subset(subset(GFF_RNA004_RATTLE, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_RATTLE_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_RATTLE_plot <- ggplot(GFF_RNA004_RATTLE,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_RATTLE_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_RATTLE_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_RATTLE_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_RNAbloom <- read.csv("GFF-RNA004_RNAbloom_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_RNAbloom$nbreadslog <- log10(GFF_RNA004_RNAbloom$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_RNAbloom$MixAlog <- log10(((GFF_RNA004_RNAbloom$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_RNAbloom_x <- subset(subset(GFF_RNA004_RNAbloom, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_RNAbloom_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_RNAbloom_plot <- ggplot(GFF_RNA004_RNAbloom,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_RNAbloom_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_RNAbloom_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_RNAbloom_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_RNAbloom2 <- read.csv("GFF-RNA004_RNAbloom2_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_RNAbloom2$nbreadslog <- log10(GFF_RNA004_RNAbloom2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_RNAbloom2$MixAlog <- log10(((GFF_RNA004_RNAbloom2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_RNAbloom2_x <- subset(subset(GFF_RNA004_RNAbloom2, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_RNAbloom2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_RNAbloom2_plot <- ggplot(GFF_RNA004_RNAbloom2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_RNAbloom2_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_RNAbloom2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_RNAbloom2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_Stringtie2 <- read.csv("GFF-RNA004_Stringtie2_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_Stringtie2$nbreadslog <- log10(GFF_RNA004_Stringtie2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_Stringtie2$MixAlog <- log10(((GFF_RNA004_Stringtie2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_Stringtie2_x <- subset(subset(GFF_RNA004_Stringtie2, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_Stringtie2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_Stringtie2_plot <- ggplot(GFF_RNA004_Stringtie2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_Stringtie2_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_Stringtie2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_Stringtie2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

GFF_RNA004_Stringtie2_noRef <- read.csv("GFF-RNA004_Stringtie2-noRef_quantif-update.csv", header=T, sep=";")
GFF_RNA004_Stringtie2_noRef$nbreadslog <- log10(GFF_RNA004_Stringtie2_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
GFF_RNA004_Stringtie2_noRef$MixAlog <- log10(((GFF_RNA004_Stringtie2_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
GFF_RNA004_Stringtie2_noRef_x <- subset(subset(GFF_RNA004_Stringtie2_noRef, Type=="TP"), TPM>0 )
x <- cor(as.matrix(GFF_RNA004_Stringtie2_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
GFF_RNA004_Stringtie2_noRef_plot <- ggplot(GFF_RNA004_Stringtie2_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_Stringtie2_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_Stringtie2_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_Stringtie2_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_TALON <- read.csv("GFF-RNA004_TALON_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_TALON$nbreadslog <- log10(GFF_RNA004_TALON$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_TALON$MixAlog <- log10(((GFF_RNA004_TALON$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_TALON_x <- subset(subset(GFF_RNA004_TALON, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_TALON_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_TALON_plot <- ggplot(GFF_RNA004_TALON,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_TALON_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_TALON_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_TALON_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#GFF_RNA004_TALON_reco <- read.csv("GFF-RNA004_TALON_reco_quantif.tsv", header=T, sep="\t")
#GFF_RNA004_TALON_reco$nbreadslog <- log10(GFF_RNA004_TALON_reco$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#GFF_RNA004_TALON_reco$MixAlog <- log10(((GFF_RNA004_TALON_reco$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#GFF_RNA004_TALON_reco_x <- subset(subset(GFF_RNA004_TALON_reco, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(GFF_RNA004_TALON_reco_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#GFF_RNA004_TALON_reco_plot <- ggplot(GFF_RNA004_TALON_reco,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=GFF_RNA004_TALON_reco_x,se=FALSE, method=lm) + stat_regline_equation(data=GFF_RNA004_TALON_reco_x,label.x=-0.8, label.y=5.75) + stat_cor(data=GFF_RNA004_TALON_reco_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

SQ3_RNA004_Bambu <- read.csv("SQ3-RNA004_Bambu_quantif-update.csv", header=T, sep=";")
SQ3_RNA004_Bambu$nbreadslog <- log10(SQ3_RNA004_Bambu$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_RNA004_Bambu$MixAlog <- log10(((SQ3_RNA004_Bambu$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_RNA004_Bambu_x <- subset(subset(SQ3_RNA004_Bambu, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_RNA004_Bambu_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_RNA004_Bambu_plot <- ggplot(SQ3_RNA004_Bambu,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_Bambu_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_Bambu_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_Bambu_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_Bambu_noRef <- read.csv("SQ3-RNA004_Bambu-noRef_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_Bambu_noRef$nbreadslog <- log10(SQ3_RNA004_Bambu_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_Bambu_noRef$MixAlog <- log10(((SQ3_RNA004_Bambu_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_Bambu_noRef_x <- subset(subset(SQ3_RNA004_Bambu_noRef, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_Bambu_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_Bambu_noRef_plot <- ggplot(SQ3_RNA004_Bambu_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_Bambu_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_Bambu_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_Bambu_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

SQ3_RNA004_isoQuant <- read.csv("SQ3-RNA004_isoQuant_quantif-update.csv", header=T, sep=";")
SQ3_RNA004_isoQuant$nbreadslog <- log10(SQ3_RNA004_isoQuant$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_RNA004_isoQuant$MixAlog <- log10(((SQ3_RNA004_isoQuant$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_RNA004_isoQuant_x <- subset(subset(SQ3_RNA004_isoQuant, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_RNA004_isoQuant_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_RNA004_isoQuant_plot <- ggplot(SQ3_RNA004_isoQuant,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_isoQuant_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_isoQuant_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_isoQuant_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_isoQuant_noRef <- read.csv("SQ3-RNA004_isoQuant-noRef_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_isoQuant_noRef$nbreadslog <- log10(SQ3_RNA004_isoQuant_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_isoQuant_noRef$MixAlog <- log10(((SQ3_RNA004_isoQuant_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_isoQuant_noRef_x <- subset(subset(SQ3_RNA004_isoQuant_noRef, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_isoQuant_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_isoQuant_noRef_plot <- ggplot(SQ3_RNA004_isoQuant_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_isoQuant_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_isoQuant_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_isoQuant_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_Flair <- read.csv("SQ3-RNA004_FLAIR_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_Flair$nbreadslog <- log10(SQ3_RNA004_Flair$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_Flair$MixAlog <- log10(((SQ3_RNA004_Flair$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_Flair_x <- subset(subset(SQ3_RNA004_Flair, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_Flair_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_Flair_plot <- ggplot(SQ3_RNA004_Flair,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_Flair_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_Flair_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_Flair_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

SQ3_RNA004_Flair_noRef <- read.csv("SQ3-RNA004_FLAIR-noRef_quantif-update.csv", header=T, sep=";")
SQ3_RNA004_Flair_noRef$nbreadslog <- log10(SQ3_RNA004_Flair_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_RNA004_Flair_noRef$MixAlog <- log10(((SQ3_RNA004_Flair_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_RNA004_Flair_noRef_x <- subset(subset(SQ3_RNA004_Flair_noRef, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_RNA004_Flair_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_RNA004_Flair_noRef_plot <- ggplot(SQ3_RNA004_Flair_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_Flair_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_Flair_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_Flair_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_FLAMES <- read.csv("SQ3-RNA004_FLAMES-filt_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_FLAMES$nbreadslog <- log10(SQ3_RNA004_FLAMES$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_FLAMES$MixAlog <- log10(((SQ3_RNA004_FLAMES$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_FLAMES_x <- subset(subset(SQ3_RNA004_FLAMES, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_FLAMES_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_FLAMES_plot <- ggplot(SQ3_RNA004_FLAMES,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_FLAMES_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_FLAMES_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_FLAMES_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_isONclust <- read.csv("SQ3-RNA004_isONclust_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_isONclust$nbreadslog <- log10(SQ3_RNA004_isONclust$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_isONclust$MixAlog <- log10(((SQ3_RNA004_isONclust$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_isONclust_x <- subset(subset(SQ3_RNA004_isONclust, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_isONclust_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_isONclust_plot <- ggplot(SQ3_RNA004_isONclust,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_isONclust_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_isONclust_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_isONclust_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_Mandalorion <- read.csv("SQ3-RNA004_Mandalorion_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_Mandalorion$nbreadslog <- log10(SQ3_RNA004_Mandalorion$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_Mandalorion$MixAlog <- log10(((SQ3_RNA004_Mandalorion$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_Mandalorion_x <- subset(subset(SQ3_RNA004_Mandalorion, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_Mandalorion_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_Mandalorion_plot <- ggplot(SQ3_RNA004_Mandalorion,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_Mandalorion_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_Mandalorion_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_Mandalorion_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_RATTLE <- read.csv("SQ3-RNA004_RATTLE_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_RATTLE$nbreadslog <- log10(SQ3_RNA004_RATTLE$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_RATTLE$MixAlog <- log10(((SQ3_RNA004_RATTLE$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_RATTLE_x <- subset(subset(SQ3_RNA004_RATTLE, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_RATTLE_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_RATTLE_plot <- ggplot(SQ3_RNA004_RATTLE,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_RATTLE_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_RATTLE_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_RATTLE_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_RNAbloom <- read.csv("SQ3-RNA004_RNAbloom_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_RNAbloom$nbreadslog <- log10(SQ3_RNA004_RNAbloom$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_RNAbloom$MixAlog <- log10(((SQ3_RNA004_RNAbloom$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_RNAbloom_x <- subset(subset(SQ3_RNA004_RNAbloom, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_RNAbloom_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_RNAbloom_plot <- ggplot(SQ3_RNA004_RNAbloom,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_RNAbloom_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_RNAbloom_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_RNAbloom_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_RNAbloom2 <- read.csv("SQ3-RNA004_RNAbloom2_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_RNAbloom2$nbreadslog <- log10(SQ3_RNA004_RNAbloom2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_RNAbloom2$MixAlog <- log10(((SQ3_RNA004_RNAbloom2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_RNAbloom2_x <- subset(subset(SQ3_RNA004_RNAbloom2, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_RNAbloom2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_RNAbloom2_plot <- ggplot(SQ3_RNA004_RNAbloom2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_RNAbloom2_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_RNAbloom2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_RNAbloom2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_Stringtie2 <- read.csv("SQ3-RNA004_Stringtie2_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_Stringtie2$nbreadslog <- log10(SQ3_RNA004_Stringtie2$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_Stringtie2$MixAlog <- log10(((SQ3_RNA004_Stringtie2$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_Stringtie2_x <- subset(subset(SQ3_RNA004_Stringtie2, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_Stringtie2_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_Stringtie2_plot <- ggplot(SQ3_RNA004_Stringtie2,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_Stringtie2_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_Stringtie2_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_Stringtie2_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

SQ3_RNA004_Stringtie2_noRef <- read.csv("SQ3-RNA004_Stringtie2-noRef_quantif-update.csv", header=T, sep=";")
SQ3_RNA004_Stringtie2_noRef$nbreadslog <- log10(SQ3_RNA004_Stringtie2_noRef$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
SQ3_RNA004_Stringtie2_noRef$MixAlog <- log10(((SQ3_RNA004_Stringtie2_noRef$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
SQ3_RNA004_Stringtie2_noRef_x <- subset(subset(SQ3_RNA004_Stringtie2_noRef, Type=="TP"), TPM>0 )
x <- cor(as.matrix(SQ3_RNA004_Stringtie2_noRef_x[,8:9]), method = c("pearson", "kendall", "spearman"))
print(x)
SQ3_RNA004_Stringtie2_noRef_plot <- ggplot(SQ3_RNA004_Stringtie2_noRef,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_Stringtie2_noRef_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_Stringtie2_noRef_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_Stringtie2_noRef_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_TALON <- read.csv("SQ3-RNA004_TALON_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_TALON$nbreadslog <- log10(SQ3_RNA004_TALON$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_TALON$MixAlog <- log10(((SQ3_RNA004_TALON$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_TALON_x <- subset(subset(SQ3_RNA004_TALON, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_TALON_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_TALON_plot <- ggplot(SQ3_RNA004_TALON,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_TALON_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_TALON_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_TALON_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#SQ3_RNA004_TALON_reco <- read.csv("SQ3-RNA004_TALON_reco_quantif.tsv", header=T, sep="\t")
#SQ3_RNA004_TALON_reco$nbreadslog <- log10(SQ3_RNA004_TALON_reco$TPM+0.1) #Take TPMs instead of nb of  reads and compare to sequins in TPMs
#SQ3_RNA004_TALON_reco$MixAlog <- log10(((SQ3_RNA004_TALON_reco$MixA*1000000)/15196.5066)+0.1) #Convert sequins relative abundances to TPMs
#SQ3_RNA004_TALON_reco_x <- subset(subset(SQ3_RNA004_TALON_reco, Type=="TP"), TPM>0 )
#x <- cor(as.matrix(SQ3_RNA004_TALON_reco_x[,8:9]), method = c("pearson", "kendall", "spearman"))
#print(x)
#SQ3_RNA004_TALON_reco_plot <- ggplot(SQ3_RNA004_TALON_reco,aes(y=nbreadslog, x=MixAlog)) + geom_abline(colour="gold", linetype = "longdash") + geom_vline(xintercept=0.82390874, colour="black", linetype = "dashed") + geom_point(aes(colour=factor(Type)),show.legend = FALSE) + geom_smooth(data=SQ3_RNA004_TALON_reco_x,se=FALSE, method=lm) + stat_regline_equation(data=SQ3_RNA004_TALON_reco_x,label.x=-0.8, label.y=5.75) + stat_cor(data=SQ3_RNA004_TALON_reco_x,aes(label=..rr.label..), label.x=-0.8, label.y=5.25) + xlab("Expected,log10(TPM+0.1)") + ylab("Observed,log10(TPM+0.1)") + theme_bw() + xlim(-1.1,6) + ylim(-1.1,6) + scale_color_manual(values = c("TP" = "midnightblue","FP" = "mediumvioletred","FPtp" = "mediumvioletred"))

#Plot LSK114
#plot_grid(GFF_LSK114_Bambu_plot,GFF_LSK114_isoQuant_plot,GFF_LSK114_Stringtie2_plot,GFF_LSK114_Flair_plot,GFF_LSK114_FLAMES_plot,GFF_LSK114_Mandalorion_plot,GFF_LSK114_TALON_plot,GFF_LSK114_TALON_reco_plot,SQ3_LSK114_Bambu_plot,SQ3_LSK114_isoQuant_plot,SQ3_LSK114_Stringtie2_plot,SQ3_LSK114_Flair_plot,SQ3_LSK114_FLAMES_plot,SQ3_LSK114_Mandalorion_plot,SQ3_LSK114_TALON_plot,SQ3_LSK114_TALON_reco_plot, nrow=2)
#PDF 20x5 inches
#plot_grid(GFF_LSK114_Bambu_noRef_plot,GFF_LSK114_isoQuant_noRef_plot,GFF_LSK114_Stringtie2_noRef_plot,GFF_LSK114_Flair_noRef_plot,SQ3_LSK114_Bambu_noRef_plot,SQ3_LSK114_isoQuant_noRef_plot,SQ3_LSK114_Stringtie2_noRef_plot,SQ3_LSK114_Flair_noRef_plot, nrow=2)
#PDF 10x5 inches
#plot_grid(GFF_LSK114_RATTLE_plot,GFF_LSK114_RNAbloom_plot,GFF_LSK114_RNAbloom2_plot,GFF_LSK114_isONclust_plot,GFF_LSK114_isONclust2_plot,SQ3_LSK114_RATTLE_plot,SQ3_LSK114_RNAbloom_plot,SQ3_LSK114_RNAbloom2_plot,SQ3_LSK114_isONclust_plot,SQ3_LSK114_isONclust2_plot, nrow=2)
#PDF 12x5 inches

#Plot RNA004
#plot_grid(GFF_RNA004_Bambu_plot,GFF_RNA004_isoQuant_plot,GFF_RNA004_Stringtie2_plot,GFF_RNA004_Flair_plot,GFF_RNA004_FLAMES_plot,GFF_RNA004_Mandalorion_plot,GFF_RNA004_TALON_plot,GFF_RNA004_TALON_reco_plot,SQ3_RNA004_Bambu_plot,SQ3_RNA004_isoQuant_plot,SQ3_RNA004_Stringtie2_plot,SQ3_RNA004_Flair_plot,SQ3_RNA004_FLAMES_plot,SQ3_RNA004_Mandalorion_plot,SQ3_RNA004_TALON_plot,SQ3_RNA004_TALON_reco_plot, nrow=2)
#PDF 20x5 inches
#plot_grid(GFF_RNA004_Bambu_noRef_plot,GFF_RNA004_isoQuant_noRef_plot,GFF_RNA004_Stringtie2_noRef_plot,GFF_RNA004_Flair_noRef_plot,SQ3_RNA004_Bambu_noRef_plot,SQ3_RNA004_isoQuant_noRef_plot,SQ3_RNA004_Stringtie2_noRef_plot,SQ3_RNA004_Flair_noRef_plot, nrow=2)
#PDF 10x5 inches
#plot_grid(GFF_RNA004_RATTLE_plot,GFF_RNA004_RNAbloom_plot,GFF_RNA004_RNAbloom2_plot,GFF_RNA004_isONclust_plot,SQ3_RNA004_RATTLE_plot,SQ3_RNA004_RNAbloom_plot,SQ3_RNA004_RNAbloom2_plot,SQ3_RNA004_isONclust_plot, nrow=2)
#PDF 10x5 inches


plot_grid(GFF_LSK114_Bambu_plot,GFF_LSK114_isoQuant_plot,GFF_LSK114_Flair_noRef_plot,GFF_LSK114_Stringtie2_noRef_plot,SQ3_LSK114_Bambu_plot,SQ3_LSK114_isoQuant_plot,SQ3_LSK114_Flair_noRef_plot,SQ3_LSK114_Stringtie2_noRef_plot,GFF_RNA004_Bambu_plot,GFF_RNA004_isoQuant_plot,GFF_RNA004_Flair_noRef_plot,GFF_RNA004_Stringtie2_noRef_plot,SQ3_RNA004_Bambu_plot,SQ3_RNA004_isoQuant_plot,SQ3_RNA004_Flair_noRef_plot,SQ3_RNA004_Stringtie2_noRef_plot,nrow=4)
#pdf 10x10