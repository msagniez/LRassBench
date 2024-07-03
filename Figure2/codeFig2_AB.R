library(ggplot2)
library(ggrepel)
library(cowplot)

setwd("D:/Assembly-Benchmarking/Fig2-PrecisionSensitivity")

#sub150k
LSK114_chrIS <- read.csv("sub150k/LSK114_chrIS-update.csv", sep=",", header = T)
LSK114_chrIS$Sensitivity <- LSK114_chrIS$TP / (LSK114_chrIS$TP + LSK114_chrIS$FN )
LSK114_chrIS$Precision <- LSK114_chrIS$TP / (LSK114_chrIS$TP + LSK114_chrIS$FP )
PCS109_chrIS <- read.csv("sub150k/PCS109_chrIS-update.csv", sep=",", header = T)
PCS109_chrIS$Sensitivity <- PCS109_chrIS$TP / (PCS109_chrIS$TP + PCS109_chrIS$FN )
PCS109_chrIS$Precision <- PCS109_chrIS$TP / (PCS109_chrIS$TP + PCS109_chrIS$FP )
RNA002 <- read.csv("sub150k/RNA002_chrIS-update.csv", sep=",", header = T)
RNA002$Sensitivity <- RNA002$TP / (RNA002$TP + RNA002$FN )
RNA002$Precision <- RNA002$TP / (RNA002$TP + RNA002$FP )
RNA004 <- read.csv("sub150k/RNA004_chrIS-update.csv", sep=",", header = T)
RNA004$Sensitivity <- RNA004$TP / (RNA004$TP + RNA004$FN )
RNA004$Precision <- RNA004$TP / (RNA004$TP + RNA004$FP )
SQ3_LSK114_chrIS <- read.csv("sub150k/SQ3_LSK114_chrIS-update.csv", sep=",", header = T)
SQ3_LSK114_chrIS$Sensitivity <- SQ3_LSK114_chrIS$TP / (SQ3_LSK114_chrIS$TP + SQ3_LSK114_chrIS$FN )
SQ3_LSK114_chrIS$Precision <- SQ3_LSK114_chrIS$TP / (SQ3_LSK114_chrIS$TP + SQ3_LSK114_chrIS$FP )
SQ3_PCS109_chrIS <- read.csv("sub150k/SQ3_PCS109_chrIS-update.csv", sep=",", header = T)
SQ3_PCS109_chrIS$Sensitivity <- SQ3_PCS109_chrIS$TP / (SQ3_PCS109_chrIS$TP + SQ3_PCS109_chrIS$FN )
SQ3_PCS109_chrIS$Precision <- SQ3_PCS109_chrIS$TP / (SQ3_PCS109_chrIS$TP + SQ3_PCS109_chrIS$FP )
SQ3_RNA002 <- read.csv("sub150k/SQ3_RNA002_chrIS-update.csv", sep=",", header = T)
SQ3_RNA002$Sensitivity <- SQ3_RNA002$TP / (SQ3_RNA002$TP + SQ3_RNA002$FN )
SQ3_RNA002$Precision <- SQ3_RNA002$TP / (SQ3_RNA002$TP + SQ3_RNA002$FP )
SQ3_RNA004 <- read.csv("sub150k/SQ3_RNA004_chrIS-update.csv", sep=",", header = T)
SQ3_RNA004$Sensitivity <- SQ3_RNA004$TP / (SQ3_RNA004$TP + SQ3_RNA004$FN )
SQ3_RNA004$Precision <- SQ3_RNA004$TP / (SQ3_RNA004$TP + SQ3_RNA004$FP )

#sub150k average
AVG_LSK114_chrIS <- LSK114_chrIS
AVG_LSK114_chrIS$Sensitivity <- (LSK114_chrIS$Sensitivity + SQ3_LSK114_chrIS$Sensitivity)/2
AVG_LSK114_chrIS$Precision <- (LSK114_chrIS$Precision + SQ3_LSK114_chrIS$Precision)/2
AVG_PCS109_chrIS <- PCS109_chrIS
AVG_PCS109_chrIS$Sensitivity <- (PCS109_chrIS$Sensitivity + SQ3_PCS109_chrIS$Sensitivity)/2
AVG_PCS109_chrIS$Precision <- (PCS109_chrIS$Precision + SQ3_PCS109_chrIS$Precision)/2
AVG_RNA002 <- RNA002
AVG_RNA002$Sensitivity <- (RNA002$Sensitivity + SQ3_RNA002$Sensitivity)/2
AVG_RNA002$Precision <- (RNA002$Precision + SQ3_RNA002$Precision)/2
AVG_RNA004 <- RNA004
AVG_RNA004$Sensitivity <- (RNA004$Sensitivity + SQ3_RNA004$Sensitivity)/2
AVG_RNA004$Precision <- (RNA004$Precision + SQ3_RNA004$Precision)/2

#sub40k
S40_GFF_chrIS <- read.csv("sub40k/chrIS-GFF-PrecisionSensitivity-update.csv", sep=",", header = T)
S40_GFF_chrIS$Sensitivity <- S40_GFF_chrIS$TP / (S40_GFF_chrIS$TP + S40_GFF_chrIS$FN )
S40_GFF_chrIS$Precision <- S40_GFF_chrIS$TP / (S40_GFF_chrIS$TP + S40_GFF_chrIS$FP )
S40_SQ3_chrIS <- read.csv("sub40k/chrIS-SQ3-PrecisionSensitivity-update.csv", sep=",", header = T)
S40_SQ3_chrIS$Sensitivity <- S40_SQ3_chrIS$TP / (S40_SQ3_chrIS$TP + S40_SQ3_chrIS$FN )
S40_SQ3_chrIS$Precision <- S40_SQ3_chrIS$TP / (S40_SQ3_chrIS$TP + S40_SQ3_chrIS$FP )
S40_SQ3_SIRV <- read.csv("sub40k/SIRV-SQ3-PrecisionSensitivity-update.csv", sep=",", header = T)
S40_SQ3_SIRV$Sensitivity <- S40_SQ3_SIRV$TP / (S40_SQ3_SIRV$TP + S40_SQ3_SIRV$FN )
S40_SQ3_SIRV$Precision <- S40_SQ3_SIRV$TP / (S40_SQ3_SIRV$TP + S40_SQ3_SIRV$FP )
S40_GFF_SIRV <- read.csv("sub40k/SIRV-GFF-PrecisionSensitivity-update.csv", sep=",", header = T)
S40_GFF_SIRV$Sensitivity <- S40_GFF_SIRV$TP / (S40_GFF_SIRV$TP + S40_GFF_SIRV$FN )
S40_GFF_SIRV$Precision <- S40_GFF_SIRV$TP / (S40_GFF_SIRV$TP + S40_GFF_SIRV$FP )

#sub40k average
S40_AVG_chrIS <- S40_GFF_chrIS
S40_AVG_chrIS$Sensitivity <- (S40_GFF_chrIS$Sensitivity + S40_SQ3_chrIS$Sensitivity)/2
S40_AVG_chrIS$Precision <- (S40_GFF_chrIS$Precision + S40_SQ3_chrIS$Precision)/2
S40_AVG_SIRV <- S40_GFF_SIRV
#Add GFF Talon results to the SQ3 table because TALON missing from the SQ3 evaluation ; showing only GFF results in final plot
tal <- S40_GFF_SIRV[which(S40_GFF_SIRV$Tool=='TALON'),]
S40_SQ3_SIRV_tal <- rbind(S40_SQ3_SIRV,tal)
S40_AVG_SIRV$Sensitivity <- (S40_GFF_SIRV$Sensitivity + S40_SQ3_SIRV_tal$Sensitivity)/2
S40_AVG_SIRV$Precision <- (S40_GFF_SIRV$Precision + S40_SQ3_SIRV_tal$Precision)/2

#Plot sub150k
chrIS_R9 <- ggplot(PCS109_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9)) 

chrIS_R10 <- ggplot(LSK114_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))

chrIS_RNA002 <- ggplot(RNA002, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))

chrIS_RNA004 <- ggplot(RNA004, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))


SQ3_chrIS_R9 <- ggplot(SQ3_PCS109_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))

SQ3_chrIS_R10 <- ggplot(SQ3_LSK114_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))

SQ3_chrIS_RNA002 <- ggplot(SQ3_RNA002, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))

SQ3_chrIS_RNA004 <- ggplot(SQ3_RNA004, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))

#Plot sub40k
GFF_chrIS_plot <- ggplot(S40_GFF_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
SQ3_chrIS_plot <- ggplot(S40_SQ3_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
SQ3_SIRV_plot <- ggplot(S40_SQ3_SIRV, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
GFF_SIRV_plot <- ggplot(S40_GFF_SIRV, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))



legend <- get_legend(chrIS_R9)
chrIS_R9 <- chrIS_R9 + theme(legend.position = "none")
chrIS_R10 <- chrIS_R10 + theme(legend.position = "none")
chrIS_RNA002 <- chrIS_RNA002 + theme(legend.position = "none")
chrIS_RNA004 <- chrIS_RNA004 + theme(legend.position = "none")
SQ3_chrIS_R9 <- SQ3_chrIS_R9 + theme(legend.position = "none")
SQ3_chrIS_R10 <- SQ3_chrIS_R10 + theme(legend.position = "none")
SQ3_chrIS_RNA002 <- SQ3_chrIS_RNA002 + theme(legend.position = "none")
SQ3_chrIS_RNA004 <- SQ3_chrIS_RNA004 + theme(legend.position = "none")
GFF_chrIS_plot <- GFF_chrIS_plot + theme(legend.position = "none")
SQ3_chrIS_plot <- SQ3_chrIS_plot + theme(legend.position = "none")
SQ3_SIRV_plot <- SQ3_SIRV_plot + theme(legend.position = "none")
GFF_SIRV_plot <- GFF_SIRV_plot + theme(legend.position = "none")

#Supp figure
SF <- plot_grid(GFF_SIRV_plot,SQ3_SIRV_plot,GFF_chrIS_plot,SQ3_chrIS_plot,chrIS_R9,chrIS_R10,chrIS_RNA002,chrIS_RNA004,SQ3_chrIS_R9,SQ3_chrIS_R10,SQ3_chrIS_RNA002,SQ3_chrIS_RNA004, nrow = 3, rel_widths=c(6,6,6,6,6))
plot_grid(SF,legend,nrow=1,rel_widths=c(18,3))
#pdf 13x8

AVG_chrIS_R9 <- ggplot(AVG_PCS109_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
AVG_chrIS_LSK114 <- ggplot(AVG_LSK114_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
AVG_chrIS_RNA002 <- ggplot(AVG_RNA002, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
AVG_chrIS_RNA004 <- ggplot(AVG_RNA004, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
AVG_chrIS_plot <- ggplot(S40_AVG_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
AVG_SIRV_plot <- ggplot(S40_AVG_SIRV, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("RawReads"="grey1","Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isONpipeline"="mediumpurple4","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
legend <- get_legend(AVG_chrIS_R9)
AVG_chrIS_R9 <- AVG_chrIS_R9 + theme(legend.position = "none")
AVG_chrIS_LSK114 <- AVG_chrIS_LSK114 + theme(legend.position = "none")
AVG_chrIS_RNA002 <- AVG_chrIS_RNA002 + theme(legend.position = "none")
AVG_chrIS_RNA004 <- AVG_chrIS_RNA004 + theme(legend.position = "none")
AVG_chrIS_plot <- AVG_chrIS_plot + theme(legend.position = "none")
AVG_SIRV_plot <- AVG_SIRV_plot + theme(legend.position = "none")

#Main figure
MF <- plot_grid(AVG_SIRV_plot,AVG_chrIS_R9,AVG_chrIS_LSK114,AVG_chrIS_plot,AVG_chrIS_RNA002,AVG_chrIS_RNA004, nrow = 2, rel_widths=c(6,6,6))
#plot_grid(MF,legend,nwow=1,rel_widths=c(18,3))
MF
#pdf 10x6

