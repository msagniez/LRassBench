library(ggplot2)
library(ggrepel)
library(cowplot)

setwd("D:/Assembly-Benchmarking/Fig2-PrecisionSensitivity")

#sub150k
Sub150k <- read.csv("sub150k/Sub150k_MainSuppB.csv", sep=";", header = T)
LSK114_chrIS <- subset(Sub150k, Evaluator == "GFFcompare" & Dataset == "LSK114_sequins")
PCS109_chrIS <- subset(Sub150k, Evaluator == "GFFcompare" & Dataset == "LSK109_sequins")
RNA002 <- subset(Sub150k, Evaluator == "GFFcompare" & Dataset == "RNA002_sequins")
RNA004 <- subset(Sub150k, Evaluator == "GFFcompare" & Dataset == "RNA004_sequins")
SQ3_LSK114_chrIS <- subset(Sub150k, Evaluator == "SQANTI3" & Dataset == "LSK114_sequins")
SQ3_PCS109_chrIS <- subset(Sub150k, Evaluator == "SQANTI3" & Dataset == "LSK109_sequins")
SQ3_RNA002 <- subset(Sub150k, Evaluator == "SQANTI3" & Dataset == "RNA002_sequins")
SQ3_RNA004 <- subset(Sub150k, Evaluator == "SQANTI3" & Dataset == "RNA004_sequins")

#sub150k average
AVG_LSK114_chrIS <- subset(Sub150k, Evaluator == "Average" & Dataset == "LSK114_sequins")
AVG_PCS109_chrIS <- subset(Sub150k, Evaluator == "Average" & Dataset == "LSK109_sequins")
AVG_RNA002 <- subset(Sub150k, Evaluator == "Average" & Dataset == "RNA002_sequins")
AVG_RNA004 <- subset(Sub150k, Evaluator == "Average" & Dataset == "RNA004_sequins")

#sub40k
Sub40k <- read.csv("sub40k/Sub40k_MainSuppA.csv", sep=";", header = T)
S40_GFF_chrIS <- subset(Sub40k, Evaluator == "GFFcompare" & Dataset == "LSK114_sequins")
S40_SQ3_chrIS <- subset(Sub40k, Evaluator == "SQANTI3" & Dataset == "LSK114_sequins")
S40_SQ3_SIRV <- subset(Sub40k, Evaluator == "SQANTI3" & Dataset == "LSK114_SIRVs")
S40_GFF_SIRV <- subset(Sub40k, Evaluator == "GFFcompare" & Dataset == "LSK114_SIRVs")

#sub40k average
S40_AVG_chrIS <- subset(Sub40k, Evaluator == "Average" & Dataset == "LSK114_sequins")
S40_AVG_SIRV <- subset(Sub40k, Evaluator == "Average" & Dataset == "LSK114_SIRVs")

#Plot sub150k
chrIS_R9 <- ggplot(PCS109_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9)) 
chrIS_R10 <- ggplot(LSK114_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
chrIS_RNA002 <- ggplot(RNA002, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
chrIS_RNA004 <- ggplot(RNA004, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))


SQ3_chrIS_R9 <- ggplot(SQ3_PCS109_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
SQ3_chrIS_R10 <- ggplot(SQ3_LSK114_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
SQ3_chrIS_RNA002 <- ggplot(SQ3_RNA002, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
SQ3_chrIS_RNA004 <- ggplot(SQ3_RNA004, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))

#Plot sub40k
GFF_chrIS_plot <- ggplot(S40_GFF_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
SQ3_chrIS_plot <- ggplot(S40_SQ3_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
SQ3_SIRV_plot <- ggplot(S40_SQ3_SIRV, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
GFF_SIRV_plot <- ggplot(S40_GFF_SIRV, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))



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

AVG_chrIS_R9 <- ggplot(AVG_PCS109_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
AVG_chrIS_LSK114 <- ggplot(AVG_LSK114_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
AVG_chrIS_RNA002 <- ggplot(AVG_RNA002, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
AVG_chrIS_RNA004 <- ggplot(AVG_RNA004, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
AVG_chrIS_plot <- ggplot(S40_AVG_chrIS, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
AVG_SIRV_plot <- ggplot(S40_AVG_SIRV, aes(fill=Tool,x=Precision, y=Sensitivity, label=Assembly)) + geom_point(aes(colour = Tool, shape = Method), size=4) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + scale_shape_manual(values=c(19, 17, 15,9))
legend <- get_legend(AVG_chrIS_R9)
AVG_chrIS_R9 <- AVG_chrIS_R9 + theme(legend.position = "none")
AVG_chrIS_LSK114 <- AVG_chrIS_LSK114 + theme(legend.position = "none")
AVG_chrIS_RNA002 <- AVG_chrIS_RNA002 + theme(legend.position = "none")
AVG_chrIS_RNA004 <- AVG_chrIS_RNA004 + theme(legend.position = "none")
AVG_chrIS_plot <- AVG_chrIS_plot + theme(legend.position = "none")
AVG_SIRV_plot <- AVG_SIRV_plot + theme(legend.position = "none")

#Main figure
plot_grid(AVG_SIRV_plot,AVG_chrIS_R9,AVG_chrIS_LSK114,AVG_chrIS_plot,AVG_chrIS_RNA002,AVG_chrIS_RNA004, nrow = 2, rel_widths=c(6,6,6))

#pdf 10x6

