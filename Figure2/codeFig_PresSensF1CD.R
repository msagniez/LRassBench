library(ggplot2)
library(ggrepel)
library(cowplot)

setwd("D:/Assembly-Benchmarking/Fig2-PrecisionSensitivity")

#sub40k
inpt40 <- read.csv("sub40k/GFF-SQ3-PresRecallF1_average.csv", sep = ";", header = T)
inpt40$Dataset <- factor(inpt40$Dataset, levels=c("LSK114_SIRVs","LSK114_sequins"))

Prec40 <- subset(inpt40,Vtype=="Precision")
Rec40 <- subset(inpt40,Vtype=="Sensitivity")
F140 <- subset(inpt40,Vtype=="F1")
Prec40_DN <- subset(Prec40,Method=="De novo")
Prec40_G <- subset(Prec40,Method=="Guided")
Prec40_AI <- subset(Prec40,Method=="Ab initio")
Rec40_DN <- subset(Rec40,Method=="De novo")
Rec40_G <- subset(Rec40,Method=="Guided")
Rec40_AI <- subset(Rec40,Method=="Ab initio")
F140_DN <- subset(F140,Method=="De novo")
F140_G <- subset(F140,Method=="Guided")
F140_AI <- subset(F140,Method=="Ab initio")

Prec40_guided <- ggplot(Prec40_G, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly != "TALON_reco")) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly == "TALON_reco")) +  scale_linetype_manual(values = c("Bambu" = "solid", "FLAIR" = "solid", "isoQuant" = "solid", "TALON" = "solid", "FLAMES" = "solid", "Mandalorion" = "solid", "Stringtie2" = "solid", "TALON_reco" = "longdash")) + ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Precision")
Prec40_DeNovo <- ggplot(Prec40_DN, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Precision")
Prec40_AbInitio <- ggplot(Prec40_AI, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Precision")

Rec40_guided <- ggplot(Rec40_G, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly != "TALON_reco")) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly == "TALON_reco")) +  scale_linetype_manual(values = c("Bambu" = "solid", "FLAIR" = "solid", "isoQuant" = "solid", "TALON" = "solid", "FLAMES" = "solid", "Mandalorion" = "solid", "Stringtie2" = "solid", "TALON_reco" = "longdash")) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Sensitivity")
Rec40_DeNovo <- ggplot(Rec40_DN, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Sensitivity")
Rec40_AbInitio <- ggplot(Rec40_AI, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Sensitivity")

F140_guided <- ggplot(F140_G, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly != "TALON_reco")) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly == "TALON_reco")) +  scale_linetype_manual(values = c("Bambu" = "solid", "FLAIR" = "solid", "isoQuant" = "solid", "TALON" = "solid", "FLAMES" = "solid", "Mandalorion" = "solid", "Stringtie2" = "solid", "TALON_reco" = "longdash")) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("F1")
F140_DeNovo <- ggplot(F140_DN, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("F1")
F140_AbInitio <- ggplot(F140_AI, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("F1")

legend <- plot_grid(get_legend(Prec40_guided),get_legend(Prec40_DeNovo),get_legend(Prec40_AbInitio),ncol=1)

Prec40_DeNovo <- Prec40_DeNovo + theme(legend.position = "none")
Prec40_guided <- Prec40_guided + theme(legend.position = "none")
Prec40_AbInitio <- Prec40_AbInitio + theme(legend.position = "none")
Rec40_DeNovo <- Rec40_DeNovo + theme(legend.position = "none")
Rec40_guided <- Rec40_guided + theme(legend.position = "none")
Rec40_AbInitio <- Rec40_AbInitio + theme(legend.position = "none")
F140_DeNovo <- F140_DeNovo + theme(legend.position = "none")
F140_guided <- F140_guided + theme(legend.position = "none")
F140_AbInitio <- F140_AbInitio + theme(legend.position = "none")



#sub150k
inpt150 <- read.csv("sub150k/GFF-SQ3-PresRecallF1_average.csv", sep = ";", header = T)
inpt150$Dataset <- factor(inpt150$Dataset, levels=c("LSK109_sequins","LSK114_sequins","RNA002_sequins","RNA004_sequins"))

Prec150 <- subset(inpt150,Vtype=="Precision")
Rec150 <- subset(inpt150,Vtype=="Sensitivity")
F1150 <- subset(inpt150,Vtype=="F1")
Prec150_DN <- subset(Prec150,Method=="De novo")
Prec150_G <- subset(Prec150,Method=="Guided")
Prec150_AI <- subset(Prec150,Method=="Ab initio")
Rec150_DN <- subset(Rec150,Method=="De novo")
Rec150_G <- subset(Rec150,Method=="Guided")
Rec150_AI <- subset(Rec150,Method=="Ab initio")
F1150_DN <- subset(F1150,Method=="De novo")
F1150_G <- subset(F1150,Method=="Guided")
F1150_AI <- subset(F1150,Method=="Ab initio")

Prec150_guided <- ggplot(Prec150_G, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly != "TALON_reco")) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly == "TALON_reco")) +  scale_linetype_manual(values = c("Bambu" = "solid", "FLAIR" = "solid", "isoQuant" = "solid", "TALON" = "solid", "FLAMES" = "solid", "Mandalorion" = "solid", "Stringtie2" = "solid", "TALON_reco" = "longdash")) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Precision")
Prec150_DeNovo <- ggplot(Prec150_DN, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Precision")
Prec150_AbInitio <- ggplot(Prec150_AI, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Precision")

Rec150_guided <- ggplot(Rec150_G, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly != "TALON_reco")) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly == "TALON_reco")) +  scale_linetype_manual(values = c("Bambu" = "solid", "FLAIR" = "solid", "isoQuant" = "solid", "TALON" = "solid", "FLAMES" = "solid", "Mandalorion" = "solid", "Stringtie2" = "solid", "TALON_reco" = "longdash")) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Sensitivity")
Rec150_DeNovo <- ggplot(Rec150_DN, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Sensitivity")
Rec150_AbInitio <- ggplot(Rec150_AI, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("Sensitivity")

F1150_guided <- ggplot(F1150_G, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly != "TALON_reco")) + geom_line(aes(colour = Tool, linetype = Assembly), data = ~ subset(., Assembly == "TALON_reco")) +  scale_linetype_manual(values = c("Bambu" = "solid", "FLAIR" = "solid", "isoQuant" = "solid", "TALON" = "solid", "FLAMES" = "solid", "Mandalorion" = "solid", "Stringtie2" = "solid", "TALON_reco" = "longdash")) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("F1")
F1150_DeNovo <- ggplot(F1150_DN, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("F1")
F1150_AbInitio <- ggplot(F1150_AI, aes(fill=Tool, x=Dataset, y=Value, group=Assembly)) + geom_point(aes(colour = Tool)) + geom_line(aes(colour = Tool)) +  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_manual(values = c("Bambu" = "brown","FLAIR"="violet","isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","RATTLE"="springgreen2","RNAbloom"="coral2","RNAbloom2"="orangered3","Stringtie2"="palegreen4","TALON"="wheat2","FLAMES"="blue2","Mandalorion"="darkgoldenrod1")) + xlab("") + ylab("F1")

legend <- plot_grid(get_legend(Prec150_guided),get_legend(Prec150_DeNovo),get_legend(Prec150_AbInitio),ncol=1)

Prec150_DeNovo <- Prec150_DeNovo + theme(legend.position = "none")
Prec150_guided <- Prec150_guided + theme(legend.position = "none")
Prec150_AbInitio <- Prec150_AbInitio + theme(legend.position = "none")
Rec150_DeNovo <- Rec150_DeNovo + theme(legend.position = "none")
Rec150_guided <- Rec150_guided + theme(legend.position = "none")
Rec150_AbInitio <- Rec150_AbInitio + theme(legend.position = "none")
F1150_DeNovo <- F1150_DeNovo + theme(legend.position = "none")
F1150_guided <- F1150_guided + theme(legend.position = "none")
F1150_AbInitio <- F1150_AbInitio + theme(legend.position = "none")

plot_grid(Prec40_guided,Prec40_DeNovo,Prec40_AbInitio,Prec150_guided,Prec150_DeNovo,Prec150_AbInitio,Rec40_guided,Rec40_DeNovo,Rec40_AbInitio,Rec150_guided,Rec150_DeNovo,Rec150_AbInitio,F140_guided,F140_DeNovo,F140_AbInitio,F1150_guided,F1150_DeNovo,F1150_AbInitio, ncol=6, rel_widths = c(3.5,3.5,3.5,6,6,6))


#pdf 12x7.5 inches (landscape)