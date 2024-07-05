library(ggplot2)
library(cowplot)

TPFP <- read.csv("GFF-PrecisionSensitivity-update.csv", sep=";", header=T)
TPFP$Sensitivity <- TPFP$TP / (TPFP$TP + TPFP$FN )
TPFP$Precision <- TPFP$TP / (TPFP$TP + TPFP$FP )

Guided <- subset(TPFP, Method=="Guided" | Method=="RawReads")
DeNovo <- subset(TPFP, Method=="De novo" | Method=="RawReads")
AbInitio <- subset(TPFP, Method=="Ab initio" | Method=="RawReads")

Gd <- ggplot(Guided, aes(fill=Tool,x=Precision, y=Sensitivity, label=Subset)) + geom_point(aes(colour = Tool, shape = Method)) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown", "isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","Stringtie2"="palegreen4")) + scale_shape_manual(values=c(15,9)) 
DN <- ggplot(DeNovo, aes(fill=Tool,x=Precision, y=Sensitivity, label=Subset)) + geom_point(aes(colour = Tool, shape = Method)) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown", "isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","Stringtie2"="palegreen4")) + scale_shape_manual(values=c(17,9))
AI <- ggplot(AbInitio, aes(fill=Tool,x=Precision, y=Sensitivity, label=Subset)) + geom_point(aes(colour = Tool, shape = Method)) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown", "isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","Stringtie2"="palegreen4")) + scale_shape_manual(values=c(19,9))

ledg <- ggplot(subset(TPFP,Subset=="Sub1k"), aes(fill=Tool,x=Precision, y=Sensitivity, label=Subset)) + geom_point(aes(colour = Tool, shape = Method)) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown", "isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","Stringtie2"="palegreen4")) + scale_shape_manual(values=c(19,17,15,9))
legend <- get_legend(ledg)

SQTPFP <- read.csv("SQ3-PrecisionSensitivity-update.csv", sep=";", header=T)
SQTPFP$Sensitivity <- SQTPFP$TP / (SQTPFP$TP + SQTPFP$FN )
SQTPFP$Precision <- SQTPFP$TP / (SQTPFP$TP + SQTPFP$FP )
SQGuided <- subset(SQTPFP, Method=="Guided" | Method=="RawReads")
SQDeNovo <- subset(SQTPFP, Method=="De novo" | Method=="RawReads")
SQAbInitio <- subset(SQTPFP, Method=="Ab initio" | Method=="RawReads")
SQGd <- ggplot(SQGuided, aes(fill=Tool,x=Precision, y=Sensitivity, label=Subset)) + geom_point(aes(colour = Tool, shape = Method)) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown", "isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","Stringtie2"="palegreen4")) + scale_shape_manual(values=c(15,9))
SQDN <- ggplot(SQDeNovo, aes(fill=Tool,x=Precision, y=Sensitivity, label=Subset)) + geom_point(aes(colour = Tool, shape = Method)) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown", "isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","Stringtie2"="palegreen4")) + scale_shape_manual(values=c(17,9))
SQAI <- ggplot(SQAbInitio, aes(fill=Tool,x=Precision, y=Sensitivity, label=Subset)) + geom_point(aes(colour = Tool, shape = Method)) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown", "isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","Stringtie2"="palegreen4")) + scale_shape_manual(values=c(19,9))

Gd <- Gd + theme(legend.position = "none")
SQGd <- SQGd + theme(legend.position = "none")
DN <- DN + theme(legend.position = "none")
SQDN <- SQDN + theme(legend.position = "none")
AI <- AI + theme(legend.position = "none")
SQAI <- SQAI + theme(legend.position = "none")

plot_grid(Gd,DN,AI,legend,SQGd,SQDN,SQAI,nrow=2,rel_widths = c(6,6,6,3))
#pdf 11x6 (supp figure)

#Plot average values for GFF and SQ3
#AVGTPFP <- read.csv("AVG-PrecisionSensitivity-update.csv", sep=";", header=T)
AVGTPFP <- TPFP
AVGTPFP$Sensitivity <- (TPFP$Sensitivity + SQTPFP$Sensitivity)/2
AVGTPFP$Precision <- (TPFP$Precision + SQTPFP$Precision)/2
AVGGuided <- subset(AVGTPFP, Method=="Guided" | Method=="RawReads")
AVGDeNovo <- subset(AVGTPFP, Method=="De novo" | Method=="RawReads")
AVGAbInitio <- subset(AVGTPFP, Method=="Ab initio" | Method=="RawReads")
AVGGd <- ggplot(AVGGuided, aes(fill=Tool,x=Precision, y=Sensitivity, label=Subset)) + geom_point(aes(colour = Tool, shape = Method)) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown", "isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","Stringtie2"="palegreen4")) + scale_shape_manual(values=c(15,9))
AVGDN <- ggplot(AVGDeNovo, aes(fill=Tool,x=Precision, y=Sensitivity, label=Subset)) + geom_point(aes(colour = Tool, shape = Method)) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown", "isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","Stringtie2"="palegreen4")) + scale_shape_manual(values=c(17,9))
AVGAI <- ggplot(AVGAbInitio, aes(fill=Tool,x=Precision, y=Sensitivity, label=Subset)) + geom_point(aes(colour = Tool, shape = Method)) + xlim(0,1) + ylim(0,1) + theme_bw() + scale_color_manual(values = c("Bambu" = "brown", "isONclust"="mediumpurple1","isONclust2"="purple","isoQuant"="turquoise2","Stringtie2"="palegreen4")) + scale_shape_manual(values=c(19,9))

AVGGd <- AVGGd + theme(legend.position = "none")
AVGDN <- AVGDN + theme(legend.position = "none")
AVGAI <- AVGAI + theme(legend.position = "none")

plot_grid(AVGGd,AVGDN,AVGAI,legend,nrow=1,rel_widths = c(6,6,6,3))
#pdf 11x3

plot_grid(Gd,DN,AI,legend,SQGd,SQDN,SQAI,"",AVGGd,AVGDN,AVGAI,legend, nrow=3,rel_widths = c(6,6,6,3))
