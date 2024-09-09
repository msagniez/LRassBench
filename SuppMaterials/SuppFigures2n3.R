## loading libraries
library(rtracklayer)
library(tidyverse)
library(stringr)
library(ggtranscript)
library(magrittr)
library(dplyr)
library(ggplot2)
library("wiggleplotr")
library("GenomicRanges")
library("GenomicFeatures")
library("biomaRt")
library(readxl)
library(gridExtra)
library(ggtext)
library(data.table)
library(reshape2)
library(cowplot)

## loading data
GFF <- read.table("~/GFF-TPnovelFPFN-update_AllCases.csv", header=F, sep=",")
colnames(GFF) <- c("Case","Assembly","Completeness","Section","TP","Novel","FP","FN")
SQ3 <- read.table("~/SQ3-TPnovelFPFN-update_AllCases.csv", header=F, sep=",")
colnames(SQ3) <- c("Case","Assembly","Completeness","Section","TP","Novel","FP","FN")

## averaging TP, novel, FP and FN values obtained from GFFcompare and SQANTI3
df <- GFF
df$TP <- (GFF$TP + SQ3$TP)/2
df$Novel <- (GFF$Novel + SQ3$Novel)/2
df$FP <- (GFF$FP + SQ3$FP)/2
df$FN <- (GFF$FN + SQ3$FN)/2
## renaming Section to S, to shorten labels
df$Section <- str_replace_all(df$Section, "Section","S")
## cases are selected to ensure the presence of all cases
dfsub <- subset(df, Case=="chrIS_full_half_null" | Case=="chrIS_half_null_full" | Case=="chrIS_null_full_half")
## rearranging dataframe for ggplot input
dfsub <- melt(dfsub, id.vars=c("Case","Assembly","Completeness","Section"))

dffull <- subset(dfsub, Completeness=="full")
dfpartial <- subset(dfsub, Completeness=="partial")
dfnull <- subset(dfsub, Completeness=="null")

## determine the order for each df
# levels=unique(dffull$Assembly[order(dffull$TP, decreasing = T)])
# levels=unique(dffull$Assembly[order(dffull$Novel, decreasing = T)]) 
level_full=c("GroundTruth","Bambu","Stringtie2","TALON","TALON_reco","isoQuant","Mandalorion","FLAIR","FLAMES")
dffull$Assembly <- factor(dffull$Assembly, levels=level_full)
level_partial=c("GroundTruth","Bambu","Stringtie2","Mandalorion","TALON","isoQuant","TALON_reco","FLAIR","FLAMES")
dfpartial$Assembly <- factor(dfpartial$Assembly, levels=level_partial)
level_null=c("GroundTruth","Mandalorion","Stringtie2","TALON","isoQuant","TALON_reco","Bambu","FLAIR","FLAMES")
dfnull$Assembly <- factor(dfnull$Assembly, levels=level_null)

## initialising the colours for ggplot2
colour_vals <- c("TP" = "midnightblue","Novel" = "mediumvioletred","FP"="turquoise2","FN" = "grey")

## arrange according to Assembly for numbering stacked bar
dffull_fill <- dffull %>%arrange(Assembly,Section, rev(variable))
dfpartial_fill <- dfpartial %>%arrange(Assembly,Section, rev(variable))
dfnull_fill <- dfnull %>%arrange(Assembly,Section, rev(variable))

## set distance of numerical label - along y-axis
dffull_fill <- dffull_fill %>% group_by(Section,Assembly,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 
dfpartial_fill <- dfpartial_fill %>% group_by(Section,Assembly,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 
dfnull_fill <- dfnull_fill %>% group_by(Section,Assembly,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 

## plotting using ggplot2
Pl1 <- ggplot(dffull_fill, aes(x=Section,y=pct,fill=variable)) + geom_bar(position="stack",stat="identity") + facet_grid(~ Assembly, switch="x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  geom_text(data=subset(dffull_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + scale_x_discrete(position = "top") + geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + scale_fill_manual(values = colour_vals) + ggtitle("Full")
Pl2 <- ggplot(dfpartial_fill, aes(x=Section,y=pct,fill=variable)) + geom_bar(position="stack",stat="identity") + facet_grid(~ Assembly, switch="x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  geom_text(data=subset(dfpartial_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + scale_x_discrete(position = "top") + geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + scale_fill_manual(values = colour_vals) + ggtitle("Partial")
Pl3 <- ggplot(dfnull_fill, aes(x=Section,y=pct,fill=variable)) + geom_bar(position="stack",stat="identity") + facet_grid(~ Assembly, switch="x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  geom_text(data=subset(dfnull_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + scale_x_discrete(position = "top") + geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + scale_fill_manual(values = colour_vals) + ggtitle("Null")
## saving legend separately
legend <- get_legend(Pl1)
## modifying background, legend etc
Pl1 <- Pl1 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
    theme(plot.background = element_rect(color = "black"))
Pl2 <- Pl2 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(plot.background = element_rect(color = "black"))
Pl3 <- Pl3 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(plot.background = element_rect(color = "black"))
## plot: Supplemental Figure 2
Plot.2 <- plot_grid(Pl1,Pl2,Pl3, nrow=1, rel_widths=c(2,2,2)) 
## saving to pdf
pdf(file="~/SupplementalFigure2.pdf", width = 20, height = 10)
plot_grid(Plot.2,legend, rel_widths=c(6,0.25)) 
dev.off()

## visualising without FP
dfsub_noFP <- dfsub %>% subset(variable!="FP")
dffull_noFP <- subset(dfsub_noFP, Completeness=="full")
dfpartial_noFP <- subset(dfsub_noFP, Completeness=="partial")
dfnull_noFP <- subset(dfsub_noFP, Completeness=="null")

dffull_noFP$Assembly <- factor(dffull_noFP$Assembly, levels=level_full)
dfpartial_noFP$Assembly <- factor(dfpartial_noFP$Assembly, levels=level_partial)
dfnull_noFP$Assembly <- factor(dfnull_noFP$Assembly, levels=level_null)

## arrange according to Assembly for numbering stacked bar
dffull_fill <- dffull_noFP %>%arrange(Assembly,Section, rev(variable))
dfpartial_fill <- dfpartial_noFP %>%arrange(Assembly,Section, rev(variable))
dfnull_fill <- dfnull_noFP %>%arrange(Assembly,Section, rev(variable))

## set distance of numerical label - along y-axis
dffull_fill <- dffull_fill %>% group_by(Section,Assembly,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 
dfpartial_fill <- dfpartial_fill %>% group_by(Section,Assembly,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 
dfnull_fill <- dfnull_fill %>% group_by(Section,Assembly,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 

## plotting using ggplot2
Pl1 <- ggplot(dffull_fill, aes(x=Section,y=pct,fill=variable)) + geom_bar(position="stack",stat="identity") + facet_grid(~ Assembly, switch="x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  geom_text(data=subset(dffull_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + scale_x_discrete(position = "top") + geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + scale_fill_manual(values = colour_vals) + ggtitle("Full")
Pl2 <- ggplot(dfpartial_fill, aes(x=Section,y=pct,fill=variable)) + geom_bar(position="stack",stat="identity") + facet_grid(~ Assembly, switch="x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  geom_text(data=subset(dfpartial_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + scale_x_discrete(position = "top") + geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + scale_fill_manual(values = colour_vals) + ggtitle("Partial")
Pl3 <- ggplot(dfnull_fill, aes(x=Section,y=pct,fill=variable)) + geom_bar(position="stack",stat="identity") + facet_grid(~ Assembly, switch="x") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  geom_text(data=subset(dfnull_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + scale_x_discrete(position = "top") + geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + scale_fill_manual(values = colour_vals) + ggtitle("Null")
## saving legend separately
legend <- get_legend(Pl1)
## modifying background, legend etc
Pl1 <- Pl1 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
    theme(plot.background = element_rect(color = "black"))
Pl2 <- Pl2 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(plot.background = element_rect(color = "black"))
Pl3 <- Pl3 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    theme(plot.background = element_rect(color = "black"))
## plot: Supplemental Figure 3
Plot.3 <- plot_grid(Pl1,Pl2,Pl3, nrow=1, rel_widths=c(2,2,2)) 
pdf(file="~/SupplementalFigure3.pdf", width = 30, height = 10)
plot_grid(Plot.3,legend, rel_widths=c(6,0.25))
dev.off()
