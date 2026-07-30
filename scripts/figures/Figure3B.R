# ---- LIBRARIES ----
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


# ---- INPUT ----
scores_novelDisc <- read_excel("Supplemental_Tables.xlsx", sheet = "SuppTable5")
GFF <- scores_novelDisc %>% subset(Evaluator=="GFFcompare")
GFF <- GFF %>% dplyr::select(-Evaluator, -Dataset, -Section, -AbundanceLevels)
# SQ3 <- scores_novelDisc %>% subset(Evaluator=="SQANTI3")
# SQ3 <- SQ3 %>% dplyr::select(-Evaluator, -Dataset, -Section, -AbundanceLevels)


df_novelDisc <- GFF

## cases are selected to ensure the presence of all cases
dfsub <- subset(df_novelDisc, Case=="chrIS_full_half_null" | Case=="chrIS_half_null_full" | Case=="chrIS_null_full_half")
## rearranging dataframe for ggplot input
dfsub <- melt(dfsub, id.vars=c("Case","Tool","Completeness","SectionLabel"))

dffull <- subset(dfsub, Completeness=="full")
dfpartial <- subset(dfsub, Completeness=="partial")
dfnull <- subset(dfsub, Completeness=="null")


## determine the order for each df
# levels=unique(dffull$Tool[order(dffull$TP, decreasing = T)])
# levels=unique(dffull$Tool[order(dffull$Novel, decreasing = T)]) 

#tool_ranking_balanced ## from simple TP-FP

tool_ranking_balanced <- dffull %>%
  mutate(var2 = ifelse(variable == "Novel", "TP", as.character(variable))) %>%
  group_by(Tool, var2) %>%
  summarise(v = sum(value, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = var2, values_from = v, values_fill = 0) %>%
  mutate(
    score = TP - FP,
    rank  = rank(-score, ties.method = "min")
  ) %>%
  arrange(rank, desc(TP), FP)
tool_ranking_balanced
##level_full=c("GroundTruth","Bambu","isoQuant3","Stringtie3","Mandalorion","FLAMES","LRAA","FLAIR","TALON_reco","TALON" ) #ver1.pdf
level_full=tool_ranking_balanced$Tool #ver2.pdf
dffull$Tool <- factor(dffull$Tool, levels=level_full)

tool_ranking_balanced <- dfpartial %>%
  mutate(var2 = ifelse(variable == "Novel", "TP", as.character(variable))) %>%
  group_by(Tool, var2) %>%
  summarise(v = sum(value, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = var2, values_from = v, values_fill = 0) %>%
  mutate(
    score = TP - FP,
    rank  = rank(-score, ties.method = "min")
  ) %>%
  arrange(rank, desc(TP), FP)
level_partial=tool_ranking_balanced$Tool
dfpartial$Tool <- factor(dfpartial$Tool, levels=level_partial)

tool_ranking_balanced <- dfnull %>%
  mutate(var2 = ifelse(variable == "Novel", "TP", as.character(variable))) %>%
  group_by(Tool, var2) %>%
  summarise(v = sum(value, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = var2, values_from = v, values_fill = 0) %>%
  mutate(
    score = TP - FP,
    rank  = rank(-score, ties.method = "min")
  ) %>%
  arrange(rank, desc(TP), FP)
tool_ranking_balanced
# level_null=c("GroundTruth","Mandalorion","isoQuant3","Stringtie3","LRAA","Bambu","FLAIR","TALON_reco","TALON","FLAMES")
level_null=tool_ranking_balanced$Tool #ver3.pdf
dfnull$Tool <- factor(dfnull$Tool, levels=level_null)

## initialising the colours for ggplot2
colour_vals <- c("TP" = "midnightblue","Novel" = "mediumvioletred","FP"="turquoise2","FN" = "grey")


#### plotting per label on X-axis ; not Section.
## arrange according to Tool for numbering stacked bar
dffull_fill <- dffull %>%arrange(Tool,SectionLabel, rev(variable))
dfpartial_fill <- dfpartial %>%arrange(Tool,SectionLabel, rev(variable))
dfnull_fill <- dfnull %>%arrange(Tool,SectionLabel, rev(variable))

## set distance of numerical label - along y-axis
dffull_fill <- dffull_fill %>% group_by(SectionLabel,Tool,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 
dfpartial_fill <- dfpartial_fill %>% group_by(SectionLabel,Tool,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 
dfnull_fill <- dfnull_fill %>% group_by(SectionLabel,Tool,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 
## plotting using ggplot2
Pl1 <- ggplot(dffull_fill, aes(x=SectionLabel,y=pct,fill=variable)) + 
  geom_bar(position="stack",stat="identity") + 
  facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=subset(dffull_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L")) + 
  geom_hline(yintercept = 1) + 
  geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Full")

Pl2 <- ggplot(dfpartial_fill, aes(x=SectionLabel,y=pct,fill=variable)) + 
  geom_bar(position="stack",stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=subset(dfpartial_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L")) + 
  geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Partial")

Pl3 <- ggplot(dfnull_fill, aes(x=SectionLabel,y=pct,fill=variable)) + 
  geom_bar(position="stack",stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=subset(dfnull_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L")) + 
  geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Null")

## saving legend separately
legend <- get_legend(Pl1)
## modifying background, legend etc
Pl1 <- Pl1 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
  theme(plot.background = element_rect(color = "black"))
Pl2 <- Pl2 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(plot.background = element_rect(color = "black"))
Pl3 <- Pl3 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(plot.background = element_rect(color = "black"))
## plot: Supplemental Figure
Plot.1 <- plot_grid(Pl1,Pl2,Pl3, nrow=1, rel_widths=c(2,2,2)) 
## saving to pdf
pdf(file="Figures/Figure3_TPFPFN_PercentStackBar_panelsWithNum_GFFcompValues_23July2026.pdf", width = 20, height = 10)
plot_grid(Plot.1,legend, rel_widths=c(6,0.25)) 
dev.off()


## plotting using ggplot2
Pl1 <- ggplot(dffull_fill, aes(x=SectionLabel,y=n,fill=variable)) + 
  geom_bar(stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=subset(dffull_fill, n != 0), aes(label = n), position = position_stack(vjust = .95), colour = "white", size=2) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L"))+ 
  geom_hline(yintercept = 275) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Full")

Pl2 <- ggplot(dfpartial_fill, aes(x=SectionLabel,y=n,fill=variable)) + 
  geom_bar(stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=subset(dfpartial_fill, n != 0), aes(label = n), position = position_stack(vjust = .95), colour = "white", size=2) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L"))+ 
  geom_hline(yintercept = 275) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Partial")

Pl3 <- ggplot(dfnull_fill, aes(x=SectionLabel,y=n,fill=variable)) + 
  geom_bar(stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=subset(dfnull_fill, n != 0), aes(label = n), position = position_stack(vjust = .95), colour = "white", size=2) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L"))+ 
  geom_hline(yintercept = 275) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Null")

## saving legend separately
legend <- get_legend(Pl1)
## modifying background, legend etc
Pl1 <- Pl1 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
  theme(plot.background = element_rect(color = "black"))
Pl2 <- Pl2 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(plot.background = element_rect(color = "black"))
Pl3 <- Pl3 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(plot.background = element_rect(color = "black"))
## plot: Figure 3
Plot.2 <- plot_grid(Pl1,Pl2,Pl3, nrow=1, rel_widths=c(2,2,2)) 
## saving to pdf
pdf(file="Figures/Figure3_TPFPFN_NumStackBar_panelsWithNum_GFFcompValues_23July2026.pdf", width = 20, height = 10)
plot_grid(Plot.2,legend, rel_widths=c(6,0.25)) 
dev.off()

## without numbers - main figure
Pl2_1 <- ggplot(dfpartial_fill, aes(x=SectionLabel,y=n,fill=variable)) + 
  geom_bar(stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L"))+ 
  geom_hline(yintercept = 275) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Partial")
Pl2_1 <- Pl2_1 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
  theme(plot.background = element_rect(color = "black"))

pdf(file="Figures/Figure3B_main_TPFPFN_NumStackBar_partial_GFFcompValues_23July2026.pdf", width = 8, height = 12)
print(Pl2_1) 
dev.off()

## with numbers - main figure
Pl2_2 <- ggplot(dfpartial_fill, aes(x=SectionLabel,y=n,fill=variable)) + 
  geom_bar(stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) + 
  geom_text(data=subset(dfpartial_fill, n != 0), aes(label = n), position = position_stack(vjust = .95), colour = "white", size=3) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L"))+ 
  geom_hline(yintercept = 275) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Partial")
Pl2_2 <- Pl2_2 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
  theme(plot.background = element_rect(color = "black"))

pdf(file="Figures/Figure3_TPFPFN_NumStackBar_partialWithNum_GFFcompValues_23July2026.pdf", width = 8, height = 12)
print(Pl2_2) 
dev.off()




## visualising without FP
dfsub_noFP <- dfsub %>% subset(variable!="FP")
dffull_noFP <- subset(dfsub_noFP, Completeness=="full")
dfpartial_noFP <- subset(dfsub_noFP, Completeness=="partial")
dfnull_noFP <- subset(dfsub_noFP, Completeness=="null")

dffull_noFP$Tool <- factor(dffull_noFP$Tool, levels=level_full)
dfpartial_noFP$Tool <- factor(dfpartial_noFP$Tool, levels=level_partial)
dfnull_noFP$Tool <- factor(dfnull_noFP$Tool, levels=level_null)

## arrange according to Tool for numbering stacked bar
dffull_fill <- dffull_noFP %>%arrange(Tool,SectionLabel, rev(variable))
dfpartial_fill <- dfpartial_noFP %>%arrange(Tool,SectionLabel, rev(variable))
dfnull_fill <- dfnull_noFP %>%arrange(Tool,SectionLabel, rev(variable))

## set distance of numerical SectionLabel - along y-axis
dffull_fill <- dffull_fill %>% group_by(SectionLabel,Tool,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 
dfpartial_fill <- dfpartial_fill %>% group_by(SectionLabel,Tool,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 
dfnull_fill <- dfnull_fill %>% group_by(SectionLabel,Tool,variable) %>% summarise(n = sum(value, na.rm = TRUE)) %>% mutate(pct = prop.table(n)) 

## plotting using ggplot2
Pl1 <- ggplot(dffull_fill, aes(x=SectionLabel,y=pct,fill=variable)) + 
  geom_bar(position="stack",stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=subset(dffull_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L"))+ 
  geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Full")

Pl2 <- ggplot(dfpartial_fill, aes(x=SectionLabel,y=pct,fill=variable)) + 
  geom_bar(position="stack",stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=subset(dfpartial_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L"))+ 
  geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Partial")

Pl3 <- ggplot(dfnull_fill, aes(x=SectionLabel,y=pct,fill=variable)) + 
  geom_bar(position="stack",stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=subset(dfnull_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L"))+ 
  geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Null")

## saving legend separately
legend <- get_legend(Pl1)
## modifying background, legend etc
Pl1 <- Pl1 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
  theme(plot.background = element_rect(color = "black"))
Pl2 <- Pl2 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(plot.background = element_rect(color = "black"))
Pl3 <- Pl3 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none",axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
  theme(plot.background = element_rect(color = "black"))
## plot: Figure 3
Plot.3 <- plot_grid(Pl1,Pl2,Pl3, nrow=1, rel_widths=c(2,2,2)) 

pdf(file="Figures/Figure3_TPFN_PercentStackBar_panelsWithNum_GFFcompValues_23July2026.pdf", width = 30, height = 10)
plot_grid(Plot.3,legend, rel_widths=c(6,0.25))
dev.off()

## main Figure - without numbers
Pl2_3 <- ggplot(dfpartial_fill, aes(x=SectionLabel,y=pct,fill=variable)) + 
  geom_bar(position="stack",stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  scale_x_discrete(position = "top", limits = c("H", "M", "L"))+ 
  geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Partial")
Pl2_3 <- Pl2_3 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
  theme(plot.background = element_rect(color = "black"))

pdf(file="Figures/Figure3C_main_TPFN_PercentStackBar_partial_GFFcompValues_23July2026.pdf", width = 8, height = 10)
print(Pl2_3) 
dev.off()

## main figure - with numbers

Pl2_4 <- ggplot(dfpartial_fill, aes(x=SectionLabel,y=pct,fill=variable)) + 
  geom_bar(position="stack",stat="identity") + facet_grid(~ Tool, switch="x") + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust=0.5, size = 8),strip.placement='outside',panel.spacing.x = unit(0,"line"),panel.border = element_blank(),strip.text.x.bottom = element_text(angle = 90),strip.background = element_blank(),plot.title = element_text(hjust = 0.5)) +  
  geom_text(data=subset(dfpartial_fill, pct!=0), aes(label = scales::percent(pct, accuracy = .1)), position = position_stack(vjust = .5), colour = "white", size=2) + 
  scale_x_discrete(position = "top", limits = c("H", "M", "L"))+ 
  geom_hline(yintercept = 1) + geom_hline(yintercept = 0) + 
  scale_fill_manual(values = colour_vals) + 
  ggtitle("Partial")
Pl2_4 <- Pl2_4 + theme(axis.title.y=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
  theme(plot.background = element_rect(color = "black"))

pdf(file="Figures/Figure3_TPFN_PercentStackBar_partialWithNum_GFFcompValues_23July2026.pdf", width = 10, height = 12)
print(Pl2_4) 
dev.off()


# > sessionInfo()
# R version 4.5.1 (2025-06-13 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 11 x64 (build 26200)

# Matrix products: default
#   LAPACK version 3.12.1

# locale:
# [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
# [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.utf8    

# time zone: America/Toronto
# tzcode source: internal

# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] reshape2_1.4.4         data.table_1.17.8      ggtext_0.1.2           gridExtra_2.3         
#  [5] biomaRt_2.64.0         GenomicFeatures_1.60.0 AnnotationDbi_1.70.0   Biobase_2.68.0        
#  [9] wiggleplotr_1.32.0     magrittr_2.0.3         ggtranscript_1.0.0     rtracklayer_1.68.0    
# [13] GenomicRanges_1.60.0   GenomeInfoDb_1.44.2    IRanges_2.42.0         S4Vectors_0.46.0      
# [17] BiocGenerics_0.54.0    generics_0.1.4         ggrepel_0.9.6          readxl_1.4.5          
# [21] cowplot_1.2.0          lubridate_1.9.4        forcats_1.0.0          stringr_1.5.2         
# [25] dplyr_1.2.0            purrr_1.1.0            readr_2.1.5            tidyr_1.3.1           
# [29] tibble_3.3.0           ggplot2_4.0.2          tidyverse_2.0.0       

# loaded via a namespace (and not attached):
#  [1] DBI_1.2.3                   bitops_1.0-9                httr2_1.2.1                
#  [4] sandwich_3.1-1              rlang_1.1.7                 multcomp_1.4-28            
#  [7] matrixStats_1.5.0           compiler_4.5.1              RSQLite_2.4.3              
# [10] png_0.1-8                   vctrs_0.7.1                 pkgconfig_2.0.3            
# [13] crayon_1.5.3                fastmap_1.2.0               dbplyr_2.5.1               
# [16] XVector_0.48.0              Rsamtools_2.24.0            tzdb_0.5.0                 
# [19] UCSC.utils_1.4.0            bit_4.6.0                   modeltools_0.2-24          
# [22] cachem_1.1.0                jsonlite_2.0.0              progress_1.2.3             
# [25] blob_1.2.4                  DelayedArray_0.34.1         BiocParallel_1.42.1        
# [28] parallel_4.5.1              prettyunits_1.2.0           R6_2.6.1                   
# [31] coin_1.4-3                  stringi_1.8.7               RColorBrewer_1.1-3         
# [34] cellranger_1.1.0            Rcpp_1.1.0                  SummarizedExperiment_1.38.1
# [37] zoo_1.8-14                  Matrix_1.7-3                splines_4.5.1              
# [40] timechange_0.3.0            tidyselect_1.2.1            rstudioapi_0.17.1          
# [43] abind_1.4-8                 yaml_2.3.10                 codetools_0.2-20           
# [46] curl_7.0.0                  plyr_1.8.9                  lattice_0.22-7             
# [49] withr_3.0.2                 KEGGREST_1.48.1             S7_0.2.0                   
# [52] survival_3.8-3              BiocFileCache_2.16.2        xml2_1.4.0                 
# [55] Biostrings_2.76.0           filelock_1.0.3              pillar_1.11.0              
# [58] MatrixGenerics_1.20.0       RCurl_1.98-1.17             hms_1.1.3                  
# [61] scales_1.4.0                glue_1.8.0                  tools_4.5.1                
# [64] BiocIO_1.18.0               GenomicAlignments_1.44.0    mvtnorm_1.3-3              
# [67] XML_3.99-0.19               grid_4.5.1                  libcoin_1.0-12             
# [70] GenomeInfoDbData_1.2.14     restfulr_0.0.16             cli_3.6.5                  
# [73] rappdirs_0.3.3              S4Arrays_1.8.1              gtable_0.3.6               
# [76] digest_0.6.37               SparseArray_1.8.1           TH.data_1.1-4              
# [79] rjson_0.2.23                farver_2.1.2                memoise_2.0.1              
# [82] lifecycle_1.0.5             httr_1.4.7                  gridtext_0.1.5             
# [85] bit64_4.6.0-1               MASS_7.3-65   



