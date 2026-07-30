##Figure 2C

# ---- LIBRARIES ----
library(tidyverse)
library(cowplot)
library(readxl)
library(ggplot2)
library(ggrepel)


# ---- INPUT ----
scores_all <- read_excel("Supplemental_Tables.xlsx", sheet = "SuppTable2")
scores_all <- scores_all %>% subset(select = -c(Method_old))

#head(scores_all)

LSK114_sub150k <- subset(scores_all, Subset == "Sub150k" & Dataset == "LSK114_sequins")
PCS109_sub150k <- subset(scores_all, Subset == "Sub150k" & Dataset == "LSK109_sequins")
RNA002_sub150k <- subset(scores_all, Subset == "Sub150k" & Dataset == "RNA002_sequins")
RNA004_sub150k <- subset(scores_all, Subset == "Sub150k" & Dataset == "RNA004_sequins")

#sub40k
chrIS_sub40k <- subset(scores_all, Subset == "Sub40k" & Dataset == "LSK114_sequins")
SIRV_sub40k <- subset(scores_all, Subset == "Sub40k" & Dataset == "LSK114_SIRVs")


#Plot scores_all
tool_cols <- c(
  "Bambu" = "brown",
  "FLAIR" = "violet",
  "isONclust" = "mediumpurple1",
  "isONclust2" = "purple",
  "isoQuant3" = "turquoise2",
  "RATTLE" = "springgreen2",
  "RNAbloom" = "coral2",
  "RNAbloom2" = "orangered3",
  "Stringtie3" = "palegreen4",
  "TALON" = "wheat2",
  "FLAMES" = "blue2",
  "Mandalorion" = "darkgoldenrod1",
  "Stringtie2"  = "darkseagreen2",
  "isoQuant"   = "lightskyblue1",
  "LRAA"        = "hotpink1",
  "LyRic"       = "firebrick1"
)
shape_vec1 <- c(15, 19, 9, 17)

## Supplementary figures ; Precision v/s sensitivity 
PCS109_sub150k_GFF <- ggplot(
  subset(PCS109_sub150k, Evaluator == "GFFcompare"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("PCS109_sequins") +
  annotate( "text", x = 0, y = 0, label = "GFFcompare", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

LSK114_sub150k_GFF <- ggplot(
  subset(LSK114_sub150k, Evaluator == "GFFcompare"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("LSK114_sequins") +
  annotate( "text", x = 0, y = 0, label = "GFFcompare", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

RNA002_sub150k_GFF <- ggplot(
  subset(RNA002_sub150k, Evaluator == "GFFcompare"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("RNA002_sequins") +
  annotate( "text", x = 0, y = 0, label = "GFFcompare", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

RNA004_sub150k_GFF <- ggplot(
  subset(RNA004_sub150k, Evaluator == "GFFcompare"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("RNA004_sequins") +
  annotate( "text", x = 0, y = 0, label = "GFFcompare", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")



PCS109_sub150k_SQ3 <- ggplot(
  subset(PCS109_sub150k, Evaluator == "SQANTI3"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("PCS109_sequins") +
  annotate( "text", x = 0, y = 0, label = "SQANTI3", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

LSK114_sub150k_SQ3 <- ggplot(
  subset(LSK114_sub150k, Evaluator == "SQANTI3"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("LSK114_sequins") +
  annotate( "text", x = 0, y = 0, label = "SQANTI3", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

RNA002_sub150k_SQ3 <- ggplot(
  subset(RNA002_sub150k, Evaluator == "SQANTI3"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("RNA002_sequins") +
  annotate( "text", x = 0, y = 0, label = "SQANTI3", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

RNA004_sub150k_SQ3 <- ggplot(
  subset(RNA004_sub150k, Evaluator == "SQANTI3"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) +
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("RNA004_sequins") +
  annotate( "text", x = 0, y = 0, label = "SQANTI3", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")




#Plot sub40k
chrIS_sub40k_GFF <- ggplot(
  subset(chrIS_sub40k, Evaluator == "GFFcompare"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("LSK114_sequins") +
  annotate( "text", x = 0, y = 0, label = "GFFcompare", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

chrIS_sub40k_SQ3 <- ggplot(
  subset(chrIS_sub40k, Evaluator == "SQANTI3"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("LSK114_sequins") +
  annotate( "text", x = 0, y = 0, label = "SQANTI3", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

SIRV_sub40k_GFF <- ggplot(
  subset(SIRV_sub40k, Evaluator == "GFFcompare"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("LSK114_SIRV") +
  annotate( "text", x = 0, y = 0, label = "GFFcompare", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

SIRV_sub40k_SQ3 <- ggplot(
  subset(SIRV_sub40k, Evaluator == "SQANTI3"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("LSK114_SIRV") +
  annotate( "text", x = 0, y = 0, label = "SQANTI3", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

plot_withLegend <- ggplot(
  subset(PCS109_sub150k, Evaluator == "GFFcompare"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1)

legend <- get_legend(plot_withLegend)


PCS109_sub150k_avg <- ggplot(
  subset(PCS109_sub150k, Evaluator == "Average"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("PCS109_sequins") +
  annotate( "text", x = 0, y = 0, label = "Average", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

LSK114_sub150k_avg <- ggplot(
  subset(LSK114_sub150k, Evaluator == "Average"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("LSK114_sequins") +
  annotate( "text", x = 0, y = 0, label = "Average", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

RNA002_sub150k_avg <- ggplot(
  subset(RNA002_sub150k, Evaluator == "Average"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("RNA002_sequins") +
  annotate( "text", x = 0, y = 0, label = "Average", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

RNA004_sub150k_avg <- ggplot(
  subset(RNA004_sub150k, Evaluator == "Average"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("RNA004_sequins") +
  annotate( "text", x = 0, y = 0, label = "Average", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

chrIS_sub40k_avg <- ggplot(
  subset(chrIS_sub40k, Evaluator == "Average"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("LSK114_sequins") +
  annotate( "text", x = 0, y = 0, label = "Average", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")

SIRV_sub40k_avg <- ggplot(
  subset(SIRV_sub40k, Evaluator == "Average"),
  aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
  geom_point(aes(colour = Tool, shape = Method), size = 4) + 
  geom_text_repel(size = 2, show.legend = FALSE) +
  xlim(0, 1) + ylim(0, 1) + 
  theme_bw() + 
  scale_color_manual(values = tool_cols) + 
  scale_shape_manual(values = shape_vec1) + 
  theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
  ggtitle("LSK114_SIRV") +
  annotate( "text", x = 0, y = 0, label = "Average", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")


#Pick main figure from following
SuppFig2_pt1 <- plot_grid(chrIS_sub40k_GFF, SIRV_sub40k_GFF, chrIS_sub40k_SQ3, SIRV_sub40k_SQ3, chrIS_sub40k_avg, SIRV_sub40k_avg, nrow = 3, rel_widths=c(6,6,6))
SuppFig2_pt2 <- plot_grid(PCS109_sub150k_GFF,LSK114_sub150k_GFF,RNA002_sub150k_GFF,RNA004_sub150k_GFF, PCS109_sub150k_SQ3,LSK114_sub150k_SQ3,RNA002_sub150k_SQ3,RNA004_sub150k_SQ3, PCS109_sub150k_avg,LSK114_sub150k_avg,RNA002_sub150k_avg,RNA004_sub150k_avg, nrow = 3, rel_widths=c(6,6,6,6,6))
ggsave("SuppFig_PrecSensitivity_EvalGFFSQ3avg_sub150k_sub40k.pdf",
       plot_grid(SuppFig2_pt1, SuppFig2_pt2,legend, nrow=1,rel_widths=c(8,16,2)), 
       width = 24, height = 11, dpi = 300)
## pdf 


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
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#  [1] ggrepel_0.9.6   readxl_1.4.5    cowplot_1.2.0   lubridate_1.9.4 forcats_1.0.0  
#  [6] stringr_1.5.2   dplyr_1.2.0     purrr_1.1.0     readr_2.1.5     tidyr_1.3.1    
# [11] tibble_3.3.0    ggplot2_4.0.2   tidyverse_2.0.0

# loaded via a namespace (and not attached):
#  [1] sandwich_3.1-1     generics_0.1.4     coin_1.4-3         stringi_1.8.7     
#  [5] lattice_0.22-7     hms_1.1.3          magrittr_2.0.3     grid_4.5.1        
#  [9] timechange_0.3.0   RColorBrewer_1.1-3 mvtnorm_1.3-3      cellranger_1.1.0  
# [13] Matrix_1.7-3       survival_3.8-3     multcomp_1.4-28    scales_1.4.0      
# [17] TH.data_1.1-4      codetools_0.2-20   modeltools_0.2-24  cli_3.6.5         
# [21] rlang_1.1.7        splines_4.5.1      withr_3.0.2        tools_4.5.1       
# [25] parallel_4.5.1     tzdb_0.5.0         vctrs_0.7.1        R6_2.6.1          
# [29] matrixStats_1.5.0  stats4_4.5.1       zoo_1.8-14         lifecycle_1.0.5   
# [33] libcoin_1.0-12     MASS_7.3-65        pkgconfig_2.0.3    pillar_1.11.0     
# [37] gtable_0.3.6       glue_1.8.0         Rcpp_1.1.0         tidyselect_1.2.1  
# [41] rstudioapi_0.17.1  farver_2.1.2       compiler_4.5.1     S7_0.2.0     