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
library(ggrepel)


# ---- INPUT ----
scores_depth <- as.data.frame(read_excel("Supplemental_Tables.xlsx", sheet = "SuppTable6"))
scores_depth <- scores_depth %>% subset(select = -c(Method_old))
scores_depth$Subset <- gsub("Sub","",scores_depth$Subset)
GFF <- scores_depth %>% subset(Evaluator=="GFF")
GFF <- GFF %>% dplyr::select(-Evaluator, -Dataset)
SQ3 <- scores_depth %>% subset(Evaluator=="SQ3")
SQ3 <- SQ3 %>% dplyr::select(-Evaluator, -Dataset)

## averaging TP, novel, FP and FN values obtained from GFFcompare and SQANTI3
df2 <- GFF
df2$key <- paste0(df2$Subset,"_",df2$Assembly,"_",df2$Method)
df3 <- SQ3
df3$key <- paste0(df3$Subset,"_",df3$Assembly,"_",df3$Method)
df3$key <- gsub("_AWKed","",df3$key)

cols_to_avg <- c("TP","FP","FN","Precision","Sensitivity","F1")
df <- df2 %>%
  dplyr::inner_join(df3, by = "key", suffix = c(".df2", ".df3")) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(paste0(cols_to_avg, ".df2")),
      ~ (. + get(sub("\\.df2$", ".df3", dplyr::cur_column()))) / 2
    )
  ) %>%
  # keep key + all df2 columns (with selected ones now averaged)
  dplyr::select(key, dplyr::ends_with(".df2")) %>%
  dplyr::rename_with(~ sub("\\.df2$", "", .))
df$Evaluator = "Average"
df <- df[,-c(1)]


## initialising the colours for ggplot2
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
  "LyRic"       = "firebrick1",
  "RawReads" = "grey40"
)
shape_vec1 <- c("AnnotationGuided"=15, "GenomeGuided"=19, "RawReads"=9, "ReferenceFree"=17)


# LSK114_depth_avg <- ggplot(df,
#   aes(fill = Tool, x = Precision, y = Sensitivity, label = Assembly) ) + 
#   geom_point(aes(colour = Tool, shape = Method), size = 4) + 
#   geom_text_repel(size = 2, show.legend = FALSE) +
#   xlim(0, 1) + ylim(0, 1) + 
#   theme_bw() + 
#   scale_color_manual(values = tool_cols) + 
#   scale_shape_manual(values = shape_vec1) + 
#   theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5, face = "bold") ) +
#   ggtitle("LSK114_sequins") +
#   annotate( "text", x = 0, y = 0, label = "Average", hjust = -0.1, vjust = -1.2, size = 4, fontface = "bold")
# 
AVGTPFP <- df
AVGTPFP$Tool   <- as.character(AVGTPFP$Tool)
AVGTPFP$Method <- as.factor(AVGTPFP$Method)

AVGGuided <- subset(AVGTPFP, Method=="AnnotationGuided" | Method=="RawReads")
AVGDeNovo <- subset(AVGTPFP, Method=="GenomeGuided" | Method=="RawReads")
AVGAbInitio <- subset(AVGTPFP, Method=="ReferenceFree" | Method=="RawReads")
AVGGd <- ggplot(AVGGuided,
                aes(
                  fill = Tool,
                  x = Precision,
                  y = Sensitivity,
                  label = Subset
                )) + geom_path(aes(group=Tool, colour = Tool), linewidth = 0.8, alpha = 0.8) + geom_point(aes(colour = Tool, shape = Method))  + geom_text_repel(aes(label = Subset), size = 3, show.legend = FALSE) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw() + scale_color_manual(values = tool_cols) + scale_shape_manual(values = c(15, 9))

AVGDN <- ggplot(AVGDeNovo,
                aes(
                  fill = Tool,
                  x = Precision,
                  y = Sensitivity,
                  label = Subset
                )) + geom_path(aes(group=Tool, colour = Tool), linewidth = 0.8, alpha = 0.8) + geom_point(aes(colour = Tool, shape = Method))  + geom_text_repel(aes(label = Subset), size = 3, show.legend = FALSE) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw() + scale_color_manual(values = tool_cols) + scale_shape_manual(values = c(19, 9))

AVGAI <- ggplot(AVGAbInitio,
                aes(
                  fill = Tool,
                  x = Precision,
                  y = Sensitivity,
                  label = Subset
                )) + geom_path(aes(group=Tool, colour = Tool), linewidth = 0.8, alpha = 0.8) + geom_point(aes(colour = Tool, shape = Method))  + geom_text_repel(aes(label = Subset), size = 3, show.legend = FALSE) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw() + scale_color_manual(values = tool_cols) + scale_shape_manual(values = c(9, 17))

AVGGd <- AVGGd + theme(legend.position = "none")
AVGDN <- AVGDN + theme(legend.position = "none")
AVGAI <- AVGAI + theme(legend.position = "none")

plot_grid(AVGGd,AVGDN,AVGAI,nrow=1,rel_widths = c(6,6,6))

#Supp Fig
ggsave("SupplementaryB_Depth_Avg_noLegend.pdf_noLegend.pdf",
       plot_grid(AVGGd,AVGDN,AVGAI,nrow=1,rel_widths = c(6,6,6)), 
       width = 16, height = 6, dpi = 300)
## pdf 

### per Evaluator

df <- GFF
df$Evaluator = "GFFcompare"
AVGTPFP <- df
AVGTPFP$Tool   <- as.character(AVGTPFP$Tool)
AVGTPFP$Method <- as.factor(AVGTPFP$Method)

AVGGuided <- subset(AVGTPFP, Method=="AnnotationGuided" | Method=="RawReads")
AVGDeNovo <- subset(AVGTPFP, Method=="GenomeGuided" | Method=="RawReads")
AVGAbInitio <- subset(AVGTPFP, Method=="ReferenceFree" | Method=="RawReads")
AVGGd <- ggplot(AVGGuided,
                aes(
                  fill = Tool,
                  x = Precision,
                  y = Sensitivity,
                  label = Subset
                )) + geom_path(aes(group=Tool, colour = Tool), linewidth = 0.8, alpha = 0.8) + geom_point(aes(colour = Tool, shape = Method))  + geom_text_repel(aes(label = Subset), size = 3, show.legend = FALSE) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw() + scale_color_manual(values = tool_cols) + scale_shape_manual(values = c(15, 9))

AVGDN <- ggplot(AVGDeNovo,
                aes(
                  fill = Tool,
                  x = Precision,
                  y = Sensitivity,
                  label = Subset
                )) + geom_path(aes(group=Tool, colour = Tool), linewidth = 0.8, alpha = 0.8) + geom_point(aes(colour = Tool, shape = Method))  + geom_text_repel(aes(label = Subset), size = 3, show.legend = FALSE) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw() + scale_color_manual(values = tool_cols) + scale_shape_manual(values = c(19, 9))

AVGAI <- ggplot(AVGAbInitio,
                aes(
                  fill = Tool,
                  x = Precision,
                  y = Sensitivity,
                  label = Subset
                )) + geom_path(aes(group=Tool, colour = Tool), linewidth = 0.8, alpha = 0.8) + geom_point(aes(colour = Tool, shape = Method))  + geom_text_repel(aes(label = Subset), size = 3, show.legend = FALSE) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw() + scale_color_manual(values = tool_cols) + scale_shape_manual(values = c(9, 17))

AVGGd <- AVGGd + theme(legend.position = "none")
AVGDN <- AVGDN + theme(legend.position = "none")
AVGAI <- AVGAI + theme(legend.position = "none")

plot_grid(AVGGd,AVGDN,AVGAI,nrow=1,rel_widths = c(6,6,6))

#Main figure

ggsave("Figure4_noLegend.pdf",
       plot_grid(AVGGd,AVGDN,AVGAI,nrow=1,rel_widths = c(6,6,6)), 
       width = 16, height = 6, dpi = 300)



df <- SQ3
df$Evaluator = "SQANTI3"
AVGTPFP <- df
AVGTPFP$Tool   <- as.character(AVGTPFP$Tool)
AVGTPFP$Method <- as.factor(AVGTPFP$Method)

AVGGuided <- subset(AVGTPFP, Method=="AnnotationGuided" | Method=="RawReads")
AVGDeNovo <- subset(AVGTPFP, Method=="GenomeGuided" | Method=="RawReads")
AVGAbInitio <- subset(AVGTPFP, Method=="ReferenceFree" | Method=="RawReads")
AVGGd <- ggplot(AVGGuided,
                aes(
                  fill = Tool,
                  x = Precision,
                  y = Sensitivity,
                  label = Subset
                )) + geom_path(aes(group=Tool, colour = Tool), linewidth = 0.8, alpha = 0.8) + geom_point(aes(colour = Tool, shape = Method))  + geom_text_repel(aes(label = Subset), size = 3, show.legend = FALSE) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw() + scale_color_manual(values = tool_cols) + scale_shape_manual(values = c(15, 9))

AVGDN <- ggplot(AVGDeNovo,
                aes(
                  fill = Tool,
                  x = Precision,
                  y = Sensitivity,
                  label = Subset
                )) + geom_path(aes(group=Tool, colour = Tool), linewidth = 0.8, alpha = 0.8) + geom_point(aes(colour = Tool, shape = Method))  + geom_text_repel(aes(label = Subset), size = 3, show.legend = FALSE) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw() + scale_color_manual(values = tool_cols) + scale_shape_manual(values = c(19, 9))

AVGAI <- ggplot(AVGAbInitio,
                aes(
                  fill = Tool,
                  x = Precision,
                  y = Sensitivity,
                  label = Subset
                )) + geom_path(aes(group=Tool, colour = Tool), linewidth = 0.8, alpha = 0.8) + geom_point(aes(colour = Tool, shape = Method))  + geom_text_repel(aes(label = Subset), size = 3, show.legend = FALSE) + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_bw() + scale_color_manual(values = tool_cols) + scale_shape_manual(values = c(9, 17))

AVGGd <- AVGGd + theme(legend.position = "none")
AVGDN <- AVGDN + theme(legend.position = "none")
AVGAI <- AVGAI + theme(legend.position = "none")

plot_grid(AVGGd,AVGDN,AVGAI,nrow=1,rel_widths = c(6,6,6))

#Supp figure

ggsave("SupplementaryB_Depth_SQANTI3_noLegend.pdf",
       plot_grid(AVGGd,AVGDN,AVGAI,nrow=1,rel_widths = c(6,6,6)), 
       width = 16, height = 6, dpi = 300)

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

