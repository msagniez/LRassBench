# ---- LIBRARIES ----
library(tidyverse)
library(cowplot)
library(readxl)
library(ggplot2)
library(ggrepel)

# ---- INPUT ----
# Subset Assembly TP FP FN Method Tool Evaluator
scores_all <- read_excel("Supplemental_Tables.xlsx", sheet = "SuppTable2")
scores_all <- scores_all %>% subset(select = -c(Method_old))

# Keep only the GFFcompare results
scores_150k <- scores_all %>% filter(Subset == "Sub150k")
scores_150k$Evaluator <- factor(scores_150k$Evaluator, levels = c("GFFcompare", "SQANTI3", "Average"))
scores_150k <- scores_150k %>%
  mutate(row_id = row_number()) %>%
  group_by(row_id) %>%
  do({
    if (.$Method[1] == "RawReads") {
      # Create three versions of this row
      out <- bind_rows(., ., .)
      out$Method <- c("AnnotationGuided", "GenomeGuided", "ReferenceFree")
      out
    } else {
      .
    }
  }) %>%
  ungroup() %>%
  select(-row_id)
scores_150k <- scores_150k %>% subset(Method != "RawReads")

plotdat <- scores_150k %>%
  select(Evaluator, Method, Dataset, Tool, Assembly,
         Precision, Sensitivity, F1) %>%
  pivot_longer(cols = c(Precision, Sensitivity, F1),
               names_to = "Metric",
               values_to = "value")


# ---- COLOURS ----
tool_cols <- c(
  "Bambu"       = "brown",
  "FLAIR"       = "violet",
  "isONclust"   = "mediumpurple1",
  "isONclust2"  = "purple",
  "isoQuant3"    = "turquoise2",
  "RATTLE"      = "springgreen2",
  "RNAbloom"    = "coral2",
  "RNAbloom2"   = "orangered3",
  "Stringtie3"  = "palegreen4",
  "TALON"       = "wheat2",
  "FLAMES"      = "blue2",
  "Mandalorion" = "darkgoldenrod1",
  "RawReads" = "lightgray",
  # new tools
  "Stringtie2"  = "darkseagreen2",
  "isoQuant"   = "lightskyblue1",
  "LRAA"        = "hotpink1",
  "LyRic"       = "firebrick1"
)
shape_map <- c(
  "RawReads_special" = 9,
  "AnnotationGuided" = 15,
  "GenomeGuided"     = 19,
  "ReferenceFree"    = 17
)

# ---- PLOTTING FUNCTION ----
make_panel <- function(dat, evaluator_label) {
  
  ggplot(
    dat,
    aes(
      x = Dataset, 
      y = value,
      colour = Tool,
      group = Assembly,
      linetype = Assembly == "TALON_reco",                     # dotted line logic
      shape = ifelse(Tool == "RawReads",                       # special RawReads shape
                     "RawReads_special",
                     Method)
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(
      Metric ~ Method,
      scales = "fixed"
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = tool_cols) +
    scale_shape_manual(values = shape_map) +
    scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "dotted")) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text.x = element_text(size = 11, face = "bold"),
      strip.text.y = element_text(size = 10, angle = 90),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    ggtitle(evaluator_label) +
    xlab("") + ylab("")
}

p_gff  <- make_panel(filter(plotdat, Evaluator == "GFFcompare"), "GFFcompare")
p_sq3  <- make_panel(filter(plotdat, Evaluator == "SQANTI3"),   "SQANTI3")
p_avg  <- make_panel(filter(plotdat, Evaluator == "Average"),   "Average")


final_fig <- plot_grid(p_gff, p_sq3, p_avg,
                       ncol = 1,
                       labels = NULL,
                       align = "v",
                       rel_heights = c(1,1,1))

# Save
ggsave("Sub150k_MultiEvaluator_MetricPanels.pdf",
       final_fig, width = 10, height = 16, dpi = 300)

final_fig


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




