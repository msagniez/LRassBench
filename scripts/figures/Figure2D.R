# ---- LIBRARIES ----
library(tidyverse)
library(cowplot)
library(readxl)
library(ggplot2)
library(ggrepel)

# ---- INPUT ----
scores_all <- read_excel("Supplemental_Tables.xlsx", sheet = "SuppTable2")

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




# --- helper: reshape to wide ---
to_wide <- function(df) {
  df %>%
    select(Assembly, Tool, Method, Dataset, Evaluator,
           Precision, Sensitivity, F1) %>%
    pivot_wider(
      id_cols = c(Assembly, Tool, Method, Dataset),
      names_from = Evaluator, 
      values_from = c(Precision, Sensitivity, F1),
      names_sep = "_",
      values_fn = mean   # collapse multiple entries to a single numeric
    )
}
to_wide_v2 <- function(df) {
  df %>%
    select(Assembly, Tool, Method, Dataset, Subset, Evaluator,
           Precision, Sensitivity, F1) %>%
    pivot_wider(
      id_cols = c(Assembly, Tool, Method, Dataset, Subset),
      names_from = Evaluator,
      values_from = c(Precision, Sensitivity, F1),
      names_sep = "_",
      values_fn = mean   # collapse multiple entries to a single numeric
    )
}

# --- helper: one metric scatter ---
metric_plot <- function(dfw, metric) {
  
  xcol <- paste0(metric, "_GFFcompare")
  ycol <- paste0(metric, "_SQANTI3")
  
  # compute absolute deviation from x=y
  dfw <- dfw %>%
    mutate(
      deviation = abs(.data[[ycol]] - .data[[xcol]]),
      label_me  = deviation > 0.10   # label only if >10% deviation
    )
  
  ggplot(dfw, aes(.data[[xcol]], .data[[ycol]],
                  colour = Tool, label = Dataset)) +
    geom_abline(slope = 1, linetype = "dashed",
                color = "goldenrod2", linewidth = 1) +
    geom_point(aes(shape = Method), size = 2.5) +
    
    # --- label only those above threshold ---
    ggrepel::geom_text_repel(
      data = subset(dfw, label_me),
      size = 2.2,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = tool_cols) +
    scale_shape_manual(values = shape_vec1) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    theme_bw() +
    theme(legend.position="none",
          plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle(metric)
}
metric_plot_v2 <- function(dfw, metric) {
  
  xcol <- paste0(metric, "_GFFcompare")
  ycol <- paste0(metric, "_SQANTI3")
  
  # compute absolute deviation from x=y
  dfw <- dfw %>%
    mutate(
      deviation = abs(.data[[ycol]] - .data[[xcol]]),
      label_me  = deviation > 0.10   # label only if >10% deviation
    )
  
  ggplot(dfw, aes(.data[[xcol]], .data[[ycol]],
                  colour = Tool, label = Subset)) +
    geom_abline(slope = 1, linetype = "dashed",
                color = "goldenrod2", linewidth = 1) +
    geom_point(aes(shape = Method), size = 2.5) +
    
    # --- label only those above threshold ---
    ggrepel::geom_text_repel(
      data = subset(dfw, label_me),
      size = 2.2,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = tool_cols) +
    scale_shape_manual(values = shape_vec1) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    theme_bw() +
    theme(legend.position="none",
          plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle(metric)
}

# --- helper: three stacked plots for one subset ---
panel_3metrics <- function(df, title) {
  dfw <- to_wide(df)
  p1 <- metric_plot(dfw, "Precision")
  p2 <- metric_plot(dfw, "Sensitivity")
  p3 <- metric_plot(dfw, "F1")
  
  plot_grid(
    ggdraw() + draw_label(title, fontface="bold", size=14, hjust=0.5),
    cowplot::plot_grid(p1, p2, p3, ncol=1, rel_heights=c(1,1,1)),
    ncol=1,
    rel_heights=c(0.10, 1)
  )
}

# --- build 3×3 grid ---
df_40k   <- scores_all %>% filter(Subset == "Sub40k")
p40   <- panel_3metrics(df_40k,  "Subset = 40k reads")
df_150k  <- scores_all %>% filter(Subset == "Sub150k")
p150  <- panel_3metrics(df_150k, "Subset = 150k reads")

# --- helper: three stacked plots for one subset ---
panel_3metrics <- function(df, title) {
  dfw <- to_wide_v2(df)
  p1 <- metric_plot_v2(dfw, "Precision")
  p2 <- metric_plot_v2(dfw, "Sensitivity")
  p3 <- metric_plot_v2(dfw, "F1")
  
  plot_grid(
    ggdraw() + draw_label(title, fontface="bold", size=14, hjust=0.5),
    cowplot::plot_grid(p1, p2, p3, ncol=1, rel_heights=c(1,1,1)),
    ncol=1,
    rel_heights=c(0.10, 1)
  )
}

df_all   <- scores_all %>% filter(Dataset == "LSK114_sequins")
pALL  <- panel_3metrics(df_all,  "All Subsets - LSK114_sequins")

final_grid <- plot_grid(p40, p150, pALL, ncol=3)

# save
ggsave("EvaluatorComparison_labelDataset_grid.pdf",
       final_grid, width=10, height=10, dpi=300)


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


