##loading libraries
library(rtracklayer)
library(tidyverse)
library(stringr)
library(ggtranscript)
library(magrittr)
library(dplyr)
library(ggplot2)
library("GenomicRanges")
library("GenomicFeatures")
library(readxl)
library(gridExtra)
library(ggtext)

##reading in all GTFs
Bambu <- readGFF("~/Documents/PhD_UDEM/work-year2/Benchmarking/GuidedAssemblies-Discovery/Bambu.gtf")
Flair <- readGFF("~/Documents/PhD_UDEM/work-year2/Benchmarking/GuidedAssemblies-Discovery/Flair.gtf")
isoQuant <- readGFF("~/Documents/PhD_UDEM/work-year2/Benchmarking/GuidedAssemblies-Discovery/isoQuant.gtf")
Stringtie2 <- readGFF("~/Documents/PhD_UDEM/work-year2/Benchmarking/GuidedAssemblies-Discovery/Stringtie2.gtf")
TALON_reco <- readGFF("~/Documents/PhD_UDEM/work-year2/Benchmarking/GuidedAssemblies-Discovery/TALON_reco.gtf")
TALON_reco$source<-gsub('TALON', 'TALONreco', TALON_reco$source)
TALON <- readGFF("~/Documents/PhD_UDEM/work-year2/Benchmarking/GuidedAssemblies-Discovery/TALON.gtf")
FLAMES <- as.data.frame(readGFF("~/Documents/PhD_UDEM/work-year2/Benchmarking/GuidedAssemblies-Discovery/FLAMES-filt.gff3"))
Mandalorion <- as.data.frame(readGFF("~/Documents/PhD_UDEM/work-year2/Benchmarking/GuidedAssemblies-Discovery/Mandalorion.gtf"))
OG_GTF <- readGFF("~/Documents/PhD_UDEM/work-year2/Benchmarking/rnasequin_annotation_2.4.gtf")
Disc_GTF <- readGFF("~/Documents/PhD_UDEM/work-year2/Benchmarking/rnasequin_annotation_2.4_TrcpDiscovery.gtf")

##modifying dataset's columns to suit transcript naming needs later
Flair$transcript_id <- gsub("\\|.*","",Flair$transcript_id)
FLAMES$geneID[FLAMES$geneID=="character(0)"] <- NA_character_
FLAMES$Parent[FLAMES$Parent=="character(0)"] <- NA_character_
FLAMES <- unite(FLAMES, trcp_id, c(ID,Parent), na.rm = TRUE, remove = F, sep = ";")
FLAMES$transcript_id <- gsub("^[a-z]*:","",FLAMES$trcp_id)
###subsetting original datasets to extract information on missing transcripts / artifically novel isoforms
novel_GTF <- OG_GTF[which(OG_GTF$transcript_id %in% c("R1_13_1", "R1_21_2", "R1_51_1", "R2_116_1", "R2_116_2", "R2_47_2", "R2_59_3", "R2_6_2", "R2_72_1")),]
altIso_GTF <- OG_GTF[which(OG_GTF$gene_id %in% c("R1_13", "R1_21", "R1_51", "R2_116", "R2_47", "R2_59", "R2_6", "R2_72")),] 
altIso_GTF <- altIso_GTF[which(!(altIso_GTF$transcript_id %in% c("R1_13_1", "R1_21_2", "R1_51_1", "R2_116_1", "R2_116_2", "R2_47_2", "R2_59_3", "R2_6_2", "R2_72_1"))),] ##known isoforms 
altIso_GTF_exons <- altIso_GTF %>% dplyr::filter(type == "exon") 
altIso_GTF_exons$start_end <- paste(altIso_GTF_exons$start, altIso_GTF_exons$end, sep = "_")
altIso_trcps <- na.omit(unique(altIso_GTF$transcript_id)) ##known isoform IDs
###filtering for exons common between novel & known
novel_GTF_start_end <- paste(novel_GTF$start, novel_GTF$end, sep = "_") ##vector of novel isoforms
altIso_GTF_exons <- altIso_GTF_exons %>% dplyr::filter(!(start_end %in% novel_GTF_start_end)) ##filtering novel exons out 

##colours for tools
ToolPalette<-c("Bambu" = "brown","FLAIR"="violet","isoQuant"="turquoise2","StringTie"="palegreen4","TALON"="wheat2","FLAMES"="blue2","TALONreco"="olivedrab3","Sequin"="mediumvioletred","Mandalorion"="darkgoldenrod1")

## function for matching the novel identified to 9 deleted transcripts. 
#### expected input: loaded GTF 
novel_matcher <- function(input_df){
  output_df <- data.frame()
  for(i in 1:nrow(input_df)){
    for (j in 1:nrow(novel_GTF)){
      if ((input_df$start[i]>= (novel_GTF$start[j]) ) && (input_df$end[i]<= (novel_GTF$end[j]) )){
        output_df <- rbind(output_df, input_df[i,])
      }
    }
  }
  #print(output_df)
  return(output_df)
}

expect_predic_all_Bambu <- novel_matcher(Bambu)
expect_predic_all_Flair <- novel_matcher(Flair)
expect_predic_all_isoQuant <- novel_matcher(isoQuant)
expect_predic_all_Stringtie2 <- novel_matcher(Stringtie2)
expect_predic_all_TALON <- novel_matcher(TALON)
expect_predic_all_TALONreco <- novel_matcher(TALON_reco)
expect_predic_all_FLAMES <- novel_matcher(FLAMES)
expect_predic_all_Mandalorion <- novel_matcher(Mandalorion)

expect_predic_all_TALON <- expect_predic_all_TALON[,-c(19)]  ##deleting extra 'source' column from TALON
expect_predic_all_TALONreco <- expect_predic_all_TALONreco[,-c(19)]  ##deleting extra 'source' column from TALONreco
expect_predic_all_TALON$source <- gsub("Sequin","TALON",expect_predic_all_TALON$source) ##changing source from Sequin to TALON to better rep in plots
expect_predic_all_TALONreco$source <- gsub("Sequin","TALONreco",expect_predic_all_TALONreco$source) ##changing source from Sequin to TALONreco to better rep in plots
expect_predic_all_isoQuant$source <- gsub("IsoQuant","isoQuant",expect_predic_all_isoQuant$source) ##changing source from IsoQuant to isoQuant to accurate rep in plots

expect_predic_all_ALL <- plyr::rbind.fill(expect_predic_all_Bambu, expect_predic_all_Flair, expect_predic_all_FLAMES, expect_predic_all_isoQuant, expect_predic_all_Mandalorion, expect_predic_all_Stringtie2, expect_predic_all_TALON, expect_predic_all_TALONreco)

##plotting funcs.
plot_transcripts_spl2 <- function(inputDF,tool_names){
  exons_inputDF <- inputDF %>% dplyr::filter(type == "exon")
  exons_inputDF <- cbind("seqnames"=exons_inputDF$seqid,exons_inputDF)
  exons_inputDF <- cbind(exons_inputDF,"novel_transcr_name"=paste(exons_inputDF$source,exons_inputDF$transcript_id,sep="_"))
  ##filter those exons that are perfect matches to alternate isoforms of same gene
  exons_inputDF$start_end <- paste(exons_inputDF$start, exons_inputDF$end, sep="_")
  exons_inputDF <- exons_inputDF %>% dplyr::filter(!(start_end %in% altIso_GTF_exons$start_end))
  ##filtering for "novel" isoforms
  for(tool_name in tool_names){
    eval_exp <- grepl(paste0("^",tool_name,".*R._(\\d+){1,2}_(\\d+)$"),exons_inputDF$novel_transcr_name)
    exons_inputDF <- exons_inputDF[which(!eval_exp),]  
  }
  ##ordering DF with Sequins on top
  level_order <- unique(exons_inputDF$novel_transcr_name)
  level_order <- level_order[str_order(level_order)]
  level_order <- level_order[c(which(grepl("Sequin_*",level_order)),which(!grepl("Sequin_*",level_order)))]
  ##chunk for colouring y axis labels
  LabelPalette <- c()
  match_str_all <- unique(gsub("_.*","",level_order)) ##level order is also novel_transcr_name
  for(i in 1:length(match_str_all)){
    num_matches <- length(which(grepl(paste0(match_str_all[i],"_"),level_order)))
    LabelPalette <- c(LabelPalette,rep(ToolPalette[match_str_all[i]],num_matches))
  }
  names(LabelPalette) <- NULL
  LabelPalette <- as.character(rev(LabelPalette))
  ##chunk for colouring y axis labels ends
  exons_inputDF %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = factor(novel_transcr_name, level = level_order)
    )) +
    geom_range(
      aes(fill = source)
    ) + 
    scale_fill_manual(values=ToolPalette) +
    scale_y_discrete(limits=rev) +
    labs(x="Genomic Coordinates", y="Transcript Name", fill="Tool") + theme_bw() +
    theme(axis.text.y = element_text(colour = LabelPalette))
  
}

plot_transcripts_intron2 <- function(inputDF,tool_names){
  exons_inputDF <- inputDF %>% dplyr::filter(type == "exon")
  exons_inputDF <- cbind("seqnames"=exons_inputDF$seqid,exons_inputDF)
  exons_inputDF <- cbind(exons_inputDF,"novel_transcr_name"=paste(exons_inputDF$source,exons_inputDF$transcript_id,sep="_"))
  for(tool_name in tool_names){
    eval_exp <- grepl(paste0("^",tool_name,".*R._(\\d+){1,2}_(\\d+)$"),exons_inputDF$novel_transcr_name)
    exons_inputDF <- exons_inputDF[which(!eval_exp),]  
  }
  level_order <- unique(exons_inputDF$novel_transcr_name)
  level_order <- level_order[str_order(level_order)]
  level_order <- level_order[c(which(grepl("Sequin_*",level_order)),which(!grepl("Sequin_*",level_order)))]
  ##chunk for colouring y axis labels
  LabelPalette <- c()
  match_str_all <- unique(gsub("_.*","",level_order)) ##level order is also
  for(i in 1:length(match_str_all)){
    num_matches <- length(which(grepl(paste0(match_str_all[i],"_"),level_order)))
    LabelPalette <- c(LabelPalette,rep(ToolPalette[match_str_all[i]],num_matches))
  }
  names(LabelPalette) <- NULL
  LabelPalette <- as.character(rev(LabelPalette))
  ##chunk for colouring y axis labels ends
  exons_inputDF %>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = factor(novel_transcr_name, level = level_order)
    )) +
    geom_range(
      aes(fill = source)
    ) + 
    geom_intron(
      data = to_intron(exons_inputDF, "novel_transcr_name"),
      arrow = grid::arrow(ends = "last", length = grid::unit(1, "mm")),
      arrow.min.intron.length = 10,
      aes(strand = strand)
    ) + 
    scale_fill_manual(values=ToolPalette) +
    scale_y_discrete(limits=rev) +
    labs(x=" ", y=" ", fill="Tool") + theme_bw() +
    theme(axis.text=element_text(size=6), axis.title=element_text(size=8), legend.position="none", 
          axis.text.y = element_text(colour = LabelPalette)) ##also for coloring y axis label
}

gene_region_selector <- function(geneID,input_df){
  output_df <- data.frame()
  input_df$start_end <- paste(input_df$start, input_df$end, sep="_")
  for(i in 1:nrow(input_df)){
    for (j in which(novel_GTF$gene_id==geneID)){
      ##select all exons matching to exons of missing isoforms
      if ((input_df$start[i]>= novel_GTF$start[j]) && (input_df$end[i]<= (novel_GTF$end[j]) )){ 
        ##it belongs in the region of desired gene
        output_df <- rbind(output_df, input_df[i,])
        ##filter those exons that are perfect matches to alternate isoforms of same gene
        output_df <- output_df %>% dplyr::filter(!(start_end %in% altIso_GTF_exons$start_end))
        
        
      }
    }
  }
  return(output_df)
}

gridPlotter_isoforms <- function(tool_name){
  nam <- paste0("expect_predic_all_",tool_name)
  tool_names <- unique(get(nam)$source)
  p1<-plot_transcripts_spl2(plyr::rbind.fill(get(nam),novel_GTF), tool_names) 
  p2<-plot_transcripts_intron2(plyr::rbind.fill(gene_region_selector("R1_13",get(nam)), OG_GTF[OG_GTF$transcript_id=="R1_13_1",]), tool_names) 
  p3<-plot_transcripts_intron2(plyr::rbind.fill(gene_region_selector("R2_59",get(nam)), OG_GTF[OG_GTF$transcript_id=="R2_59_3",]), tool_names) 
  p4<-plot_transcripts_intron2(plyr::rbind.fill(gene_region_selector("R1_21",get(nam)), OG_GTF[OG_GTF$transcript_id=="R1_21_2",]), tool_names) 
  p5<-plot_transcripts_intron2(plyr::rbind.fill(gene_region_selector("R2_6",get(nam)), OG_GTF[OG_GTF$transcript_id=="R2_6_2",]), tool_names) 
  p6<-plot_transcripts_intron2(plyr::rbind.fill(gene_region_selector("R2_116",get(nam)), OG_GTF[OG_GTF$transcript_id=="R2_116_1",]), tool_names) 
  p7<-plot_transcripts_intron2(plyr::rbind.fill(gene_region_selector("R2_47",get(nam)), OG_GTF[OG_GTF$transcript_id=="R2_47_2",]), tool_names) 
  p8<-plot_transcripts_intron2(plyr::rbind.fill(gene_region_selector("R2_72",get(nam)), OG_GTF[OG_GTF$transcript_id=="R2_72_1",]), tool_names) 
  p9<-plot_transcripts_intron2(plyr::rbind.fill(gene_region_selector("R2_116",get(nam)), OG_GTF[OG_GTF$transcript_id=="R2_116_2",]), tool_names)
  p10<-plot_transcripts_intron2(plyr::rbind.fill(gene_region_selector("R1_51",get(nam)), OG_GTF[OG_GTF$transcript_id=="R1_51_1",]), tool_names)
  gridExtra::grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow = 4, widths = c(1, 1, 1), heights = c(4,1,1,1),
                          layout_matrix = rbind(c(1),
                                                c(2, 3, 4),
                                                c(5, 6, 7),
                                                c(8, 9, 10)))
}

gridPlotter_isoforms("ALL")

##could be used for individual tools as well:
# gridPlotter_isoforms("Bambu")
