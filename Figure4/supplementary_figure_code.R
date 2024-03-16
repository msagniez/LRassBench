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

## extracting the beginning (of Sequins) predicted novel transcripts
new_allTools_beginPred <- new_allTools[which(new_allTools$end<97000),]

plot_transcripts_intron <- function(inputDF){
  exons_inputDF <- inputDF %>% dplyr::filter(type == "exon")
  exons_inputDF <- cbind("seqnames"=exons_inputDF$seqid,exons_inputDF)
  exons_inputDF <- cbind(exons_inputDF,"novel_transcr_name"=paste(exons_inputDF$source,exons_inputDF$transcript_id,sep="_"))
  level_order <- unique(exons_inputDF$novel_transcr_name)
  level_order <- level_order[str_order(level_order)]
  level_order <- level_order[c(which(grepl("Sequin_*",level_order)),which(!grepl("Sequin_*",level_order)))]
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
    labs(x="Genomic Coordinates", y="Transcript Name", fill="Tool") + theme_bw()
}

###print plot
plot_transcripts_intron(new_allTools_beginPred) 


