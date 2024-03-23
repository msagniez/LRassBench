## Plot figure 6 ##

setwd("D:/Assembly-Benchmarking/Fig6-Time")
library(ggplot2)
library(cowplot)

Tfig1 <- read.csv("Rfig-Guided.csv",sep=";",header=T)
Tfig2 <- read.csv("Rfig-DeNovo.csv",sep=";",header=T)
Tfig3 <- read.csv("Rfig-AbInitio.csv",sep=";",header=T)
Tfig1$Step <- factor(Tfig1$Step, levels=c("Alignment","Pre-process","Assembler","Post-process"))
Tfig2$Step <- factor(Tfig2$Step, levels=c("Alignment","Pre-process","Assembler","Post-process"))
Tfig3$Step <- factor(Tfig3$Step, levels=c("Alignment","Pre-process","Assembler","Post-process"))

F1 <- ggplot(Tfig1,aes(fill=Step,y=Assembly,x=Value)) + geom_bar(position=position_stack(reverse = TRUE), stat="identity") + scale_fill_manual(values = c("Alignment"="mediumvioletred","Pre-process"="midnightblue","Assembler"="mediumturquoise","Post-process"="mediumspringgreen")) + theme_bw() + xlab("") + ylab("Guided") + scale_x_continuous(limits = c(0,60), expand = c(0, 0))
F2 <- ggplot(Tfig2,aes(fill=Step,y=Assembly,x=Value)) + geom_bar(position=position_stack(reverse = TRUE), stat="identity") + scale_fill_manual(values = c("Alignment"="mediumvioletred","Pre-process"="midnightblue","Assembler"="mediumturquoise","Post-process"="mediumspringgreen")) + theme_bw() + xlab("") + ylab("De novo") + scale_x_continuous(limits = c(0,60), expand = c(0, 0))
F3 <- ggplot(Tfig3,aes(fill=Step,y=Assembly,x=Value)) + geom_bar(position=position_stack(reverse = TRUE), stat="identity") + scale_fill_manual(values = c("Alignment"="mediumvioletred","Pre-process"="midnightblue","Assembler"="mediumturquoise","Post-process"="mediumspringgreen")) + theme_bw() + xlab("Time (min)") + ylab("Ab initio") + scale_x_continuous(limits = c(0,60), expand = c(0, 0))

plot_grid(F1,F2,F3,ncol=1,rel_heights=c(7.5,4.5,5))
#Pb with De novo plot which is wider than the other 2 because of Assembly labels

F1 <- F1 + theme(axis.text.y = element_blank())
F2 <- F2 + theme(axis.text.y = element_blank())
F3 <- F3 + theme(axis.text.y = element_blank())
plot_grid(F1,F2,F3,ncol=1,rel_heights=c(7.5,4.5,5))

#pdf 10x5 with and without labels
