#Volcano plot for 50K

library(tidyverse)
library(dplyr)

p1 <- read.csv("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_50000gen_400bp_foldchanges_mark1.csv")
p1$diffexpressed <- "not significant"
p1$diffexpressed[p1$log2FoldChange > 0 & p1$padj < 0.05] <- "upregulated"
p1$diffexpressed[p1$log2FoldChange < 0 & p1$padj < 0.05] <- "downregulated"
p1$log2FoldChange_rounded <- round(p1$log2FoldChange,3)
p_dedupe <- p1[!duplicated(p1$log2FoldChange_rounded),]

p_dedupe %>%
  #filter(rt_status!="readthrough") %>%
  filter(rt_status!="readthrough") %>%
  ggplot(aes(x=log2FoldChange_rounded, y=-log10(padj),col=diffexpressed,shape=type,size=type)) + 
  geom_point() +
  scale_shape_manual(values = c(1,19))+ 
  scale_size_manual(values = c(1,4))+ 
  theme_minimal() + 
  xlim(-15,15) + 
  ylim(0,120) + 
  scale_x_continuous(limits = c(-12,12),breaks=c(-16,-12,-8,-4,0,4,8,12,16)) +
  scale_color_manual(values=c("darkorange3","azure3","deepskyblue3")) + 
  #  geom_vline(xintercept=c(-2, 2), col="azure4") + 
  geom_hline(yintercept=-log10(0.05), col="azure4")+
  theme(#legend.key.size = unit(4, 'cm'), #change legend key size
    #legend.key.height = unit(0.75, 'cm'), #change legend key height
    #legend.key.width = unit(4, 'cm'), #change legend key width
    #legend.title = element_text(size=18), #change legend title font size
    #legend.text = element_text(size=16),
    axis.title = element_text(size=16),
    axis.text = element_text(size=18),
    legend.position="none") +
  labs(x = "Log2 fold change (from ancestor)")
  
#Time series bar graph
  
cat htseq_time_series_400bp_namesorted_foldchanges_mark1.csv | sed "s/gen5000,ZDB409/5000/g" | sed "s/gen10000,ZDB425/10000/g" | sed "s/gen15000,ZDB445/15000/g" | sed "s/gen20000,ZDB467/20000/g" | sed "s/gen25000,ZDB483/25000_1/g" | sed "s/gen25000,ZDB486/25000_2/g" | sed "s/gen25000,ZDB488/25000_3/g" | sed "s/gen27000,ZDB309/27000/g" | sed "s/gen30000,ZDB17/30000/g" | sed "s/gen31500,ZDB199/31500_1/g" | sed "s/gen31500,ZDB200/31500_2/g" | sed "s/gen31500,ZDB25/31500_3/g" | sed "s/gen31500,ZDB564/31500_4/g" | sed "s/gen33000,CZB154/33000/g" | sed "s/strain,//g" > testy && mv testy htseq_time_series_400bp_namesorted_foldchanges_mark1.csv

library(tidyverse)
library(dplyr)

p1 <- read.csv("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_time_series_400bp_namesorted_foldchanges_mark1.csv")
p1$diffexpressed <- "Not significant"
p1$diffexpressed[p1$log2FoldChange > 0 & p1$padj < 0.05] <- "upregulated"
p1$diffexpressed[p1$log2FoldChange < 0 & p1$padj < 0.05] <- "downregulated"
p1_nons <- p1 %>% filter(diffexpressed!="Not significant") %>% filter(type=="noncoding") %>% filter(rt_status!="readthrough")
p1_nons$generation <- factor(p1_nons$generation,levels = c("5000","10000","15000","20000","25000_1","25000_2","25000_3","27000","30000","31500_1","31500_2","31500_3","31500_4","33000"))
#p1_nons$diffexpressed <- factor(p1_nons$diffexpressed,levels = c("upregulated","downregulated","Not significant"))
ggplot(p1_nons, aes(x=generation,fill=diffexpressed)) +
  geom_bar(position="dodge", color = "azure4",width=0.75) +
  scale_fill_manual(values=c("darkorange","deepskyblue")) +
  ylim(0,50) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none") +
  xlab("Evolved strains (Ara-3)")+
  ylab("Non-coding windows with changed expression compared to ancestor")+
  ggsave("figure_1b_1.pdf",device="pdf",bg="transparent",width=12,height=6)
