#I extracted the generation - allele freq pairs manually with this code:

#For SNP or DEL:

grep "1881296" Ara-2*K/output/*gd | grep "SNP.*frequency=" | sed "s/K\//_/g" | sort -nk 2 -t "_" | cut -f 1,7 | sed "s/\t/_/g" | sed "s/=/_/g" | cut -f 2,5 -d "_" | sed "s/_/\t/g" | awk -F'\t' 'FNR == NR { data[$1] = $0; next } $1 in data { print data[$1] } !($1 in data) { print $1 "\t0" }' - 0_to_60

#For MOB:

grep "1881296" Ara-2*K/output/*gd | grep "MOB.*frequency=" | sed "s/K\//_/g" | sort -nk 2 -t "_" | cut -f 1,9 | sed "s/\t/_/g" | sed "s/=/_/g" | cut -f 2,5 -d "_" | sed "s/_/\t/g" | awk -F'\t' 'FNR == NR { data[$1] = $0; next } $1 in data { print data[$1] } !($1 in data) { print $1 "\t0" }' - 0_to_60

#Then pasted them in excel.

The R code (wip):

library(tidyverse)
library(ggplot2)

getwd()
setwd("/stor/scratch/Ochman/hassan/Good2017_LTEE_data/trimmed")

p1 <- read.csv("breseq.tsv",sep='\t')

p1 %>%
  #filter(target_id=="Ara-2_440_SNP") %>%
  #filter(target_id=="Ara-2_610_SNP") %>%
  #filter(target_id=="Ara-2_731_SNP") %>%
  #filter(target_id=="Ara-2_991_SNP") %>%
  #filter(target_id=="Ara-3_547_MOB") %>%
  #filter(target_id=="Ara-6_24_MOB") %>%
  #filter(target_id=="Ara-6_33_MOB") %>%
  #filter(target_id=="Ara-6_59_DEL") %>%
  #filter(target_id=="Ara+1_31_MOB") %>%
  #filter(target_id=="Ara+1_35_MOB") %>%
  #filter(target_id=="Ara+1_76_DEL") %>%
  #filter(target_id=="Ara+1_94_DEL") %>%
  #filter(target_id=="Ara+1_99_MOB") %>%
  #filter(target_id=="Ara+5_30_MOB") %>%
  #filter(target_id=="Ara+5_65_MOB") %>%
  #filter(target_id=="Ara+5_76_MOB") %>%
  ggplot(aes(x=Generation,y=Allele_freq)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 60000, 500)) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    #    panel.grid.major = element_blank(), #remove major gridlines
    #    panel.grid.minor = element_blank(), #remove minor gridlines
    #    panel.grid.major = element_line(size = 0.5, linetype = 'dashed',colour = "gray"),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size=10, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.key.size = unit(1, 'cm'),
    legend.text = element_text(size=15),
    legend.position="top")
