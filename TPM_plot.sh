#for 100bp

cat *ancestor*htseqready.gtf | cut -f 4-7,9 | sort -u | cut -f 1 -d ';' | sed "s/transcript_id \"//g" | sed "s/\"//g" | awk -F '\t' '{OFS=FS}{print $5,$2-$1}' | sort -k1 > 100bp_lengths
for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do
sort -k 4 -o "$i"_500_htseq.tsv "$i"_500_htseq.tsv
join -1 1 -2 4 100bp_lengths "$i"_500_htseq.tsv | sed "s/ /\t/g" | sed -e '/REL60/s/$/\tancestor/' -e '/REL60/!s/$/\tevolved/' | sed "1s/^/target_id\tlength\treplicate\tseqtype\tline\tcount\tstatus\n/g" > "$i"_htseq_lengths.tsv
done

#On to R

library(ggplot2)
library(dplyr)
library(tidyverse)

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/100bp")

variable_names <- c("Ara+1", "Ara+2", "Ara+3", "Ara+4", "Ara+5", "Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara-5", "Ara-6")

for (variable in variable_names) {
  
  p1 <- read_tsv(paste0(variable, "_htseq_lengths.tsv"))
  p1 <- p1 %>% unite("replseqline",c(replicate,seqtype,line),sep="_") %>%
    group_by(replseqline) %>%
    summarise(replseqline=replseqline,target_id=target_id,count=count,countsum=sum(count/length),status=status,length=length) %>%
    separate(replseqline, into = c("replicate", "seqtype","line"),sep="_")
  p1$tpm <- (((p1$count/p1$length)/(p1$countsum))*1000000)
  p1 <- p1 %>%
    group_by(target_id, status) %>%
    summarise(average_tpm = mean(tpm))
  write_csv(p1, paste0(variable, "_htseq_tpm.tsv"))
}

#Back to bash again

cat *tpm.tsv | grep -f ../50K_unambiguous_protogenes.txt | sed "1s/^/target_id,status,meantpm\n/g" > 50K_protogenes_TPMs.csv

#Back to R again

p3 <- read.csv("50K_protogenes_TPMs.csv")

p3 <- p3 %>%
  mutate(mutation_type = case_when(
    grepl("MOB", target_id, ignore.case = TRUE) ~ "IS150",
    grepl("DEL", target_id, ignore.case = TRUE) ~ "Deletion",
    grepl("SNP", target_id, ignore.case = TRUE) ~ "SNP",
    TRUE ~ NA_character_
  ))

p3 <- p3 %>%
  pivot_wider(names_from = status, values_from = meantpm, names_prefix = "tpm_", values_fill = NA)

p3 <- p3 %>%
  mutate(tpmdiff = tpm_evolved - tpm_ancestor)

p3$target_id <- sapply(str_split(p3$target_id, "_"), function(x) paste(x[1:3], collapse = "_"))

p3 %>%
  mutate(target_id = fct_reorder(target_id, tpmdiff)) %>%
  ggplot(aes(x=target_id,y=tpm_ancestor, color=mutation_type)) +
  geom_point(aes(y = tpm_ancestor),size=3,shape=16) +
  geom_point(aes(y = tpm_evolved),size=5,shape=17) +
  geom_segment(aes(x=target_id,y=tpm_ancestor,xend=target_id,yend=tpm_evolved)) +
  #  scale_color_manual(values=c("darkblue","steelblue4","lightskyblue","lightblue1","gray86")) +
  scale_color_manual(values=c("darkblue","steelblue4","blue","powderblue")) +
  scale_y_continuous(breaks = seq(0, 500, len = 11)) +
  #geom_hline(yintercept=5, col="red",linetype= "solid",size=0.5)+
  xlab("Non-genic regions with increased expression downstream of a mutation")+
  ylab("Normalized RNA counts (TPM)")+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    #    panel.grid.major = element_blank(), #remove major gridlines
    #    panel.grid.minor = element_blank(), #remove minor gridlines
    #    panel.grid.major = element_line(size = 0.5, linetype = 'dashed',colour = "gray"),
    panel.grid.minor = element_line(size = 0.5, linetype = 'dashed',colour = "gray"),
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent'),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size=16, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 16),
    legend.key.size = unit(1, 'cm'),
    legend.text = element_text(size=15),
    legend.position="top")
