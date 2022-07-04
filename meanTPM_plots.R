#50000gen

test1 <- read.csv("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_50000gen_400bp_count_mark1.csv")
test1$countsum <- test1$est_counts
test1 <- test1 %>% unite("replseqline",c(replicate,seqtype,line),sep="_") %>%
  group_by(replseqline) %>%
  summarise(replseqline=replseqline,target_id=target_id,est_counts=est_counts,countsum=sum(countsum/400),type=type,rt_status=rt_status) %>%
  separate(replseqline, into = c("replicate", "seqtype","line"),sep="_")
test1$tpm <- (((test1$est_counts/400)/(test1$countsum))*1000000)
test2 <- test1 %>%
  unite("seqline",c(seqtype,line),sep="_") %>%
  unite("sample",c(seqline,target_id),sep="%") %>%
  select(sample,replicate,tpm,type,rt_status) %>%
  group_by(sample) %>%
  summarise(sample=sample,tpm=mean(tpm),type=type,rt_status=rt_status) %>%
  unique() %>%
  separate(sample,into=c("seqline","target_id"),sep="%") %>%
  separate(seqline,into=c("seqtype","line"),sep="_")

write.csv(test2,"/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_50000gen_400bp_count_mark1_tpm.csv")

test2$above5PM <- "NO"
test2$above5PM[test2$tpm>5] <- "YES"
test2$above1PM <- "NO"
test2$above1PM[test2$tpm>1] <- "YES"

test2_nons <- test2 %>%
  filter(seqtype=="ribo") %>%
  #filter(seqtype=="rna") %>%
  filter(above5PM=="YES")
test2_nons$line <- factor(test2_nons$line,levels = c("REL606","REL607","Ara+1","Ara+2","Ara+3","Ara+4","Ara+5","Ara-1","Ara-2","Ara-3","Ara-4","Ara-5","Ara-6"))
test2_nons$seqtype <- factor(test2_nons$seqtype,levels = c("rna","ribo"))
test2_nons$rt_status <- factor(test2_nons$rt_status,levels = c("readthrough","noncoding"))

test2_nons %>%
  filter(rt_status!="annotated") %>%
  #filter(rt_status!="readthrough") %>%
  ggplot(aes(x=line,fill=rt_status)) +
  geom_bar(position="stack", color = "azure4",width=0.75) +
  #scale_fill_manual(values=c("grey67","cornflowerblue")) +
  scale_fill_manual(values=c("lightcyan1","cornflowerblue")) +
  ylim(0,300) +
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
  xlab("Lines")+
  ylab("Non-coding windows with normalized read counts (TPM) > 5")
  
#timeseries

test1 <- read.csv("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_time_series_400bp_namesorted_count_mark2.csv")
test1$countsum <- test1$est_counts
test1 <- test1 %>% unite("genstrainrepl",c(generation,strain,replicate),sep="%") %>% 
  group_by(genstrainrepl) %>% 
  summarise(genstrainrepl=genstrainrepl,target_id=target_id,est_counts=est_counts,countsum=sum(countsum/400),type=type,rt_status=rt_status) %>%
  separate(genstrainrepl, into = c("generation", "strain","replicate"),sep="%")
test1$tpm <- (((test1$est_counts/400)/(test1$countsum))*1000000)
test2 <- test1 %>% 
  unite("genstrain",c(generation,strain),sep="@") %>% 
  unite("sample",c(genstrain,target_id),sep="%") %>% 
  select(sample,replicate,tpm,type,rt_status) %>% 
  group_by(sample) %>% 
  summarise(sample=sample,tpm=mean(tpm),type=type,rt_status=rt_status) %>% 
  unique() %>% 
  separate(sample,into=c("genstrain","target_id"),sep="%") %>% 
  separate(genstrain,into=c("generation","strain"),sep="@")

write.csv(test2,"/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_time_series_400bp_namesorted_count_mark2_tpm.csv")

test2$above5PM <- "NO"
test2$above5PM[test2$tpm>5] <- "YES"
test2$above1PM <- "NO"
test2$above1PM[test2$tpm>1] <- "YES"

test2_nons <- test2 %>%
  #filter(seqtype=="ribo") %>%
  #filter(seqtype=="rna") %>%
  filter(above5PM=="YES")
test2_nons$strain <- factor(test2_nons$strain,levels = c("Ancestor_1","Ancestor_2","5000","10000","15000","20000","25000_1","25000_2","25000_3","27000","30000","31500_1","31500_2","31500_3","31500_4","33000"))
#test2_nons$seqtype <- factor(test2_nons$seqtype,levels = c("rna","ribo"))
test2_nons$rt_status <- factor(test2_nons$rt_status,levels = c("readthrough","noncoding"))

test2_nons %>%
  filter(rt_status!="annotated") %>%
  #filter(rt_status!="readthrough") %>%
  ggplot(aes(x=strain,fill=rt_status)) +
  geom_bar(position="stack", color = "azure4",width=0.75) +
  #scale_fill_manual(values=c("grey67","cornflowerblue")) +
  scale_fill_manual(values=c("lightcyan1","cornflowerblue")) +
  #ylim(0,300) +
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
  xlab("Lines")+
  ylab("Non-coding windows with normalized read counts (TPM) > 5")
