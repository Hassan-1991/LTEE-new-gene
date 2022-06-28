#Time_series_dataset

library(tidyverse)
library(DESeq2)
library(matrixStats)
library(apeglm)
kdf3 <- read_csv("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_time_series_400bp_count.csv") %>% 
  select(generation, strain, replicate, target_id, est_counts) %>% 
  unite("sample", c(generation, strain, replicate), sep = "_")
line.names2 <- unique(str_extract(unique(kdf3$sample),"gen[0-9]{2,5}_[A-Z]{3}[0-9]{2,4}"))
line.names2 <- line.names2[!is.na(line.names2)]
names(line.names2) <- line.names2
wide.kdf3 <- kdf3 %>% 
  pivot_wider(names_from = sample, values_from = est_counts)
wide.kdf4 <- wide.kdf3 %>%
  rowwise %>%
  mutate(row.sum = sum(c_across(starts_with("gen")))) %>%
  filter(row.sum>5) %>%
  select(1:49)
df.list2 <- lapply(line.names2, function(x){
  anc.to.pick <- paste("gen0", sep = "_")
  wide.kdf4 %>%
    select(target_id, starts_with(x), starts_with(anc.to.pick)) %>% 
    mutate_if(is.double, as.integer) %>%                            
    column_to_rownames("target_id") %>%                             
    as.matrix()                                                     
})
(conds.df2 <- data.frame(condition = factor(c(rep("evo", 3), rep("anc", 6)), levels = c("anc", "evo"))))
deseq.list <- lapply(df.list2, function(x){
  d1 <- DESeqDataSetFromMatrix(countData = x,
                               colData = conds.df2,
                               design = ~condition)
  d2 <- DESeq(d1)
  d3 <- lfcShrink(d2,
                  coef = 2,
                  type = "apeglm")
  d4 <- as_tibble(d3, rownames = "target_id")
  return(d4)
})
deseq.df <- bind_rows(deseq.list, .id ="sample") %>% 
  separate(sample, into = c("generation", "strain"), sep = "_")
write_csv(deseq.df, "/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_time_series_400bp_foldchanges.csv")

#50000gen_dataset

kdf2 <- read_csv("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_50000gen_400bp_count.csv") %>% 
  select(replicate, seqtype, line, target_id, est_counts) %>% 
  unite("sample", c(seqtype, line, replicate), sep = "_")
line.names <- unique(str_extract(unique(kdf2$sample),"[a-z]{3,4}_Ara[MP][1-6]"))
line.names <- line.names[!is.na(line.names)]
names(line.names) <- line.names
wide.kdf2 <- kdf2 %>% 
  pivot_wider(names_from = sample, values_from = est_counts)
wide.kdf3 <- wide.kdf2 %>%
  rowwise %>%
  mutate(row.sum = sum(c_across(starts_with("r")))) %>%
  filter(row.sum>10) %>%
  select(1:53)
df.list <- lapply(line.names, function(x){
  seqtype <- unlist(str_split(x, "_"))[1]
  anc.to.pick <- paste(seqtype, "REL", sep = "_")
  wide.kdf3 %>%
    select(target_id, starts_with(x), starts_with(anc.to.pick)) %>% 
    mutate_if(is.double, as.integer) %>%                            
    column_to_rownames("target_id") %>%                             
    as.matrix()                                                     
})
(conds.df <- data.frame(condition = factor(c(rep("evo", 2), rep("anc", 4)), levels = c("anc", "evo"))))
deseq.list2 <- lapply(df.list, function(x){
  d1 <- DESeqDataSetFromMatrix(countData = x,
                               colData = conds.df,
                               design = ~condition)  
  d2 <- DESeq(d1)
  d3 <- lfcShrink(d2,
                  coef = 2,
                  type = "apeglm")
  d4 <- as_tibble(d3, rownames = "target_id")
  return(d4)
})
deseq.df2 <- bind_rows(deseq.list2, .id ="sample") %>% 
  separate(sample, into = c("seqtype", "line"), sep = "_")
write_csv(deseq.df2, "/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_50000gen_400bp_foldchanges.csv")
