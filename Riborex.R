#Getting differentially expressed translation while accounting for transcription

#Riborex

library(riborex)

getwd()

# load the kallisto results and reshape it slightly, filter tRNA, ERCC, viruses
kdf <- read_csv("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step2_gmapping/400bp/htseq_50000gen_400bp_count.csv") %>% 
  select(repl, seqtype, line, target_id, est_counts) %>% 
  unite("sample", c(seqtype, line, repl), sep = "_")

line.names <- unique(str_extract(unique(kdf$sample), "Ara[-+][1-6]"))
line.names <- line.names[!is.na(line.names)]
names(line.names) <- line.names
line.names

wide.kdf <- kdf %>% 
  pivot_wider(names_from = sample, values_from = est_counts)

head(wide.kdf)

conds <- factor(c(rep("evo", 2), rep("anc", 4)), levels = c("anc", "evo"))

rr.list <- lapply(line.names, function(x){
  ## get the RNAseq datasets
  # the rna line we're looking for
  rna.line <- paste0("rna_", x)
  
  # the RNAseq counts table, convert to matrix
  rna <- wide.kdf %>% 
    select(target_id, starts_with(rna.line), starts_with("rna_REL")) %>%   # pick the correct columns, one evo and both anc
    mutate_if(is.double, as.integer) %>%                                   # convert cols to integers
    column_to_rownames("target_id") %>%                                    # gene ids to rownames
    as.matrix()                                                            # change to a matrix
  
  ## get the riboseq datasets
  # the ribo line we're looking for
  ribo.line <- paste0("ribo_", x)
  
  # the riboseq counts table, convert to matrix
  ribo <- wide.kdf %>% 
    select(target_id, starts_with(ribo.line), starts_with("ribo_REL")) %>% #
    mutate_if(is.double, as.integer) %>% 
    column_to_rownames("target_id") %>% 
    as.matrix()
  
  ## run riborex
  rrr <- riborex(rnaCntTable = rna,
                 riboCntTable = ribo,
                 rnaCond = conds,
                 riboCond = conds)
  
  # convert to data frame with no rownames
  rdf <- as_tibble(rrr, rownames = "target_id")
  
  return(rdf)
})

rrdf <- bind_rows(rr.list, .id = "line")

#write_csv(rrdf, "/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step2_gmapping/400bp/htseq_50000gen_400bp_count_RiboRex.csv")


#MUTATIONS#

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/100bp/")

variable_names <- c("Ara+5")

for (variable in variable_names) {
  kdf2 <- read_tsv(paste0(variable, "_500_htseq.tsv")) %>% 
    select(replicate, seqtype, line, target_id, est_counts) %>% 
    unite("sample", c(seqtype, line, replicate), sep = "_")
  
  line.names <- unique(str_extract(unique(kdf2$sample), "Ara[-+][1-6]"))
  line.names <- line.names[!is.na(line.names)]
  names(line.names) <- line.names
  
  wide.kdf2 <- kdf2 %>% 
    pivot_wider(names_from = sample, values_from = est_counts)
  
  conds <- factor(c(rep("evo", 2), rep("anc", 4)), levels = c("anc", "evo"))
  
  rr.list <- lapply(line.names, function(x){
    ## get the RNAseq datasets
    # the rna line we're looking for
    rna.line <- paste0("rna_", x)
    
    # the RNAseq counts table, convert to matrix
    rna <- wide.kdf2 %>% 
      select(target_id, starts_with(rna.line), starts_with("rna_REL")) %>%   # pick the correct columns, one evo and both anc
      mutate_if(is.double, as.integer) %>%                                   # convert cols to integers
      column_to_rownames("target_id") %>%                                    # gene ids to rownames
      as.matrix()                                                            # change to a matrix
    
    ## get the riboseq datasets
    # the ribo line we're looking for
    ribo.line <- paste0("ribo_", x)
    
    # the riboseq counts table, convert to matrix
    ribo <- wide.kdf2 %>% 
      select(target_id, starts_with(ribo.line), starts_with("ribo_REL")) %>% #
      mutate_if(is.double, as.integer) %>% 
      column_to_rownames("target_id") %>% 
      as.matrix()
    
    ## run riborex
    rrr <- riborex(rnaCntTable = rna,
                   riboCntTable = ribo,
                   rnaCond = conds,
                   riboCond = conds)
    
    # convert to data frame with no rownames
    rdf <- as_tibble(rrr, rownames = "target_id")
    
    return(rdf)
    
  })
  
  rrdf <- bind_rows(rr.list, .id = "line")
  
  write_csv(rrdf, paste0(variable, "_riborex_results.csv"))
}


setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/1010_coverage_plots/del_25_AraP1/")

p1 <- read.csv("AraP1_bothrep_genomecov_relevant_tidy_mod.csv")

get_average_high_low <- function(data) {
  result <- data %>%
    group_by(position) %>%
    summarize(
      average = mean(count),
      high_value = max(count),
      low_value = min(count),
      sd = sd(count),
      n = n(),
      deletion_status = deletion_status,
    )
  return(result)
}

p1_evolved <- p1 %>% filter(evo_status == "evolved")

ggplot(get_average_high_low(p1_evolved), aes(x = position, y = average)) +
  geom_point(aes(y = high_value, color = deletion_status), size = 2, position = position_dodge(width = 0.2)) +
  geom_point(aes(y = low_value, color = deletion_status), size = 2, position = position_dodge(width = 0.2)) +
  geom_segment(aes(x = position, xend = position, y = high_value, yend = low_value, color = deletion_status), 
               size = 0.3, linetype = "dashed", position = position_dodge(width = 0.2)) +
  geom_line(aes(color = deletion_status), size = 1.5, position = position_dodge(width = 0.2)) +
  labs(x = "Position", y = "RNAseq coverage (read count)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  ylim(0,50)+
  xlim(375,1000) +
  scale_color_manual(values = c("pre-deletion" = "#C8A2C8", "post-deletion" = "#FFDAB9"), 
                     guide = guide_legend(title = "Deletion Status")) +
  scale_color_manual(values = c("pre-deletion" = "darkviolet", "post-deletion" = "darkorange1"), 
                     aesthetics = c("geom_line"), 
                     guide = FALSE)

ggsave("AraP1_evolved.pdf", width = 15, height = 8, units = "in", dpi = 300)

p1_ancestor <- p1 %>% filter(evo_status == "ancestor")

ggplot(get_average_high_low(p1_ancestor), aes(x = position, y = average)) +
  geom_point(aes(y = high_value, color = deletion_status), size = 2, position = position_dodge(width = 0.2)) +
  geom_point(aes(y = low_value, color = deletion_status), size = 2, position = position_dodge(width = 0.2)) +
  geom_segment(aes(x = position, xend = position, y = high_value, yend = low_value, color = deletion_status), 
               size = 0.3, linetype = "dashed", position = position_dodge(width = 0.2)) +
  geom_line(aes(color = deletion_status), size = 1.5, position = position_dodge(width = 0.2)) +
  labs(x = "Position", y = "RNAseq coverage (read count)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  ylim(0,50)+
  xlim(375,1024) +
  scale_color_manual(values = c("pre-deletion" = "#C8A2C8", "post-deletion" = "#FFDAB9"), 
                     guide = guide_legend(title = "Deletion Status")) +
  scale_color_manual(values = c("pre-deletion" = "darkviolet", "post-deletion" = "darkorange1"), 
                     aesthetics = c("geom_line"), 
                     guide = FALSE)

ggsave("AraP1_ancestor.pdf", width = 15, height = 8, units = "in", dpi = 300)
