###FIGURE 6A###

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/")

p1 <- read.csv("Ara-2_731_SNP_minus_coverage_casestudy", sep = '\t')

p1_filtered <- p1 %>% filter(status=="ancestor")

plot1 <- ggplot(get_average_high_low(p1_filtered), aes(x = position, y = average)) +
  geom_point(aes(y = high_value), color = "#C8A2C8", size = 2, position = position_dodge(width = 0.2)) +
  geom_point(aes(y = low_value), color = "#C8A2C8", size = 2, position = position_dodge(width = 0.2)) +
  geom_segment(aes(x = position, xend = position, y = high_value, yend = low_value), 
               color = "#C8A2C8", size = 0.8, linetype = "dashed", position = position_dodge(width = 0.2)) +
  geom_line(color = "darkviolet", size = 1.5, position = position_dodge(width = 0.2)) +
  labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, upstream") +
  theme_minimal() +
  ylim(0,120) +
  scale_x_reverse() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  )

plot1

p2_filtered <- p1 %>% filter(status=="evolved")

plot2 <- ggplot(get_average_high_low(p2_filtered), aes(x = position, y = average)) +
  geom_point(aes(y = high_value), color = "#FFDAB9", size = 2, position = position_dodge(width = 0.2)) +
  geom_point(aes(y = low_value), color = "#FFDAB9", size = 2, position = position_dodge(width = 0.2)) +
  geom_segment(aes(x = position, xend = position, y = high_value, yend = low_value), 
               color = "#FFDAB9", size = 0.8, linetype = "dashed", position = position_dodge(width = 0.2)) +
  geom_line(color = "darkorange1", size = 1.5, position = position_dodge(width = 0.2)) +
  labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, upstream") +
  theme_minimal() +
  ylim(0,120) +
  scale_x_reverse() +
theme(
  legend.position = "none",
  plot.title = element_text(hjust = 0.5, size = 20),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text.x = element_text(size = 16),
  axis.text.y = element_text(size = 16)
)

plot2

ggsave("021324_6A_anc.pdf", plot = plot1, device = "pdf", bg = "transparent", width=15, height=10)
ggsave("021324_6A_evo.pdf", plot = plot2, device = "pdf", bg = "transparent", width=15, height=10)


###FIGURE 6B###

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/1010_coverage_plots/del_25_AraP1/")

#How was this file generated?

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
