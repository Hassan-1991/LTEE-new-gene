library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(grid)

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/all_coverages/")

values <- c(
  "Ara+1_99_MOB",
  "Ara+5_30_MOB",
  "Ara+5_65_MOB",
  "Ara-3_547_MOB",
  "Ara-6_24_MOB")  

get_average_high_low <- function(data) {
  result <- data %>%
    filter(seqtype == "rna") %>%
    group_by(position) %>%
    summarize(
      average = mean(count),
      high_value = max(count),
      low_value = min(count),
      sd = sd(count),
      n = n()
    )
  return(result)
}

# Create a function to generate the plots
generate_plots <- function(X) {
  p1 <- read.csv(paste0(X, "_coverage_ancestor_left"), sep = '\t')
  p2 <- read.csv(paste0(X, "_coverage_ancestor_right"), sep = '\t')
  p3 <- read.csv(paste0(X, "_coverage_left"), sep = '\t')
  p4 <- read.csv(paste0(X, "_coverage_right"), sep = '\t')
  
  plot1 <- ggplot(get_average_high_low(p1), aes(x = position, y = average)) +
    geom_point(aes(y = high_value), color = "#C8A2C8", size = 2, position = position_dodge(width = 0.2)) +
    geom_point(aes(y = low_value), color = "#C8A2C8", size = 2, position = position_dodge(width = 0.2)) +
    geom_segment(aes(x = position, xend = position, y = high_value, yend = low_value), 
                 color = "#C8A2C8", size = 0.3, linetype = "dashed", position = position_dodge(width = 0.2)) +
    geom_line(color = "darkviolet", size = 1.5, position = position_dodge(width = 0.2)) +
    labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, upstream") +
    ylim(0,20) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  
  plot2 <- ggplot(get_average_high_low(p2), aes(x = position, y = average)) +
    geom_point(aes(y = high_value), color = "#C8A2C8", size = 2, position = position_dodge(width = 0.2)) +
    geom_point(aes(y = low_value), color = "#C8A2C8", size = 2, position = position_dodge(width = 0.2)) +
    geom_segment(aes(x = position, xend = position, y = high_value, yend = low_value), 
                 color = "#C8A2C8", size = 0.3, linetype = "dashed", position = position_dodge(width = 0.2)) +
    geom_line(color = "darkviolet", size = 1.5, position = position_dodge(width = 0.2)) +
    labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, downstream") +
    ylim(0,20) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  
  plot3 <- ggplot(get_average_high_low(p3), aes(x = position, y = average)) +
    geom_point(aes(y = high_value), color = "#FFDAB9", size = 2, position = position_dodge(width = 0.2)) +
    geom_point(aes(y = low_value), color = "#FFDAB9", size = 2, position = position_dodge(width = 0.2)) +
    geom_segment(aes(x = position, xend = position, y = high_value, yend = low_value), 
                 color = "#FFDAB9", size = 0.3, linetype = "dashed", position = position_dodge(width = 0.2)) +
    geom_line(color = "darkorange1", size = 1.5, position = position_dodge(width = 0.2)) +
    labs(x = "Position", y = "Coverage (read count)", title = "Evolved, upstream") +
    ylim(0,20) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  
  plot4 <- ggplot(get_average_high_low(p4), aes(x = position, y = average)) +
    geom_point(aes(y = high_value), color = "#FFDAB9", size = 2, position = position_dodge(width = 0.2)) +
    geom_point(aes(y = low_value), color = "#FFDAB9", size = 2, position = position_dodge(width = 0.2)) +
    geom_segment(aes(x = position, xend = position, y = high_value, yend = low_value), 
                 color = "#FFDAB9", size = 0.3, linetype = "dashed", position = position_dodge(width = 0.2)) +
    geom_line(color = "darkorange1", size = 1.5, position = position_dodge(width = 0.2)) +
    labs(x = "Position", y = "Coverage (read count)", title = "Evolved, downstream") +
    ylim(0,20) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  collage <- grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, top=textGrob(X, gp=gpar(fontsize=25)))
  return(collage)
}

plot_list <- lapply(values, generate_plots)
final_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))
print(final_plot)

ggsave("MOB_2.pdf", final_plot, width = 80, height = 90, units = "cm")
