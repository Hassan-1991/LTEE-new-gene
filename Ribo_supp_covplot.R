#####Ribocoverage_minus#####

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(grid)

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/")

values <- c(
  "Ara+3_1110_SNP_minus_ribo",
  "Ara+3_1238_SNP_minus_ribo",
  "Ara+3_967_SNP_minus_ribo",
  "Ara+3_988_INS_minus_ribo",
  "Ara+4_32_SNP_minus_ribo",
  "Ara+5_76_MOB_minus_ribo",
  "Ara-2_659_SNP_minus_ribo",
  "Ara-2_777_SNP_minus_ribo",
  "Ara-3_454_DEL_minus_ribo",
  "Ara-3_468_SNP_minus_ribo",
  "Ara-3_725_MOB_minus_ribo",
  "Ara-3_776_MOB_minus_ribo",
  "Ara-4_772_INS_minus_ribo",
  "Ara-4_871_SNP_minus_ribo",
  "Ara-5_54_DEL_minus_ribo"
)  

get_average_high_low <- function(data) {
  result <- data %>%
    filter(seqtype == "ribo") %>%
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
generate_plots <- function(X, ylim_value) {
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
    labs(x = "Position", y = "Coverage (read count)") +
    scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    scale_x_reverse() +
    ylim(0, ylim_value) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  collage <- grid.arrange(plot1, plot3, ncol = 2, top=textGrob(X, gp=gpar(fontsize=25)))
  return(collage)
}

# Define the ylim values
ylim_values <- c(350, 20, 500, 35, 150, 400, 15, 400, 40, 45, 30, 30, 50, 600, 25)

# Create plots with different ylim values
plot_list <- mapply(generate_plots, values, ylim_values, SIMPLIFY = FALSE)

# Arrange and print the final plot
final_plot <- do.call(grid.arrange, c(plot_list, ncol = 3))
print(final_plot)

ggsave("Ribocoverage_minus1.pdf", final_plot, width = 80, height = 80, units = "cm",limitsize=FALSE)

#####RiboCoverage_plus#####

values <- c(
  "Ara+1_11_SNP_plus_ribo",
  "Ara+2_4_SNP_plus_ribo",
  "Ara+3_206_SNP_plus_ribo",
  "Ara+3_664_SNP_plus_ribo",
  "Ara-1_644_SNP_plus_ribo",
  "Ara-2_202_SNP_plus_ribo",
  "Ara-2_402_INS_plus_ribo",
  "Ara-3_462_SNP_plus_ribo",
  "Ara-3_71_SNP_plus_ribo",
  "Ara-4_151_SNP_plus_ribo",
  "Ara-4_255_SNP_plus_ribo",
  "Ara-6_13_DEL_plus_ribo",
  "Ara-6_59_DEL_plus_ribo",
  "Ara-6_84_SNP_plus_ribo")  

get_average_high_low <- function(data) {
  result <- data %>%
    #filter(seqtype == "rna") %>%
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
generate_plots <- function(X, ylim_value) {
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
    labs(x = "Position", y = "Coverage (read count)") +
    #scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    #scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    #scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    #scale_x_reverse() +
    ylim(0, ylim_value) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  collage <- grid.arrange(plot2, plot4, ncol = 2, top=textGrob(X, gp=gpar(fontsize=25)))
  return(collage)
}

# Define the ylim values
ylim_values <- c(45, 45, 80, 60, 120, 800, 150, 450, 150, 55, 800,3000,250,300)

# Create plots with different ylim values
plot_list <- mapply(generate_plots, values, ylim_values, SIMPLIFY = FALSE)

# Arrange and print the final plot
final_plot <- do.call(grid.arrange, c(plot_list, ncol = 3))
print(final_plot)

ggsave("Ribocoverage_plus1.pdf", final_plot, width = 80, height = 80, units = "cm",limitsize=FALSE)


#####RNAseq_plus#####

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(grid)

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/all_coverages")

values <- c(
  "Ara+1_31_MOB",
  "Ara+1_35_MOB",
  "Ara+1_102_MOB",
  "Ara+5_30_MOB",
  "Ara-5_84_MOB",
  "Ara-6_24_MOB",
  "Ara+1_11_SNP",
  "Ara+3_206_SNP",
  "Ara+3_253_SNP",
  "Ara+4_52_DEL",
  "Ara-1_704_SNP",
  "Ara-1_85_SNP",
  "Ara-2_749_SNP",
  "Ara-3_71_SNP",
  "Ara-5_67_SNP",
  "Ara-6_25_SNP",
  "Ara+4_47_SNP",
  "Ara-1_644_SNP",
  "Ara-2_920_SNP",
  "Ara+3_664_SNP",
  "Ara+5_54_SNP",
  "Ara-2_1037_DEL",
  "Ara-2_489_INS",
  "Ara-2_118_SNP",
  "Ara-2_402_INS",
  "Ara-6_59_DEL",
  "Ara-2_165_INS",
  "Ara-2_437_INS",
  "Ara-3_462_SNP",
  "Ara-4_288_DEL",
  "Ara-6_13_DEL"
)  

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
generate_plots <- function(X, ylim_value) {
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
    labs(x = "Position", y = "Coverage (read count)") +
    #scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    #scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    #scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    #scale_x_reverse() +
    ylim(0, ylim_value) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  collage <- grid.arrange(plot2, plot4, ncol = 2, top=textGrob(X, gp=gpar(fontsize=25)))
  return(collage)
}

# Define the ylim values
ylim_values <- c(30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 60,
                 60,
                 60,
                 60,
                 60,
                 60,
                 60,
                 250,
                 250,
                 250,
                 250,
                 250,
                 600,
                 600,
                 2500)

# Create plots with different ylim values
plot_list <- mapply(generate_plots, values, ylim_values, SIMPLIFY = FALSE)

# Arrange and print the final plot
final_plot <- do.call(grid.arrange, c(plot_list, ncol = 3))
print(final_plot)

ggsave("RNAcoverage_plus1.pdf", final_plot, width = 80, height = 120, units = "cm",limitsize=FALSE)

#####RNAseq_minus#####

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(grid)

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/all_coverages")

values <- c(
  "Ara+1_99_MOB",
  "Ara+1_120_MOB",
  "Ara+5_65_MOB",
  "Ara-3_547_MOB",
  "Ara-3_776_MOB",
  "Ara-6_33_MOB",
  "Ara+1_75_DEL",
  "Ara+1_76_DEL",
  "Ara+3_988_INS",
  "Ara-2_440_SNP",
  "Ara-2_497_SNP",
  "Ara-2_610_SNP",
  "Ara-2_991_SNP",
  "Ara-3_454_DEL",
  "Ara-4_308_INS",
  "Ara-4_772_INS",
  "Ara+5_76_MOB",
  "Ara+1_94_DEL",
  "Ara+1_47_SNP",
  "Ara+3_1403_SNP",
  "Ara-2_409_DEL",
  "Ara-2_716_SNP",
  "Ara-2_842_DEL",
  "Ara-3_136_SNP",
  "Ara-3_387_DEL",
  "Ara-4_1061_DEL",
  "Ara-4_48_SNP",
  "Ara-2_583_INS",
  "Ara-2_731_SNP",
  "Ara+2_35_SNP",
  "Ara+4_49_SNP",
  "Ara-1_37_SNP"
)  

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
generate_plots <- function(X, ylim_value) {
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
    labs(x = "Position", y = "Coverage (read count)") +
    scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    scale_x_reverse() +
    ylim(0, ylim_value) +
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
    labs(x = "Position", y = "Coverage (read count)") +
    scale_x_reverse() +
    ylim(0, ylim_value) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  collage <- grid.arrange(plot1, plot3, ncol = 2, top=textGrob(X, gp=gpar(fontsize=25)))
  return(collage)
}

# Define the ylim values
ylim_values <- c(30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 30,
                 60,
                 60,
                 60,
                 60,
                 60,
                 60,
                 60,
                 60,
                 60,
                 60,
                 60,
                 200,
                 200,
                 200,
                 200,
                 200)

# Create plots with different ylim values
plot_list <- mapply(generate_plots, values, ylim_values, SIMPLIFY = FALSE)

# Arrange and print the final plot
final_plot <- do.call(grid.arrange, c(plot_list, ncol = 3))
print(final_plot)

ggsave("RNAcoverage_minus1.pdf", final_plot, width = 80, height = 120, units = "cm",limitsize=FALSE)

#####MURI_plus#####

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)

setwd("/stor/work/Ochman/hassan/Fall_2022/RNAseq/")

# List of X values
X_values <- c(
  "Ara+1_31_MOB",
  "Ara-5_84_MOB",
  "Ara+1_35_MOB",
  "Ara-6_59_DEL",
  "Ara-1_704_SNP",
  "Ara-2_118_SNP"
)

ylim_values <- c(15, 15, 15, 150, 15, 100)

# Function to generate plots
generate_plot <- function(X,ylim) {
  p1 <- read.csv(paste0(X, "_right"), sep = '\t') %>% filter(!startsWith(line, "Glucose"))
  
  # Filter out lines where all points have a value of 0
  p1_filtered <- p1 %>%
    group_by(line) %>%
    filter(any(count != 0))
  
  plot1 <- ggplot(p1_filtered, aes(x = position, y = count, color = line)) +
    geom_line() +
    labs(x = "Position", y = "Coverage (read count)", title = X) +
    ylim(0, ylim) +
    theme_minimal() +
    theme(
      legend.position = "right",  # Show legend
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  
  return(plot1)
}

# Create a list of plots
plot_list <- mapply(generate_plot, X_values, ylim_values, SIMPLIFY = FALSE)

# Arrange the plots in four columns
collage <- grid.arrange(grobs = plot_list, ncol = 4)

# Print the collage
print(collage)


########MURI_minus#####

# List of X values
X_values <- c(
  "Ara-6_33_MOB",
  "Ara+5_76_MOB",
  "Ara+1_76_DEL",
  "Ara-2_440_SNP",
  "Ara-2_731_SNP",
  "Ara-2_991_SNP"
)

ylim_values <- c(15, 15, 15, 25, 125, 25)

# Function to generate plots
generate_plot <- function(X, ylim_value) {
  p1 <- read.csv(paste0(X, "_left"), sep = '\t') %>% filter(!startsWith(line, "Glucose"))
  
  # Filter out lines where all points have a value of 0
  p1_filtered <- p1 %>%
    group_by(line) %>%
    filter(any(count != 0))
  
  plot1 <- ggplot(p1_filtered, aes(x = position, y = count, color = line)) +
    geom_line() +
    labs(x = "Position", y = "Coverage (read count)", title = X) +
    ylim(0, ylim_value) +
    scale_x_reverse() +
    theme_minimal() +
    theme(
      legend.position = "right",  # Show legend
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16)
    )
  
  return(plot1)
}

# Create a list of plots
plot_list <- mapply(generate_plot, X_values, ylim_values, SIMPLIFY = FALSE)

# Arrange the plots in four columns
collage <- grid.arrange(grobs = plot_list, ncol = 4)

# Print the collage
print(collage)

ggsave("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/collage.pdf", collage, width = 80, height = 80, units = "cm", limitsize=FALSE)






############TIMESERIES###############

#This is a lot of manual work

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/timeseries")

#X <- "MOB_1776434_IS150"
#X <- "MOB_4110237_IS150"
#X <- "MOB_4415710_IS150"
#X <- "SNP_2446984_C"
#X <- "SNP_3475057_T"
#X <- "MOB_3015771_IS150"
#X <- "MOB_1101970_IS150"
#X <- "MOB_1462266_IS150"
X <- "INS_2103915_CAGCCAGCCAGCCAGCCAGCCAGC"

Y <- c("REL606",
       "Ara-3_5000gen_ZDB409",
       "Ara-3_10000gen_ZDB425",
       "Ara-3_15000gen_ZDB445",
       "Ara-3_20000gen_ZDB467",
       #"Ara-3_25000gen_ZDB483",
       #"Ara-3_25000gen_ZDB486",
       "Ara-3_25000gen_ZDB488",
       #"Ara-3_27000gen_ZDB309",
       #"Ara-3_30000gen_ZDB17",
       "Ara-3_31500gen_ZDB199",
       "Ara-3_31500gen_ZDB200")#,
       #"Ara-3_31500gen_ZDB25",
       #"Ara-3_31500gen_ZDB564")#,
       #"Ara-3_33000gen_CZB154"#,
       #"Ara-3_50000gen_11364")

for (value in Y) {
  p1 <- read.csv(paste0(X, "_down_coverage_", value), sep = '\t')
  
  create_plot <- function(data) {
    result <- data %>%
      group_by(position) %>%
      summarize(
        average = mean(count),
        high_value = max(count),
        low_value = min(count),
        sd = sd(count),
        n = n()
      )
    
    ggplot(result, aes(x = position, y = average)) +
      geom_point(aes(y = high_value), color = "#FFDAB9", size = 2, position = position_dodge(width = 0.2)) +
      geom_point(aes(y = low_value), color = "#FFDAB9", size = 2, position = position_dodge(width = 0.2)) +
      geom_segment(aes(x = position, xend = position, y = high_value, yend = low_value), 
                   color = "#FFDAB9", size = 0.3, linetype = "dashed", position = position_dodge(width = 0.2)) +
      geom_line(color = "darkorange1", size = 1.5, position = position_dodge(width = 0.2)) +
      labs(x = "Position", y = "Coverage (read count)") +
      #ylim(0,400) + #MOB_1776434_IS150.pdf
      #ylim(0,30) + #MOB_4110237_IS150.pdf
      #ylim(0,60) + #MOB_4415710_IS150.pdf
      #ylim(0,30) + #SNP_2446984_C.pdf
      #ylim(0,30) + #SNP_3475057_T.pdf
      #ylim(0,60) + #MOB_3015771_IS150.pdf
      #ylim(0,30) + #MOB_1101970_IS150.pdf
      #ylim(0,60) + #MOB_1462266_IS150.pdf
      ylim(0,30) + #INS_2103915_CAGCCAGCCAGCCAGCCAGCCAGC.pdf
      theme_minimal() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))
  }
  
  plot <- create_plot(p1)
  plots_list[[value]] <- plot
}

#spacing <- unit(5, "lines")  # Adjust this value as needed

# Create the final collage with adjusted spacing
#collage <- grid.arrange(grobs = plots_list, ncol = 6, top = X, heights = rep(unit(2, "null"), 3), widths = rep(unit(2, "null"), 6), padding = spacing)

collage <- grid.arrange(grobs = plots_list[Y], ncol = 5)

# Display the collage
print(collage)

#ggsave("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/newTS_plots/MOB_1776434_IS150.pdf", collage, width = 60, height = 40, units = "cm", limitsize=FALSE)
#ggsave("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/newTS_plots/MOB_4110237_IS150.pdf", collage, width = 60, height = 40, units = "cm", limitsize=FALSE)
#ggsave("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/newTS_plots/MOB_4415710_IS150.pdf", collage, width = 60, height = 40, units = "cm", limitsize=FALSE)
#ggsave("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/newTS_plots/SNP_2446984_C.pdf", collage, width = 60, height = 40, units = "cm", limitsize=FALSE)
#ggsave("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/newTS_plots/SNP_3475057_T.pdf", collage, width = 60, height = 40, units = "cm", limitsize=FALSE)
#ggsave("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/newTS_plots/MOB_3015771_IS150.pdf", collage, width = 60, height = 40, units = "cm", limitsize=FALSE)
#ggsave("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/newTS_plots/MOB_1101970_IS150.pdf", collage, width = 60, height = 40, units = "cm", limitsize=FALSE)
#ggsave("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/newTS_plots/MOB_1462266_IS150.pdf", collage, width = 60, height = 40, units = "cm", limitsize=FALSE)
#ggsave("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/newTS_plots/INS_2103915_CAGCCAGCCAGCCAGCCAGCCAGC.pdf", collage, width = 60, height = 40, units = "cm", limitsize=FALSE)
