#Coverage plots with riboseq data
#Lots of manual work, especially to figure out and assign ylim values

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(grid)

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/")

#First ten for minus:

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
  "Ara-3_468_SNP_minus_ribo"
)  

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
    labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, upstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, downstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Evolved, upstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Evolved, downstream") +
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
ylim_values <- c(350, 20, 500, 35, 150, 400, 15, 12000, 40, 45)

# Create plots with different ylim values
plot_list <- mapply(generate_plots, values, ylim_values, SIMPLIFY = FALSE)

# Arrange and print the final plot
final_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))
print(final_plot)

ggsave("ribo_supp_minus1.pdf", final_plot, width = 80, height = 80, units = "cm",limitsize=FALSE)

#Next 5 for minus:

values <- c(
  "Ara-3_725_MOB_minus_ribo",
  "Ara-3_776_MOB_minus_ribo",
  "Ara-4_772_INS_minus_ribo",
  "Ara-4_871_SNP_minus_ribo",
  "Ara-5_54_DEL_minus_ribo"
)  

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
    labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, upstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, downstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Evolved, upstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Evolved, downstream") +
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
ylim_values <- c(10, 10, 50, 800, 30)

# Create plots with different ylim values
plot_list <- mapply(generate_plots, values, ylim_values, SIMPLIFY = FALSE)

# Arrange and print the final plot
final_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))
print(final_plot)

ggsave("ribo_supp_minus2.pdf", final_plot, width = 80, height = 48, units = "cm",limitsize=FALSE)

##First ten plus:

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
  "Ara-4_151_SNP_plus_ribo")  

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
    labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, upstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, downstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Evolved, upstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Evolved, downstream") +
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
ylim_values <- c(45, 45, 80, 60, 120, 800, 150, 350, 150, 55)

# Create plots with different ylim values
plot_list <- mapply(generate_plots, values, ylim_values, SIMPLIFY = FALSE)

# Arrange and print the final plot
final_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))
print(final_plot)

#Last four plus:

values <- c(
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
    labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, upstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Ancestor, downstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Evolved, upstream") +
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
    labs(x = "Position", y = "Coverage (read count)", title = "Evolved, downstream") +
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
ylim_values <- c(800,3000,250,300)

# Create plots with different ylim values
plot_list <- mapply(generate_plots, values, ylim_values, SIMPLIFY = FALSE)

# Arrange and print the final plot
final_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))
print(final_plot)
