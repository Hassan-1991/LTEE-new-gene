#TO BE FIXED

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(grid)

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/timeseries/")

X <- "DEL_3015771_19352"  # Replace with the desired value of X
#X <- "MOB_4415710_IS150"  # Replace with the desired value of X
#X <- "MOB_1776434_IS150"  # Replace with the desired value of X
Y <- c("REL606",
       #"Ara-3_5000gen_ZDB409",
       "Ara-3_10000gen_ZDB425",
       #"Ara-3_15000gen_ZDB445",
       "Ara-3_20000gen_ZDB467",
       #"Ara-3_25000gen_ZDB483",
       #"Ara-3_25000gen_ZDB486",
       #"Ara-3_25000gen_ZDB488",
       #"Ara-3_27000gen_ZDB309",
       "Ara-3_30000gen_ZDB17",
       #"Ara-3_31500gen_ZDB199",
       #"Ara-3_31500gen_ZDB200",
       "Ara-3_31500gen_ZDB25"
       #"Ara-3_31500gen_ZDB564",
       #"Ara-3_33000gen_CZB154",
       #"Ara-3_50000gen_11364"
)

plots_list <- list()

#test <- read.csv("MOB_1776434_IS150_down_coverage_REL606",sep = '\t')

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
      geom_point(aes(y = high_value), color = "#C8A2C8", size = 2, position = position_dodge(width = 0.2)) +
      geom_point(aes(y = low_value), color = "#C8A2C8", size = 2, position = position_dodge(width = 0.2)) +
      geom_segment(aes(x = position, xend = position, y = high_value, yend = low_value), 
                   color = "#C8A2C8", size = 0.3, linetype = "dashed", position = position_dodge(width = 0.2)) +
      geom_line(color = "darkviolet", size = 1.5, position = position_dodge(width = 0.2)) +
      labs(x = "Position", y = "Coverage (read count)") +
      ylim(0,60) +
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

spacing <- unit(5, "lines")  # Adjust this value as needed

# Create the final collage with adjusted spacing
#collage <- grid.arrange(grobs = plots_list, ncol = 6, top = X, heights = rep(unit(2, "null"), 3), widths = rep(unit(2, "null"), 6), padding = spacing)

collage <- grid.arrange(
  grobs = plots_list,
  ncol = 5,  # Set the number of columns to 5
  heights = rep(unit(2, "null"), length(plots_list)),  # Adjust heights if needed
  widths = rep(unit(2, "null"), 5),  # Set the width of each column
  padding = spacing  # Optional: Adjust padding between plots
)

# Display the collage
print(collage)

ggsave(paste0(X, "_collage1.pdf"), collage, device = "pdf", width = 50, height = 18, units = "in", limitsize = FALSE)
