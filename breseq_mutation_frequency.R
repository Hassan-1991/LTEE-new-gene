#I extracted the generation - allele freq pairs manually with this code:

#For SNP or DEL:

grep "1881296" Ara-2*K/output/*gd | grep "SNP.*frequency=" | sed "s/K\//_/g" | sort -nk 2 -t "_" | cut -f 1,7 | sed "s/\t/_/g" | sed "s/=/_/g" | cut -f 2,5 -d "_" | sed "s/_/\t/g" | awk -F'\t' 'FNR == NR { data[$1] = $0; next } $1 in data { print data[$1] } !($1 in data) { print $1 "\t0" }' - 0_to_60

#For MOB:

grep "1881296" Ara-2*K/output/*gd | grep "MOB.*frequency=" | sed "s/K\//_/g" | sort -nk 2 -t "_" | cut -f 1,9 | sed "s/\t/_/g" | sed "s/=/_/g" | cut -f 2,5 -d "_" | sed "s/_/\t/g" | awk -F'\t' 'FNR == NR { data[$1] = $0; next } $1 in data { print data[$1] } !($1 in data) { print $1 "\t0" }' - 0_to_60

#Then pasted them in excel.

#R code (WIP)

library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)  # Load the gridExtra package for grid.arrange()

setwd("/stor/scratch/Ochman/hassan/Good2017_LTEE_data/trimmed")

p1 <- read.csv("breseq_10protogenes_curated.tsv", sep = '\t')

values_to_plot <- c("Ara-2_118_SNP","Ara-2_610_SNP","Ara-2_731_SNP","Ara-3_547_MOB","Ara-6_24_MOB","Ara-6_59_DEL","Ara+1_94_DEL","Ara+1_99_MOB","Ara+5_30_MOB","Ara+5_65_MOB")

plots_list <- list()  # Create an empty list to store the plots

for (value in values_to_plot) {
  plot <- p1 %>%
    filter(target_id == value) %>%
    ggplot(aes(x = Generation, y = Allele_freq / 100)) +
    geom_line(size = 1.5) +
    scale_x_continuous(breaks = seq(0, 60000, 5000)) +
    labs(title = value, x = "Generation", y = "Allele frequency") +
    ylim(0, 1) +
    theme(
      plot.background = element_blank(),  # Transparent plot background
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.background = element_rect(fill = 'transparent'),
      legend.box.background = element_rect(fill = 'transparent'),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 10, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 16),
      legend.key.size = unit(1, 'cm'),
      legend.text = element_text(size = 15),
      legend.position = "top",
      plot.margin = margin(1, 1, 1, 1),  # Adjust the margins around the plot
      panel.background = element_rect(color = "black", size = 1)  # Add a black border around the plot area
    )
  
  plots_list[[value]] <- plot  # Add the plot to the list
}

# Create a collage by arranging all the plots in a grid
collage <- do.call(grid.arrange, c(plots_list, ncol = 1))

# Save the collage
ggsave("collage.pdf", collage, device="pdf", width = 20, height = 60, units = "in", limitsize=FALSE)
