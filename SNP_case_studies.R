#Two examples were selected manually

egrep "_118_|_731_" 63_coverage_candidates.txt > SNP_case_studies
mapfile -t values < SNP_case_studies
for value in "${values[@]}"; do     X_ancestor=$(grep "$value" ancestor_upregs.bed | cut -f 3);     Y_ancestor=$(grep "$value" ancestor_upregs.bed | cut -f 6);     awk -F '\t' -v X="$X_ancestor" -v Y="$Y_ancestor" '$3>X-200 && $3<X+200 && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X+200}' coverage/*REL60* | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" | sed "s/^/ancestor\t/g" >> "${value}"_coverage_casestudy; done
for value in "${values[@]}"; do     X_evolved=$(grep "$value" evolved_upregs.bed | cut -f 3);     Y_evolved=$(grep "$value" evolved_upregs.bed | cut -f 6);     awk -F '\t' -v X="$X_evolved" -v Y="$Y_evolved" '$3>X-200 && $3<X+200 && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X+200}' coverage/*Ara-2* | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" | sed "s/^/evolved\t/g" >> "${value}"_coverage_casestudy; done
awk -i inplace -F '\t' '{OFS=FS}{$7=400-$7}1' Ara-2_118_SNP_minus_coverage_casestudy
for i in *study; do sed -i "1s/^/status\treplicate\tseqtype\tline\tcount\tposition\tposition_abs\n/g" $i; done

#On to R

library(tidyverse)

getwd()
setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection")

p1 <- read.csv("Ara-2_731_SNP_minus_coverage_casestudy", sep = '\t')

get_average_sd <- function(data) {
  averages <- data %>%
    filter(seqtype == "rna") %>%
    unite("sample",c(status,position),sep="_") %>%
    group_by(sample) %>%
    summarize(sample=sample,
              position_abs=position_abs,
              average = mean(count),
              sd = sd(count),
              n = n())
}

p2 <- get_average_sd(p1)

p2 <- p2 %>% separate(sample, into = c("status","position"),sep="_")

p2 %>%
  #filter(position_abs>100) %>% filter(position_abs<300) %>%
  ggplot(aes(x = position_abs, y = average, color=status)) +
    geom_line() +
    geom_errorbar(aes(ymin = average - sd, ymax = average + sd), width = 0.2) +
    labs(x = "Position", y = "Count") +
    ylim(0, 125) +
    scale_color_manual(values = c("darkviolet", "darkorange1")) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
