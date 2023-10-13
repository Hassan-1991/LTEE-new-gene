#50K dataset

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do grep -v "#" "$i"_50000gen_[0-9][0-9][0-9][0-9][0-9].gd | sed "s/^/"$i"\t/g" | awk -F '\t' '{OFS=""}{print $1,"_",$3,"_",$2,"\t",$0}'; done | cut -f 1,6- > 50k_ancestorcoords.gd
for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do grep -v "#" "$i"_50000gen*applied.gd | sed "s/^/"$i"\t/g" | awk -F '\t' '{OFS=""}{print $1,"_",$3,"_",$2,"\t",$0}'; done > 50k_applied.gd

#Making the coordinates easier to extract by putting them in bed files
#The easiest way to do this is specify mutation start and end coordinates on both strands, and then just retain the proto-gene-relevant strand

grep -f 63_coverage_candidates_withoutstrand.txt 50k_ancestorcoords.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | awk -F '\t' '{OFS=""}{print $1,"\t",$4,"\t",$4,"\t",$2,"_plus\t",0,"\t+"}' | grep -v "DEL" > ancestor_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_ancestorcoords.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | awk -F '\t' '{OFS=""}{print $1,"\t",$4,"\t",$4,"\t",$2,"_minus\t",0,"\t-"}' | grep -v "DEL" >> ancestor_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_ancestorcoords.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | grep "DEL" | awk -F '\t' '{OFS=""}{print $1,"\t",$4,"\t",$4+$5,"\t",$2,"_plus\t",0,"\t+"}' >> ancestor_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_ancestorcoords.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | grep "DEL" | awk -F '\t' '{OFS=""}{print $1,"\t",$4,"\t",$4+$5,"\t",$2,"_minus\t",0,"\t-"}' >> ancestor_upregs
grep -f 63_coverage_candidates.txt ancestor_upregs > ancestor_upregs.bed

grep -f 63_coverage_candidates_withoutstrand.txt 50k_applied.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | egrep "SNP|DEL" | awk -F '\t' '{OFS=""}{print $7,"\t",$8,"\t",$8,"\t",$2,"_plus\t0\t+"}' > evolved_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_applied.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | egrep "SNP|DEL" | awk -F '\t' '{OFS=""}{print $7,"\t",$8,"\t",$8,"\t",$2,"_minus\t0\t-"}' >> evolved_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_applied.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | egrep -v "SNP|DEL" | awk '{OFS="\t"}{ for (i=1; i<=NF; i++) { if (i == 1 || i == 2 || $i ~ /(applied_start=|applied_end=)/) printf "%s ", $i } printf "\n" }' | sed "s/applied_end=//g" | sed "s/applied_start=//g" | awk '{OFS=""}{print $1,"\t",$4,"\t",$3,"\t",$2,"_plus\t0\t+"}' >> evolved_upregs
grep -f 63_coverage_candidates_withoutstrand.txt 50k_applied.gd | awk -F "_" '{OFS=""}{print $1,"\t",$0}' | egrep -v "SNP|DEL" | awk '{OFS="\t"}{ for (i=1; i<=NF; i++) { if (i == 1 || i == 2 || $i ~ /(applied_start=|applied_end=)/) printf "%s ", $i } printf "\n" }' | sed "s/applied_end=//g" | sed "s/applied_start=//g" | awk '{OFS=""}{print $1,"\t",$4,"\t",$3,"\t",$2,"_minus\t0\t-"}' >> evolved_upregs
grep -f 63_coverage_candidates.txt evolved_upregs > evolved_upregs.bed

#Make genome coverage files for each RNA and Riboseq run

for i in $(ls /stor/work/Ochman/hassan/LTEE_analysis/50k_bowtie_indices_bamfiles/*sorted.bam | rev | cut -f 1 -d '/' | cut -f 2- -d "_" | rev | sort -u);
do
bedtools genomecov -d -strand + -ibam /stor/work/Ochman/hassan/LTEE_analysis/50k_bowtie_indices_bamfiles/"$i"*sorted.bam | sed "s/^/"$i"\t/g" | sed "s/$/\t+/g" | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" > "$i"_plus
bedtools genomecov -d -strand - -ibam /stor/work/Ochman/hassan/LTEE_analysis/50k_bowtie_indices_bamfiles/"$i"*sorted.bam | sed "s/^/"$i"\t/g" | sed "s/$/\t-/g" | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" > "$i"_minus
cat "$i"_plus "$i"_minus > "$i"_coverage.tsv
done

#Make coverage plots for each mutation

mapfile -t values < 63_coverage_candidates.txt

#ancestor_right_200

for value in "${values[@]}"
do
    X_ancestor=$(grep "$value" ancestor_upregs.bed | cut -f 3)
    Y_ancestor=$(grep "$value" ancestor_upregs.bed | cut -f 6)
    awk -F '\t' -v X="$X_ancestor" -v Y="$Y_ancestor" '$3>X && $3<X+200 && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X}' coverage/*REL60* | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" >> "${value}_coverage_ancestor"_right
done

#ancestor_left_200

for value in "${values[@]}"
do
    X_ancestor=$(grep "$value" ancestor_upregs.bed | cut -f 2)
    Y_ancestor=$(grep "$value" ancestor_upregs.bed | cut -f 6)
    awk -F '\t' -v X="$X_ancestor" -v Y="$Y_ancestor" '$3>X-200 && $3<X && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X+200}' coverage/*REL60* | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" >> "${value}_coverage_ancestor_left"
done

#evolved_right_200

for value in "${values[@]}"
do
    X_evolved=$(grep "$value" evolved_upregs.bed | cut -f 3)
    Y_evolved=$(grep "$value" evolved_upregs.bed | cut -f 6)
    awk -F '\t' -v X="$X_evolved" -v Y="$Y_evolved" '$3>X && $3<X+200 && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X}' coverage/*"$(cut -f 1 -d "_" <<< "$value")"* | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" >> "${value}_coverage_right"
done

#evolved_Left_200

for value in "${values[@]}"
do
    X_evolved=$(grep "$value" evolved_upregs.bed | cut -f 2)
    Y_evolved=$(grep "$value" evolved_upregs.bed | cut -f 6)
    awk -F '\t' -v X="$X_evolved" -v Y="$Y_evolved" '$3>X-200 && $3<X && $5==Y { print $1 "\t" $4 "\t" $3 "\t" $3-X+200}' coverage/*"$(cut -f 1 -d "_" <<< "$value")"* | sed "s/-rna-/\trna\t/g" | sed "s/-ribo-/\tribo\t/g" >> "${value}_coverage_left"
done

#For coverage plots in the minus strand, one more fix is required

for i in $(grep "minus" ../63_coverage_candidates.txt | cut -f 1-3 -d "_"); do awk -i inplace '{OFS='\t'}{$6 = 200 - $6}1' "$i"*coverage*; done
sed -i "1s/^/replicate\tseqtype\tline\tcount\tposition_abs\tposition\n/g" *coverage

#R code to generate coverage plots

library(dplyr)
library(ggplot2)
library(gridExtra)
library(stringr)
library(grid)

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/all_coverages")

function(data) {
  averages <- data %>%
    filter(seqtype == "rna") %>%
    group_by(position) %>%
    summarize(average = mean(count),
              sd = sd(count),
              n = n())
}

values <- c(
  "Ara+1_31_MOB",
  "Ara+1_35_MOB",
  "Ara+1_102_MOB",
  "Ara+5_30_MOB",
  "Ara-5_84_MOB",
  "Ara-6_24_MOB")

# Create a function to generate the plots
generate_plots <- function(X) {
  p1 <- read.csv(paste0(X, "_coverage_ancestor_left"), sep = '\t')
  p2 <- read.csv(paste0(X, "_coverage_ancestor_right"), sep = '\t')
  p3 <- read.csv(paste0(X, "_coverage_left"), sep = '\t')
  p4 <- read.csv(paste0(X, "_coverage_right"), sep = '\t')
  
  plot1 <- ggplot(get_average_sd(p1), aes(x = position, y = average)) +
    geom_line(color = "darkviolet") +
    geom_errorbar(aes(ymin = average - sd, ymax = average + sd), width = 0.2, color = "darkviolet") +
    labs(x = "Position", y = "Coverage (read count)", title ="Ancestor, upstream") +
    #ylim(0, 15) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size=20),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))
  
  plot2 <- ggplot(get_average_sd(p2), aes(x = position, y = average)) +
    geom_line(color = "darkviolet") +
    geom_errorbar(aes(ymin = average - sd, ymax = average + sd), width = 0.2, color = "darkviolet") +
    labs(x = "Position", y = "Coverage (read count)", title ="Ancestor, downstream") +
    #ylim(0, 15) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size=20),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))
  
  plot3 <- ggplot(get_average_sd(p3), aes(x = position, y = average)) +
    geom_line(color = "darkorange1") +
    geom_errorbar(aes(ymin = average - sd, ymax = average + sd), width = 0.2, color = "darkorange1") +
    labs(x = "Position", y = "Coverage (read count)", title ="Evolved, upstream") +
    #ylim(0, 15) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size=20),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))
  
  plot4 <- ggplot(get_average_sd(p4), aes(x = position, y = average)) +
    geom_line(color = "darkorange1") +
    geom_errorbar(aes(ymin = average - sd, ymax = average + sd), width = 0.2, color = "darkorange1") +
    labs(x = "Position", y = "Coverage (read count)", title ="Evolved, downstream") +
    #ylim(0, 15) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size=20),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16))
  
  collage <- grid.arrange(plot1, plot2, plot3, plot4, ncol = 2, top=textGrob(X, gp=gpar(fontsize=25)))
  return(collage)
}

plot_list <- lapply(values, generate_plots)
final_plot <- do.call(grid.arrange, c(plot_list, ncol = 2))
print(final_plot)

ggsave(paste0("firstsix.pdf"), final_plot, width = 22, height = 30, units = "in")
