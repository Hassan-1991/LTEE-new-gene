#This code describes how to make the transcriptome comprised of 400-bp windows
#And mark the windows with their annotation status
#The final file is fed into the code of figure 2

#Step-1: Make an initial file with all 400-bp windows in the genome

awk '!/^>/ { printf "%s", $0; n="\n"} /^>/ {print n $0; n = "" } END {printf "%s", n }' REL606.fasta | grep -v "^>" | wc

seq 1 400 4629813 > col1 #Start coordinate of each window
seq 1 400 4629813 | tail -n+2 > col2 #Stop coordinate of each window
paste col1 col2 | cat -n | sed "s/ //g" | awk -F '\t' '{OFS=FS}{print "REL606\t.\tCDS",$2,$3-1,".\t+\t0\tgene_id \"window_plus_",$1,"\";transcript_id \"window_plus_",$1,"\";"}' | head -n-1 | sed "s/_\t/_/g" | sed "s/\t\"/\"/g" > REL606_slidwindow400.gtf #Windows on plus strand
paste col1 col2 | cat -n | sed "s/ //g" | awk -F '\t' '{OFS=FS}{print "REL606\t.\tCDS",$2,$3-1,".\t-\t0\tgene_id \"window_minus_",$1,"\";transcript_id \"window_minus_",$1,"\";"}' | head -n-1 | sed "s/_\t/_/g" | sed "s/\t\"/\"/g" >> REL606_slidwindow400.gtf #Windows on minus strand
rm col* #remove intermediate files
bedtools intersect -a REL606_slidwindow400.gtf -b REL606.gtf -wo | egrep "repeat_region|IS" | cut -f9 | sort -u | grep -vf - REL606_slidwindow400.gtf > test && mv test REL606_slidwindow400.gtf
gffread -E -w REL606_slidwindow400_CDS.faa -g REL606.fasta REL606_slidwindow400.gtf #get their sequences

#Step-2: From this file, get rid of all windows that map more than once in any genome, experienced insertions or deletions, diverged substantially, or got deleted from any evolved genome

#We first map the window sequences to each evolved genome with gmap

variables=("Ara+1" "Ara+2" "Ara+3" "Ara+4" "Ara+5" "Ara-1" "Ara-2" "Ara-3" "Ara-4" "Ara-5" "Ara-6")
for var in "${variables[@]}"; do
echo "/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d $var -f 2 --gff3-fasta-annotation=1 REL606_slidwindow400_CDS.faa > ${var}_slidwindow400_CDS.gff3"
echo "grep \"mRNA\" ${var}_slidwindow400_CDS.gff3 | grep \"path2\" | cut -f 2 -d '=' | cut -f 1 -d '.' | sed 's/^/=/g' | sed 's/\$/;/g' > ${var}_400_doublemaps"
echo "grep -v -f ${var}_400_doublemaps ${var}_slidwindow400_CDS.gff3 | grep \"mRNA\" > ${var}_400_step1"
echo "awk -F '\t' '(\$5-\$4>410)||(\$5-\$4<390)' ${var}_400_step1 > ${var}_400_longshortmaps"
echo "grep -v -f ${var}_400_longshortmaps ${var}_400_step1 > ${var}_400_step2"; done | split -l 5 -

#Run the resulting x-files in separate screens

rm *longshortmaps
rm *doublemaps #Remove interim files

#check to see if the remaining "common sequences" map with a suspiciously low identity or query cover by scanning their ID and qCov values: 

cat *step2 | grep "mRNA" | grep "path1" | sed "s/=/;/g" | cut -f 12-15 -d ';' | sort -u | sort -t ';' -nk 2 | less
#lowest coverage - 98 for 400
#lowest identity - 78.1

#These values are not very low, so we don't need to subtract anything.

#From this set, pick out the windows that appear 11 times - meaning they haven't been deleted from any of the evolved genomes

cat Ara*400*step2 | cut -f 2,9 | cut -f 1 -d ';' | sed "s/.mrna1//g" | sed "s/ID=//g" | cut -f 2 | sort | uniq -c | sed "s/ //g" | grep "11window" | sed "s/11window/window/g" | sed "s/^/=/g" | sed "s/$/;/g" > common_400bp_windows.txt

#Extract these windows from the gff3s constructed earlier, convert them into gtfs

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6;
do grep -f common_400bp_windows.txt "$i"_400_step2 | sed "s/ID=/transcript_id=\"/g" | sed "s/.mrna1;Name=/\";gene_id \"/g" | sed "s/;Parent/\";Parent/g" | sed "s/\tmRNA\t/\tCDS\t/g" | awk -F '\t' '{OFS=FS}{print $1,$2=".",$3,$4,$5,$6,$7,$8,$9}' > "$i"_400_commonwindows.gtf; done

#Generate the feature lists for two ancestors

sed "s/=/\"/g" common_400bp_windows.txt | sed "s/;/\"/g" > common_400bp_windows_2.txt 
grep -f common_400bp_windows_2.txt REL606_slidwindow400.gtf > REL606_400_commonwindows.gtf 
cp REL606_400_commonwindows.gtf REL607_400_commonwindows.gtf

#BOWTIE AND HTSEQ STEPS#

#htseq generates one file for each experiment, cat them into a .csv file. Then:

sed "s/$/,400/g" htseq_50000gen_400bp_count.csv | sed "1s/^/replicate,seqtype,line,target_id,est_counts,length\n/g" > htseq_50000gen_400bp_count_lengths.csv

#Now we analyze this with R code for generating TPM counts:

library(ggplot2)
library(dplyr)
library(tidyverse)

setwd("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0622_post_rejection/100bp")

variable_names <- c("Ara+1", "Ara+2", "Ara+3", "Ara+4", "Ara+5", "Ara-1", "Ara-2", "Ara-3", "Ara-4", "Ara-5", "Ara-6")

for (variable in variable_names) {
  
  p1 <- read_tsv(paste0(variable, "_htseq_lengths.tsv"))
  p1 <- p1 %>% unite("replseqline",c(replicate,seqtype,line),sep="_") %>%
    group_by(replseqline) %>%
    summarise(replseqline=replseqline,target_id=target_id,count=count,countsum=sum(count/length),status=status,length=length) %>%
    separate(replseqline, into = c("replicate", "seqtype","line"),sep="_")
  p1$tpm <- (((p1$count/p1$length)/(p1$countsum))*1000000)
  p1 <- p1 %>%
    group_by(target_id, status) %>%
    summarise(average_tpm = mean(tpm))
  write_csv(p1, paste0(variable, "_htseq_tpm.tsv"))
}

#Back to bash:

sed -i "s/\"//g" htseq_50000gen_400bp_count_lengths_tpm.csv

#Marking as annotated, readthrough, non-coding, antisense...

bedtools intersect -s -a REL606_400_commonwindows.gtf -b REL606.gtf -wo | awk -F '\t' '($19>10)' | cut -f 9 | cut -f 2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" > annot
grep -vf annot REL606_400_commonwindows.gtf | bedtools intersect -S -wo -a - -b REL606.gtf  | awk -F '\t' '($19>10)' | cut -f9 | cut -f 2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" > antisense
cat annot antisense | grep -vf - REL606_400_commonwindows.gtf | cut -f9 | cut -f 2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" > intergenic
cat antisense intergenic > non_coding
grep -f non_coding REL606_400_commonwindows.gtf | bedtools closest -s -D b -a - -b REL606.gtf | awk -F '\t' '(($19>-300)&&($19<0))' | cut -f 9 | cut -f 2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" > upstream_300
grep -f non_coding REL606_400_commonwindows.gtf | bedtools closest -s -D b -a - -b REL606.gtf | awk -F '\t' '(($19>-99)&&($19<0))' | cut -f 9 | cut -f 2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" > upstream_99
grep -f non_coding REL606_400_commonwindows.gtf | bedtools closest -s -D b -a - -b REL606.gtf | awk -F '\t' '(($19>0)&&($19<99))' | cut -f 9 | cut -f 2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" > readthrough

for i in annot antisense intergenic non_coding upstream_300 upstream_99 readthrough; do sed -i "s/\"/,/g" $i; done

grep -vf upstream_99 upstream_300 > upstream_300_only
grep -vf upstream_300 non_coding > non_upstream
grep -vf readthrough non_coding > non_readthrough

grep -f annot htseq_50000gen_400bp_count_lengths_tpm.csv | sed "s/$/,annot/g" > test
grep -f antisense htseq_50000gen_400bp_count_lengths_tpm.csv | sed "s/$/,antisense/g" >> test
grep -f intergenic htseq_50000gen_400bp_count_lengths_tpm.csv | sed "s/$/,intergenic/g" >> test

grep -f annot test | sed "s/$/,annot/g" > test2
grep -f upstream_300_only test | sed "s/$/,upstream_300/g" >> test2
grep -f upstream_99 test | sed "s/$/,upstream_99/g" >> test2
grep -f non_upstream test | sed "s/$/,nonreadthrough/g" >> test2

grep -f annot test2 | sed "s/$/,annot/g" > htseq_50000gen_400bp_count_lengths_tpm_marked.csv
grep -f readthrough test2 | sed "s/$/,readthrough/g" >> htseq_50000gen_400bp_count_lengths_tpm_marked.csv
grep -f non_readthrough test2 | sed "s/$/,nonreadthrough/g" >> htseq_50000gen_400bp_count_lengths_tpm_marked.csv

sed -i "1s/^/replicate,seqtype,line,target_id,est_counts,countsum,length,tpm,annot_status,upstream_status,readthrough_status\n/g" htseq_50000gen_400bp_count_lengths_tpm_marked.csv

#Get rid of the few cases that are difficult to classify - ones that do overlap with annotated features but with <10bp of overlap

grep -f non_coding REL606_400_commonwindows.gtf | bedtools closest -s -D b -a - -b REL606.gtf | awk -F '\'t '($19==0)' | cut -f 9 | cut -f 2 -d "\"" | sort -u | sed "s/^/,/g" | sed "s/$/,/g" > difficult_to_classify
grep -vf difficult_to_classify htseq_50000gen_400bp_count_lengths_tpm_marked.csv > htseq_50000gen_400bp_count_lengths_tpm_marked_final.csv
