 #Goal: Write scripts to
#1 - Get relative fraction of antisense transcription to annotated
#  (The window approach gave me number of windows transcribed, but this is total amount - would be complementary)
#2 - Relative fraction of transcription in gene-antigene pairs

#let's keep it simple and focus on protein-coding genes.

#Step-1: Get coordinates of all antisense elements. Just flip the strands

#We'll just have two gtfs instead of one, but I think it's imp to make sure there are no duplicates.
#Also start with insremoved versions

cp ../../insremoved_gtf/REL606_insremoved.gtf .

awk -F '\t' '($3=="CDS")' ../../insremoved_gtf/REL606_insremoved.gtf | sed "s/\t+\t/\t%\t/g" | sed "s/\t-\t/\t+\t/g" | sed "s/\t%\t/\t-\t/g" | sed "s/id \"/id \"antisense_/g" > REL606_insremoved_antisense.gtf

#Do any of these antisense regions encroach on other gene territories?

#I need to do subset analyses with 30bp, 50bp, 100bp removal to see if results change.

#let's start with 30bp, just REL606+REL607

bedtools closest -a REL606_insremoved_antisense.gtf -b REL606_insremoved.gtf -s -d | awk -F '\t' '($19<30)' | cut -f 9 | sort -u | cut -f 2 -d "\"" | sort -u | cut -f 2- -d "_" | sort -u > REL606_tooclose_antisense.txt
cat REL606_insremoved.gtf REL606_insremoved_antisense.gtf | grep -v -f REL606_tooclose_antisense.txt > REL606_sense_antisense.gtf

#GENERALIZE

awk -F '\t' '($3=="CDS")' something_insremoved.gtf | sed "s/\t+\t/\t%\t/g" | sed "s/\t-\t/\t+\t/g" | sed "s/\t%\t/\t-\t/g" | sed "s/id \"/id \"antisense_/g" > something_insremoved_antisense.gtf
bedtools closest -a something_insremoved_antisense.gtf -b something_insremoved.gtf -s -d | awk -F '\t' '($19<30)' | cut -f 9 | sort -u | cut -f 2 -d "\"" | sort -u | cut -f 2- -d "_" | sort -u > something_tooclose_antisense.txt
cat something_insremoved.gtf something_insremoved_antisense.gtf | grep -v -f something_tooclose_antisense.txt > something_sense_antisense.gtf

sed "s/^/sed \"s\/something\//g" somethings | sed "s/$/\/g\" basic_code.sh >> final_code.sh/g" > interim_code.sh 


#I wrote the HTseq code by changing the file names of pre-existing htseq scripts, generated elsewhere in this repo

sed "s/200_commonwindows/sense_antisense/g" htseq_50000gen_200bp_code.sh | sed "s/htseq_50000gen_200bp_count/htseq_sensenantisense/g" | sed "s/all \.\.\//all /g" | sed "s/REL607_sense/REL606_sense/g" | sed -z "s/\n/ \& /g" | sed "s/ \& $//g" > htseq_ALLSENSEANTISENSE.sh

#To get lengths for each gene per line, I first fused gene names with line names:

sed "s/,ECB/_ECB/g" htseq_sensenantisense.csv | sed "s/,antisense/_antisense/g" > htseq_sensenantisense_interim.csv

#I then ran the code for getting lengths:

awk -F '\t' '{OFS=FS}{print $4,$5,$9}' something*sense_antisense.gtf | cut -f 1 -d ';' | awk -F '\t' '{OFS=FS}{print $0, $2-$1}' | cut -f 3,4 | cut -f 2- -d ' ' | sed "s/\"//g" | sed "s/^/something_/g" > something_gene_lengths.txt

#The "somethings" are Ara+1,Ara+2 etc; did further manual edits in the last version of the code
#sort both files and join

#final cleanup

sed "s/_ECB/,ECB/g" htseq_sensenantisense_lengths.csv | sed "s/_antisense/,antisense/g" | sed "s/antisense,/antisense_/g" > test && mv test htseq_sensenantisense_lengths.csv

#On to R!!!

#Rscript

library(tidyverse)
library(dplyr)

test1 <- read.csv("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0306_exploratory_antisense_intergenic_readthrough/htseq_sensenantisense_lengths.csv")
test1$countsum <- test1$est_counts
test1 <- test1 %>% unite("replseqline",c(replicate,seqtype,line),sep="_") %>%
  group_by(replseqline) %>%
  summarise(replseqline=replseqline,length=length,target_id=target_id,est_counts=est_counts,countsum=sum(countsum/length)) %>%
  separate(replseqline, into = c("replicate", "seqtype","line"),sep="_")

test1$tpm <- (((test1$est_counts/test1$length)/(test1$countsum))*1000000)

test2 <- test1 %>%
  unite("seqline",c(seqtype,line),sep="_") %>%
  unite("sample",c(seqline,target_id),sep="%") %>%
  dplyr::select(sample,replicate,tpm) %>%
  group_by(sample) %>%
  summarise(sample=sample,tpm=mean(tpm)) %>%
  unique() %>%
  separate(sample,into=c("seqline","target_id"),sep="%") %>%
  separate(seqline,into=c("seqtype","line"),sep="_")

write.csv(test2,"/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0306_exploratory_antisense_intergenic_readthrough/htseq_senseantisense_tpm.csv")

#marking

cut -f 2- -d ',' htseq_senseantisense_tpm.csv | tr -d "\"" | grep "antisense" | sed "s/antisense_//g" | sed "s/$/,antisense/g" > htseq_senseantisense_tpm_status.csv
cut -f 2- -d ',' htseq_senseantisense_tpm.csv | tr -d "\"" | tail -n+2 | grep -v "antisense" | sed "s/$/,sense/g" >> htseq_senseantisense_tpm_status.csv

#Marking with evo/anc status:

sed '/REL60/{;s/$/,ancestor/}' htseq_senseantisense_tpm_status.csv | sed '/Ara/{;s/$/,evolved/}' > htseq_senseantisense_tpm_status_evostat.csv
