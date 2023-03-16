#For the 400bp windows, I want to test these hypotheses:

1 - Amount of annotated:non-annotated stuff goes up/down with time
  a) Sum(TPM)? - Basically, what fraction of the 10^6 are unannotated, can be represented as relative
  b) Mean(TPM)? Average TPM of non-annotated, idk if this is relative to annotated
2 - Number of annotated and/or no-annotated windows goes up/down with time

#Converting htseq count files (marked with gene_types) to tpm values (without averaging for replicates)

test1 <- read.csv("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_50000gen_400bp_count_mark1.csv")
test1$countsum <- test1$est_counts
test1_mod <- test1 %>% unite("replseqline",c(replicate,seqtype,line),sep="_") %>%
  group_by(replseqline) %>%
  summarise(replseqline=replseqline,target_id=target_id,est_counts=est_counts,countsum=sum(countsum/400),type=type,rt_status=rt_status) %>%
  separate(replseqline, into = c("replicate", "seqtype","line"),sep="_")
test1_mod$tpm <- (((test1_mod$est_counts/400)/(test1_mod$countsum))*1000000)

write.csv(test1_mod,"/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_50000gen_400bp_count_mark1_tpm.csv")
sed '/REL60/{;s/$/,ancestor/}' htseq_50000gen_400bp_count_mark1_tpm.csv | sed '/Ara/{;s/$/,evolved/}' | cut -f 2- -d ',' | tr -d "\"" | sed "1s/$/,evostat/g" > htseq_50000gen_400bp_count_mark1_tpm_evostat.csv

test2 <- read.csv("/stor/work/Ochman/hassan/LTEE_analysis/LTEE_data/post_committee_meeting/0620_slidingwindow/step4_csvmarking/htseq_50000gen_400bp_count_mark1_tpm_evostat.csv")

test2$above10PM <- "NO"
test2$above10PM[test2$tpm>10] <- "YES"
test2$above5PM <- "NO"
test2$above5PM[test2$tpm>5] <- "YES"
test2$above1PM <- "NO"
test2$above1PM[test2$tpm>1] <- "YES"

#Convert this to a table - the stupid way

#Number of windows >10 TPM per experiment

for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$11=="YES")' | grep "$i" | grep "rep1" | cut -f 4 -d ',' | sort -u | wc -l; done > test1
for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$11=="YES")' | grep "$i" | grep "rep2" | cut -f 4 -d ',' | sort -u | wc -l; done >> test1

#Just annotated

for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$11=="YES")' | awk -F ',' '($7=="annotated")' | grep "$i" | grep "rep1" | cut -f 4 -d ',' | sort -u | wc -l; done > test2
for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$11=="YES")' | awk -F ',' '($7=="annotated")' | grep "$i" | grep "rep2" | cut -f 4 -d ',' | sort -u | wc -l; done >> test2

#Just unannotated

for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$11=="YES")' | awk -F ',' '($7!="annotated")' | grep "$i" | grep "rep1" | cut -f 4 -d ',' | sort -u | wc -l; done > test3
for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$11=="YES")' | awk -F ',' '($7!="annotated")' | grep "$i" | grep "rep2" | cut -f 4 -d ',' | sort -u | wc -l; done >> test3

#Number of windows >5 TPM per experiment

for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$12=="YES")' | grep "$i" | grep "rep1" | cut -f 4 -d ',' | sort -u | wc -l; done > test4
for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$12=="YES")' | grep "$i" | grep "rep2" | cut -f 4 -d ',' | sort -u | wc -l; done >> test4

#Just annotated

for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$12=="YES")' | awk -F ',' '($7=="annotated")' | grep "$i" | grep "rep1" | cut -f 4 -d ',' | sort -u | wc -l; done > test5
for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$12=="YES")' | awk -F ',' '($7=="annotated")' | grep "$i" | grep "rep2" | cut -f 4 -d ',' | sort -u | wc -l; done >> test5

#Just unannotated

for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$12=="YES")' | awk -F ',' '($7!="annotated")' | grep "$i" | grep "rep1" | cut -f 4 -d ',' | sort -u | wc -l; done > test6
for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$12=="YES")' | awk -F ',' '($7!="annotated")' | grep "$i" | grep "rep2" | cut -f 4 -d ',' | sort -u | wc -l; done >> test6

#Number of windows >1 TPM per experiment

for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$13=="YES")' | grep "$i" | grep "rep1" | cut -f 4 -d ',' | sort -u | wc -l; done > test7
for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$13=="YES")' | grep "$i" | grep "rep2" | cut -f 4 -d ',' | sort -u | wc -l; done >> test7

#Just annotated

for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$13=="YES")' | awk -F ',' '($7=="annotated")' | grep "$i" | grep "rep1" | cut -f 4 -d ',' | sort -u | wc -l; done > test8
for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$13=="YES")' | awk -F ',' '($7=="annotated")' | grep "$i" | grep "rep2" | cut -f 4 -d ',' | sort -u | wc -l; done >> test8

#Just unannotated

for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$13=="YES")' | awk -F ',' '($7!="annotated")' | grep "$i" | grep "rep1" | cut -f 4 -d ',' | sort -u | wc -l; done > test9
for i in REL606 REL607 AraP1 AraP2 AraP3 AraP4 AraP5 AraM1 AraM2 AraM3 AraM4 AraM5 AraM6; do cut -f 2- -d ',' htseq_50000gen_400bp_count_mark1_tpm_evostat_cutoffs.csv | tr -d "\"" | awk -F ',' '($2=="rna"&&$13=="YES")' | awk -F ',' '($7!="annotated")' | grep "$i" | grep "rep2" | cut -f 4 -d ',' | sort -u | wc -l; done >> test9

printf 'REL606\nREL607\nAraP1\nAraP2\nAraP3\nAraP4\nAraP5\nAraM1\nAraM2\nAraM3\nAraM4\nAraM5\nAraM6\n%.0s' {1..2} > test10

printf '18801,10006,8796\n%.0s' {1..26} > test11

printf 'rep1\n%.0s' {1..13} > test12
printf 'rep2\n%.0s' {1..13} >> test12

paste test10 test11 test12 test1 test2 test3 test4 test5 test6 test7 test8 test9 | sed "s/\t/,/g" | sed '/REL60/{;s/$/,ancestor/}' | sed '/Ara/{;s/$/,evolved/}' | sed -z "1s/^/line,replicate,above10_all,above10_annotated,above10_unannotated,above5_all,above5_annotated,above5_unannotated,above1_all,above1_annotated,above1_unannotated,totalwindows,totalannotatedwindows,totanunannotatedwindows,evostat\n/g" > htseq_50000gen_400bp_cutofftables.csv

#Add in labels (too lazy, should be ez with sed)

#Is there a significant difference between the number of transcribed windows between anc and evo?

ggplot(test3, aes(x=evostat,y=above5_all)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test3, aes(x=evostat,y=above5_annotated)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test3, aes(x=evostat,y=above5_unannotated)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test3, aes(x=evostat,y=above1_all)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test3, aes(x=evostat,y=above1_annotated)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test3, aes(x=evostat,y=above1_unannotated)) + geom_boxplot() + stat_compare_means(method = "t.test")

#Particular lines:

test4 <- test3 %>% filter(line=="REL60"|line=="REL607"|line=="AraP1")

ggplot(test4, aes(x=evostat,y=above5_all)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test4, aes(x=evostat,y=above5_annotated)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test4, aes(x=evostat,y=above5_unannotated)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test4, aes(x=evostat,y=above1_all)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test4, aes(x=evostat,y=above1_annotated)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test4, aes(x=evostat,y=above1_unannotated)) + geom_boxplot() + stat_compare_means(method = "t.test")

#Is there a difference between the tpm values between ancestor and evolved?

ggplot(test1 %>% filter(type=="annotated"), aes(x=evostat,y=tpm)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test1 %>% filter(type=="unannotated"), aes(x=evostat,y=tpm)) + geom_boxplot() + stat_compare_means(method = "t.test")

#Between ancestor and each line?

ggplot(test1 %>% filter(type=="annotated") %>% filter(line=="REL606"|line=="REL607"|line=="AraP5"), aes(x=evostat,y=tpm)) + geom_boxplot() + stat_compare_means(method = "t.test")
ggplot(test1 %>% filter(type!="annotated") %>% filter(line=="REL606"|line=="REL607"|line=="AraP5"), aes(x=evostat,y=tpm)) + geom_boxplot() + stat_compare_means(method = "t.test")

