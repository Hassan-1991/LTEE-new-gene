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

sed "s/200_commonwindows/sense_antisense/g" htseq_50000gen_200bp_code.sh | sed "s/htseq_50000gen_200bp_count/htseq_sensenantisense/g" | sed "s/all \.\.\//all /g" | sed "s/REL607_sense/REL606_sense/g" | sed -z "s/\n/ \& /g" | sed "s/ \& $//g" > htseq_ALLSENSEANTISENSE.sh


#Htseq
#R stuff for tpm

awk -F '\t' '{OFS=FS}{print $4,$5,$9}' REL606_sense_antisense.gtf | cut -f 1 -d ';' | awk -F '\t' '{OFS=FS}{print $0, $2-$1}' | cut -f 3,4 | cut -f 2- -d ' ' | sed "s/\"//g" > REL606_gene_lengths.txt
join -1 1 -2 4 -t ',' REL606_gene_lengths.txt htseq_test.csv | grep -v "ECB_02723" > htseq_test_lengths.csv

#for some reason there are two ECB_02723s, removed it for now

