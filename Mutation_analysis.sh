####Step-1: Make individualized GTF files for each generation, containing 500bp plus_down and minus_up####

ls *.gd > filenames
sed "s/^/cat /g" filenames | sed "s/$/ \| grep \-v \"\#\" \| cut \-f 1,5,6,7 \> /g" | 
awk -F ' ' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$2}' | 
sed "s/$/t/g" | 
sed "s/.gdt/_mutations.tsv/g" |
bash

#nano the following (marked with double ##s) into a file called basic_code.sh:

##

#MOB:

sed "s/\t/,/g" SOMETHING_mutations.tsv | grep "MOB" | sed "s/,/_/g" | sed "s/^/transcript_id \"SOMETHING_/g" | sed "s/$/\";gene_id \"SOMETHING_/g" | awk -F "\"" '{OFS=FS}{print $1,$2,$3,$2}' | sed "s/$/\";/g" | sort -t "_" -nk 3 > col9_MOB
sed "s/\t/,/g" SOMETHING_mutations.tsv | grep "MOB" | sed "s/^/REL606\t.\tCDS\t/g" | grep -v "\-1" | sed "s/MOB,//g" | sed "s/,/\t/g" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$4+500,".\t+\t0"}' > col1_8_MOBplus
sed "s/\t/,/g" SOMETHING_mutations.tsv | grep "MOB" | sed "s/^/REL606\t.\tCDS\t/g" | grep "\-1" | sed "s/MOB,//g" | sed "s/,/\t/g" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$4-500,".\t-\t0"}' | sed "s/,/\t/g" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8}' > col1_8_MOBminus
cat col1_8_MOBplus col1_8_MOBminus | sort -t '	' -nk 4 > col1_8_MOB
paste col1_8_MOB col9_MOB > all_MOBs.gtf
grep "tabby-tabby" all_MOBs.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$4,$6,$7,$8,$9}' > all_MOBs2.gtf
grep "tabby+tabby" all_MOBs.gtf >> all_MOBs2.gtf
sort -t 'tabby' -nk 4 all_MOBs2.gtf > all_MOBs.gtf

#SNP:

grep "SNP" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/,/_/g" | sed "s/^/transcript_id \"SOMETHING_/g" | sed "s/$/\";gene_id \"SOMETHING_/g" | awk -F "\"" '{OFS=FS}{print $1,$2,$3,$2}' | sed "s/$/\";/g" > col9_SNP 
grep "SNP" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/^/REL606\t.\tCDS\t/g" | awk -F ',' '{OFS=FS}{print $1,$2,$2+500}' | sed "s/SNP,//g" | sed "s/,/\t/g" | sed "s/$/\t.\t+\t0/g" > col1-8_SNP_plus_down 
grep "SNP" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/^/REL606\t.\tCDS\t/g" | awk -F ',' '{OFS=FS}{print $1,$3=$2,$2=$2-500}' | awk -F ',' '{OFS=FS}{print $1,$3,$2}' | sed "s/SNP,//g" | sed "s/,/\t/g" | sed "s/$/\t.\t-\t0/g" > col1-8_SNP_minus_up 
paste col1-8_SNP_plus_down col9_SNP | sed "s/\";/_plus_down\";/g" > all_SNP_regions.gtf 
paste col1-8_SNP_minus_up col9_SNP | sed "s/\";/_minus_up\";/g" >> all_SNP_regions.gtf 
sort -t 'tabby' -nk 4 all_SNP_regions.gtf > all_SNP_regions2.gtf 
mv all_SNP_regions2.gtf all_SNP_regions.gtf

#INS

grep "INS" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/,/_/g" | sed "s/^/transcript_id \"SOMETHING_/g" | sed "s/$/\";gene_id \"SOMETHING_/g" | awk -F "\"" '{OFS=FS}{print $1,$2,$3,$2}' | sed "s/$/\";/g" > col9_INS  
grep "INS" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/^/REL606\t.\tCDS\t/g" | awk -F ',' '{OFS=FS}{print $1,$2,$2+500}' | sed "s/INS,//g" | sed "s/,/\t/g" | sed "s/$/\t.\t+\t0/g" > col1-8_INS_plus_down  
grep "INS" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/^/REL606\t.\tCDS\t/g" | awk -F ',' '{OFS=FS}{print $1,$3=$2,$2=$2-500}' | awk -F ',' '{OFS=FS}{print $1,$3,$2}' | sed "s/INS,//g" | sed "s/,/\t/g" | sed "s/$/\t.\t-\t0/g" > col1-8_INS_minus_up  
paste col1-8_INS_plus_down col9_INS | sed "s/\";/_plus_down\";/g" > all_INS_regions.gtf  
paste col1-8_INS_minus_up col9_INS | sed "s/\";/_minus_up\";/g" >> all_INS_regions.gtf  
sort -t 'tabby' -nk 4 all_INS_regions.gtf > all_INS_regions2.gtf  
mv all_INS_regions2.gtf all_INS_regions.gtf

#DEL

grep "DEL" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/,/_/g" | sed "s/^/transcript_id \"SOMETHING_/g" | sed "s/$/\";gene_id \"SOMETHING_/g" | awk -F "\"" '{OFS=FS}{print $1,$2,$3,$2}' | sed "s/$/\";/g" > col9_DEL 
grep "DEL" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/^/REL606\t.\tCDS\t/g" | awk -F ',' '{OFS=FS}{print $1,$2,$2+$3}' | sed "s/DEL,//g" | sed "s/,/\t/g" | sed "s/$/\t.\t+\t0/g" > col1-8_DEL 
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+500,$6,$7,$8}' col1-8_DEL > col1-8_DEL_plus_down 
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$4-500,$6,$7,$8}' col1-8_DEL | awk -F '\t' '{OFS=FS}{print $1, $2, $3, $5,$4,$6,$7,$8}' | sed "s/+/-/g" > col1-8_DEL_minus_up 
paste col1-8_DEL_plus_down col9_DEL | sed "s/\";/_plus_down\";/g" > all_DEL_regions.gtf 
paste col1-8_DEL_minus_up col9_DEL | sed "s/\";/_minus_up\";/g" >> all_DEL_regions.gtf 
sort -t 'tabby' -nk 4 all_DEL_regions.gtf > all_DEL_regions2.gtf 
mv all_DEL_regions2.gtf all_DEL_regions.gtf

#CON:

grep "CON" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/,/_/g" | sed "s/^/transcript_id \"SOMETHING_/g" | sed "s/$/\";gene_id \"SOMETHING_/g" | awk -F "\"" '{OFS=FS}{print $1,$2,$3,$2}' | sed "s/$/\";/g" > col9_CON 
grep "CON" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/^/REL606\t.\tCDS\t/g" | awk -F ',' '{OFS=FS}{print $1,$2,$2+$3}' | sed "s/CON,//g" | sed "s/,/\t/g" | sed "s/$/\t.\t+\t0/g" > col1-8_CON 
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+500,$6,$7,$8}' col1-8_CON > col1-8_CON_plus_down 
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$4-500,$6,$7,$8}' col1-8_CON | awk -F '\t' '{OFS=FS}{print $1, $2, $3, $5,$4,$6,$7,$8}' | sed "s/+/-/g" > col1-8_CON_minus_up 
paste col1-8_CON_plus_down col9_CON | sed "s/\";/_plus_down\";/g" > all_CON_regions.gtf 
paste col1-8_CON_minus_up col9_CON | sed "s/\";/_minus_up\";/g" >> all_CON_regions.gtf 
sort -t 'tabby' -nk 4 all_CON_regions.gtf > all_CON_regions2.gtf 
mv all_CON_regions2.gtf all_CON_regions.gtf

#SUB:

grep "SUB" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/,/_/g" | sed "s/^/transcript_id \"SOMETHING_/g" | sed "s/$/\";gene_id \"SOMETHING_/g" | awk -F "\"" '{OFS=FS}{print $1,$2,$3,$2}' | sed "s/$/\";/g" > col9_SUB
grep "SUB" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/^/REL606\t.\tCDS\t/g" | awk -F ',' '{OFS=FS}{print $1,$2,$2+$3}' | sed "s/SUB,//g" | sed "s/,/\t/g" | sed "s/$/\t.\t+\t0/g" > col1-8_SUB
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+500,$6,$7,$8}' col1-8_SUB > col1-8_SUB_plus_down
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$4-500,$6,$7,$8}' col1-8_SUB | awk -F '\t' '{OFS=FS}{print $1, $2, $3, $5,$4,$6,$7,$8}' | sed "s/+/-/g" > col1-8_SUB_minus_up
paste col1-8_SUB_plus_down col9_SUB | sed "s/\";/_plus_down\";/g" > all_SUB_regions.gtf  
paste col1-8_SUB_minus_up col9_SUB | sed "s/\";/_minus_up\";/g" >> all_SUB_regions.gtf  
sort -t 'tabby' -nk 4 all_SUB_regions.gtf > all_SUB_regions2.gtf  
mv all_SUB_regions2.gtf all_SUB_regions.gtf

#AMP:

grep "AMP" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/,/_/g" | sed "s/^/transcript_id \"SOMETHING_/g" | sed "s/$/\";gene_id \"SOMETHING_/g" | awk -F "\"" '{OFS=FS}{print $1,$2,$3,$2}' | sed "s/$/\";/g" > col9_AMP 
grep "AMP" SOMETHING_mutations.tsv | sed "s/\t/,/g" | sed "s/^/REL606\t.\tCDS\t/g" | awk -F ',' '{OFS=FS}{print $1,$2,$2+$3}' | sed "s/AMP,//g" | sed "s/,/\t/g" | sed "s/$/\t.\t+\t0/g" > col1-8_AMP
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+500,$6,$7,$8}' col1-8_AMP > col1-8_AMP_plus_down 
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$4-500,$6,$7,$8}' col1-8_AMP | awk -F '\t' '{OFS=FS}{print $1, $2, $3, $5,$4,$6,$7,$8}' | sed "s/+/-/g" > col1-8_AMP_minus_up 
paste col1-8_AMP_plus_down col9_AMP | sed "s/\";/_plus_down\";/g" > all_AMP_regions.gtf   
paste col1-8_AMP_minus_up col9_AMP | sed "s/\";/_minus_up\";/g" >> all_AMP_regions.gtf   
sort -t 'tabby' -nk 4 all_AMP_regions.gtf > all_AMP_regions2.gtf
mv all_AMP_regions2.gtf all_AMP_regions.gtf

cat all_SNP_regions.gtf all_SUB_regions.gtf all_INS_regions.gtf all_DEL_regions.gtf all_CON_regions.gtf all_AMP_regions.gtf all_MOBs.gtf | sort -t 'tabby' -nk 4 > SOMETHING_mutation_adj_regions.gtf

rm col*
rm all*

##

sed "s/.gd//g" filenames > SOMETHINGS
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" basic_code.sh >> final_code.sh/g" > interim_code.sh 
bash interim_code.sh 
sed -i "s/tabby/\\t/g" final_code.sh
bash final_code.sh
rm *.sh

#manually remove records that start with values <0, or replace them with "1".
#Check to see which files contain such records:

cat *.gtf | sort -k4,4n | less

####Step-2(a): Convert the .gtf files to .faa files using gffread####

ls *.gtf | cut -f 1-3 -d '_' > filenames
sed "s/^/..\/..\/..\/..\/..\/RNAseq_denovo\/tools\/gffread-0.12.7.Linux_x86_64\/gffread \-E \-w /g" filenames | sed "s/$/_mutation_adj_regions.faa \-g ..\/..\/..\/all_genomes\/REL606.fasta/g" | awk -F ' ' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$4}' | sed "s/$/t/g" | sed "s/_mutation_adj_regions.faat/.fasta/g" | sed "s/genomes\/ /genomes\//g" | awk -F ' ' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$4}' | sed "s/$/t/g" | sed "s/.faat/.gtf/g" | bash

#############NEXT TIME: STEP 2B ONWARDS!!!################

