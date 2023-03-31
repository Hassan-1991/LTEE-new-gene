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

####Step-2(b) map these to evolved populations and make gtf files unique to each population####

#gmap code to map these sequences to their evolved counterparts. Must be run in a directory that already has the indices

for i in Ara+1_50000gen_11392 Ara+2_50000gen_11342 Ara+3_50000gen_11345 Ara+4_50000gen_11348 Ara+5_50000gen_11367 Ara-1_50000gen_11330 Ara-2_50000gen_11333 Ara-3_50000gen_11364 Ara-4_50000gen_11336 Ara-5_50000gen_11339 Ara-6_50000gen_11389; do gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 "$i"_mutation_adj_regions.faa > "$i"_mutation_adj_regions.gff3; done

#Get rid of all cases where a sequence maps more than once. Skeletal code:

grep "mRNA" SOMETHING_mutation_adj_regions.gff3 | grep "path2" | cut -f 1 -d ';' | cut -f 9 | sed "s/ID=//g" | sed "s/.mrna2//g" > SOMETHING_duplicates.txt
grep -v -f SOMETHING_duplicates.txt SOMETHING_mutation_adj_regions.gff3 | grep "path1" | grep "mRNA" | sed "s/mRNA/CDS/g" | sed "s/ID=/transcript_id \"/g" | sed "s/.mrna1;Name=/\";gene_id \"/g" | sed "s/;Parent/\";Parent/g"| cut -f 1,2 -d ";" | sed "s/$/;/g" > SOMETHING_mutation_adj_regions.gtf
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" 4b_basic_code.sh >> 4b_final_code.sh/g" > 4b_interim_code.sh
bash 4b_interim_code.sh
bash 4b_final_code.sh

#Add in upstream control strands for each case:

for i in Ara+1_50000gen_11392 Ara+2_50000gen_11342 Ara+3_50000gen_11345 Ara+4_50000gen_11348 Ara+5_50000gen_11367 Ara-1_50000gen_11330 Ara-2_50000gen_11333 Ara-3_50000gen_11364 Ara-4_50000gen_11336 Ara-5_50000gen_11339 Ara-6_50000gen_11389; do
sed "s/+/%/g" ../step3_gmapgff3_to_gtf/"$i"_mutation_adj_regions.gtf | sed "s/\t-\t/\t+\t/g" | sed "s/%/-/g" | sed "s/id \"/id \"control_/g" >> "$i"_mutation_adj_regions.gtf

#For some reason, these doubled the control number - so just sort -u all files.

for file in *.gtf; do sort -u -o $file $file; done

#Get rid of large maps (>600):

awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10=$5-$4}' SOMETHING_mutation_adj_regions.gtf | awk -F '\t' '($10<600)' > SOMETHING_mutation_adj_regions_4b.gtf

sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" 4b_basic_code.sh >> 4b_final_code.sh/g" > 4b_interim_code.sh
bash 4b_interim_code.sh
bash 4b_final_code.sh

#Generate corresponding ancestor file for each population, which contains just the final mutations

#tag initial gtfs with the ancestor name

for file in *.gtf; do mv $file ancestor_$file; done

grep -v "control" SOMETHING_mutation_adj_regions_4b.gtf | cut -f 2 -d "\"" > SOMETHING_IDs.txt
grep -f SOMETHING_IDs.txt ancestor_SOMETHING_mutation_adj_regions.gtf | sort -u > ancestor_SOMETHING_mutation_adj_regions_4b.gtf
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" 4b_basic_code.sh >> 4b_final_code.sh/g" > 4b_interim_code.sh
bash 4b_interim_code.sh
bash 4b_final_code.sh

#Now top up these ancestral 4b files with controls

sed "s/\t+\t/\t%\t/g" ancestor_SOMETHING_mutation_adj_regions_4b.gtf |  sed "s/\t-\t/\t+\t/g" | sed "s/%/-/g" | sed "s/id \"/id \"control_/g" | sort -u > ancestor_SOMETHING_mutation_adj_regions_5.gtf
cat ancestor_SOMETHING_mutation_adj_regions_4b.gtf >> ancestor_SOMETHING_mutation_adj_regions_5.gtf
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" 4b_basic_code.sh >> 4b_final_code.sh/g" > 4b_interim_code.sh
bash 4b_interim_code.sh
bash 4b_final_code.sh

for file in ancestor*_5.gtf; do sort -u -o $file $file; done

more fixes:

awk -F '\t' '{OFS=FS}($2=".")' SOMETHING*.gtf | cut -f 1-9 > SOMETHING_fixed.gtf
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" 4b_basic_code.sh >> 4b_final_code.sh/g" > 4b_interim_code.sh
bash 4b_interim_code.sh
bash 4b_final_code.sh

#The transcriptomes should consist not just of the mutation-adjacent regions, but also the annotated genes in each line (that don't overlap with these regions)

../../../../../RNAseq_denovo/tools/bedtools intersect -a ../../../insremoved_gtf/SOMETHING*gtf -b SOMETHING*gtf -v -s | cut -f 2 -d "\"" | sort -u > SOMETHING_coding_ids.txt
../../../../../RNAseq_denovo/tools/bedtools intersect -a ../../../insremoved_gtf/REL606_insremoved.gtf -b ancestor_SOMETHING*gtf -v -s | cut -f 2 -d "\"" | sort -u  > ancestor_SOMETHING_coding_ids.txt
cat ancestor_SOMETHING_coding_ids.txt SOMETHING_coding_ids.txt | sort | uniq -c | sed "s/ /,/g" | sed "s/,,,,,,//g" | grep "2," | cut -f 2 -d ',' | sed "s/^/\"/g" | sed "s/$/\"/g" > SOMETHING_common_IDs.txt
grep -f SOMETHING_common_IDs.txt ../../../insremoved_gtf/SOMETHING*gtf | sort -u >> SOMETHING_mutation_adj_regions_4b.gtf
grep -f SOMETHING_common_IDs.txt ../../../insremoved_gtf/REL606_insremoved.gtf | sort -u >> ancestor_SOMETHING*gtf

#manually fix in case of Ara-3, only pick out 50k gen

bash 4b_interim_code.sh
bash 4b_final_code.sh

#For Ara-3:

../../../../../RNAseq_denovo/tools/bedtools intersect -a ../../../insremoved_gtf/REL606_insremoved.gtf -b REL606_33k_mutation.gtf -v -s | cut -f 2 -d "\"" | sort -u | sed "s/^/\"/g" | sed "s/$/\"/g" > coding_ids.txt
grep -f coding_ids.txt ../../../insremoved_gtf/SOMETHING_insremoved.gtf >> SOMETHING_33k_mutation.gtf
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" 4b_basic_code.sh >> 4b_final_code.sh/g" > 4b_interim_code.sh
bash 4b_interim_code.sh
bash 4b_final_code.sh

#Run HTSEQ

ls ancestor_SOMETHING*gtf >> SOMETHING_test1
ls ancestor_SOMETHING*gtf >> SOMETHING_test1
ls ancestor_SOMETHING*gtf >> SOMETHING_test1
ls ancestor_SOMETHING*gtf >> SOMETHING_test1
ls ancestor_SOMETHING*gtf >> SOMETHING_test1
ls ancestor_SOMETHING*gtf >> SOMETHING_test1
ls ancestor_SOMETHING*gtf >> SOMETHING_test1
ls ancestor_SOMETHING*gtf >> SOMETHING_test1
ls SOMETHING*.gtf >> SOMETHING_test1
ls SOMETHING*.gtf >> SOMETHING_test1
ls SOMETHING*.gtf >> SOMETHING_test1
ls SOMETHING*.gtf >> SOMETHING_test1
ls ../../../../50k_bowtie_indices_bamfiles/*R*_sorted.bam | sort -t '_' -nk 4 > SOMETHING_test2
ls ../../../../50k_bowtie_indices_bamfiles/*SOMETHING_sorted.bam | sort -t '_' -nk 4 >> SOMETHING_test2
paste SOMETHING_test2 SOMETHING_test1 | sed "s/^/time htseq-count -f bam -a 0 -t CDS --nonunique all /g" | sed "s/\t/ /g" | sed "s/$/ \| head -n -5 /g" > SOMETHING_test3
rev SOMETHING_test2 | cut -f 1 -d '/' | rev | cut -f 1 -d '_' | sed "s/-/,/g" | sed "s/am//g" | sed "s/ap//g" > SOMETHING_test4
paste SOMETHING_test3 SOMETHING_test4 | sed "s/\t/ \| sed \"s\/\^\//g" | sed "s/$/,\/g\"/g" | sed "s/$/ \| sed \"s\/TABBY\/\,\/g\" >> SOMETHING_htseq_50K_mutation_count.csv/g" > SOMETHING_htseq_mutation_code.sh

sed -i "s/TABBY/\t/g" 4b_basic_code.sh

sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" 4b_basic_code.sh >> 4b_final_code.sh/g" > 4b_interim_code.sh
bash 4b_interim_code.sh
bash 4b_final_code.sh

#after running:

for file in *.csv; do sed -i "1s/^/repl,seqtype,line,target_id,est_counts\n/g" $file; done
for file in *.csv; do sed -i "s/AraR6/REL606/g" $file; done
for file in *.csv; do sed -i "s/AraR7/REL607/g" $file; done

#Run DESeq2 as in the other code

tail -n+2 SOMETHING_htseq_50K_mutation_foldchanges.csv | cut -f 3 -d ',' | sort -u | sed "s/^/\"/g" | sed "s/$/\"/g" > SOMETHING_test1 
grep -v -f SOMETHING_test1 SOMETHING_mutation_adj_regions_4b.gtf | cut -f 2 -d "\"" | sed "s/^/ribo,SOMETHING,/g" | sed "s/$/,0,NA,NA,NA,NA/g" >> SOMETHING_htseq_50K_mutation_foldchanges.csv 
grep -v -f SOMETHING_test1 SOMETHING_mutation_adj_regions_4b.gtf | cut -f 2 -d "\"" | sed "s/^/rna,SOMETHING,/g" | sed "s/$/,0,NA,NA,NA,NA/g" >> SOMETHING_htseq_50K_mutation_foldchanges.csv 
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" FCcleanup_basic_code.sh >> FCcleanup_final_code.sh/g" > FCcleanup_interim_code.sh 
bash FCcleanup_interim_code.sh 
bash FCcleanup_final_code.sh

#for 33k:

#Get all founds:

cut -f 2 -d "\"" SOMETHING_33k_mutation_fixed.gtf | sort -u > SOMETHING_found
grep -v -f SOMETHING_found all | sed "s/^/SOMETHING,BiolRep1,/g" | sed "s/$/,0/g" >> test
grep -v -f SOMETHING_found all | sed "s/^/SOMETHING,BiolRep2,/g" | sed "s/$/,0/g" >> test
grep -v -f SOMETHING_found all | sed "s/^/SOMETHING,BiolRep3,/g" | sed "s/$/,0/g" >> test

