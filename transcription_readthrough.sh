#Step-2 (post-gmapping):

#Identify bad maps, multiple maps, and deletions for each evolved genome. Convert remainder to gtf

grep "^>" ../REL606_transterm_motifs.faa | tr -d "^>" | cut -f 1 -d ' ' | sed "s/^/=/g" | sed "s/$/.mrna/g" | sort -u > all_transterm_motif_names.txt

grep "mrna1" something*gff3 | awk -F '\t' '($3=="mRNA")' | sed "s/identity=/identity;/g" | sed "s/coverage=/coverage;/g" | awk -F ';' '($7<80||$9<80)' | cut -f 2 -d "=" | cut -f 1 -d '.' | sed "s/^/=/g" | sed "s/$/.mrna/g" > something_bad_maps.txt
grep "mrna2" something*gff3 | awk -F '\t' '($3=="mRNA")' | cut -f 2 -d "=" | cut -f 1 -d '.' | sed "s/^/=/g" | sed "s/$/.mrna/g" | sort -u > something_multi_maps.txt
cat something_bad_maps.txt something_multi_maps.txt | sort -u > something_bad_multi_maps
grep -v -f something_bad_multi_maps something*gff3 | awk -F '\t' '($3=="mRNA")' | cut -f 1 -d ';' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,$8,"gene_id \""$9,"\";transcript_id \""$9,"\";"}' | sed "s/ID=//g" | sed "s/\t\";/\";/g" > something_transterm_100offset.gtf
awk -F '\t' '($3=="mRNA")' something*gff3 | cut -f 2 -d '=' | cut -f 1 -d '.' | sed "s/^/=/g" | sed "s/$/.mrna/g" | sort -u > all_something_mapped_motifs
grep -v -f all_something_mapped_motifs all_transterm_motif_names.txt > something_deleted_motifs.txt

ls Ara*gff3 | cut -f 1-3 -d "_" > somethings

sed "s/^/sed \"s\/something\//g" somethings | sed "s/$/\/g\" basic_code_gff3processing.sh >> final_code_gff3processing.sh/g" > interim_code_gff3processing.sh 
bash interim_code_gff3processing.sh
bash final_code_gff3processing.sh

#Step-3 (make gtfs out of applied-gd files):

grep -v "#" something*applied.gd | awk -F "applied_start=" '{print $2}' | cut -f 1 > something_col3
grep -v "#" something*applied.gd | awk -F "applied_end=" '{print $2}' | cut -f 1 > something_col4
grep -v "#" something*applied.gd | awk -F '\t' '{OFS="_"}{print $2,$1}' > something_col9
paste something_col3 something_col4 something_col9 | sed "s/^/REL606\t.\tCDS\t/g" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,".\t+\t.\tgene_id \"",$6,"\";transcript_id \"",$6,"\";"}' | sed "s/\t\";/\";/g" | sed "s/id \"\t/id \"/g" > something_applied.gtf
rm something_col*

sed "s/^/sed \"s\/something\//g" somethings | sed "s/$/\/g\" basic_code.sh >> final_code.sh/g" > interim_code.sh
bash interim_code.sh
bash final_code.sh

#Step-4: Intersect between the two gtf files to determine the mutated motifs

bedtools intersect -a ../step2_gmapping_to_evolved/something*offset.gtf -b ../step3_appliedgd_to_appliedgtf/something_applied.gtf > something_mutated_motifs.gtf

sed "s/^/sed \"s\/something\//g" somethings | sed "s/$/\/g\" basic_code_intersect.sh >> final_code_intersect.sh/g" > interim_code_intersect.sh
bash interim_code_intersect.sh
bash final_code_intersect.sh

