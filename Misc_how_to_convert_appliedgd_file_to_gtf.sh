#step-1: get the start and end positions

grep -v "#" Ara+1_50000gen-applied.gd | awk -F "applied_start=" '{print $2}' | cut -f 1 > col3
grep -v "#" Ara+1_50000gen-applied.gd | awk -F "applied_end=" '{print $2}' | cut -f 1 > col4

#step-2: get the name of each mutation (type+serial)

grep -v "#" Ara+1_50000gen-applied.gd | awk -F '\t' '{OFS="_"}{print $2,$1}' > col9

#step-3: add additional information to make the gtf:

paste col4 col3 col9 | 
sed "s/^/REL606\t.\tCDS\t/g" | 
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,".\t+\t.\tgene_id \"",$6,"\";transcript_id \"",$6,"\";"}' | 
sed "s/\t\";/\";/g" | 
sed "s/id \"\t/id \"/g" > Ara+1_50000gen-applied.gtf
rm col*

#how to do this for all

#in the directory where all the applied.gd files are:

ls *applied.gd | sed "s/-applied/_applied/g" | rev | cut -f 2- -d "_" | rev | less > somethings

grep -v "#" something*applied.gd | awk -F "applied_start=" '{print $2}' | cut -f 1 > something_col3
grep -v "#" something*applied.gd | awk -F "applied_end=" '{print $2}' | cut -f 1 > something_col4
grep -v "#" something*applied.gd | awk -F '\t' '{OFS="_"}{print $2,$1}' > something_col9
paste something_col3 something_col4 something_col9 | 
sed "s/^/REL606\t.\tCDS\t/g" | 
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,".\t+\t.\tgene_id \"",$6,"\";transcript_id \"",$6,"\";"}' | 
sed "s/\t\";/\";/g" | 
sed "s/id \"\t/id \"/g" > something_applied.gtf
rm something_col*

sed "s/^/sed \"s\/something\//g" somethings | sed "s/$/\/g\" basic_code.sh >> final_code.sh/g" > interim_code.sh 
bash interim_code.sh
bash final_code.sh
