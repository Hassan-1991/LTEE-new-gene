#get genome length first

awk '!/^>/ { printf "%s", $0; n="\n"} /^>/ {print n $0; n = "" } END {printf "%s", n }' REL606.fasta | grep -v "^>" | wc

#The second number below is the size of the sliding window
#The third number is genome length, extracted from the first step

seq 1 400 4629813 > col1
seq 1 400 4629813 | tail -n+2 > col2
paste col1 col2 | cat -n | sed "s/ //g" | awk -F '\t' '{OFS=FS}{print "REL606\t.\tCDS",$2,$3-1,".\t+\t0\tgene_id \"window_plus_",$1,"\";transcript_id \"window_plus_",$1,"\";"}' | head -n-1 | sed "s/_\t/_/g" | sed "s/\t\"/\"/g" > REL606_slidwindow400.gtf
paste col1 col2 | cat -n | sed "s/ //g" | awk -F '\t' '{OFS=FS}{print "REL606\t.\tCDS",$2,$3-1,".\t-\t0\tgene_id \"window_minus_",$1,"\";transcript_id \"window_minus_",$1,"\";"}' | head -n-1 | sed "s/_\t/_/g" | sed "s/\t\"/\"/g" >> REL606_slidwindow400.gtf
rm col*

#Use this gtf file to extract corresponding sequences from the genome

gffread -E -w REL606_slidwindow400_CDS.faa -g REL606.fasta REL606_slidwindow400.gtf

#Generate gmap indices for each line to which these sequences will be mapped

-------------

#Map to evolved genomes (indices should be in the same folder)

echo 'gmap -D . -d SOMETHING -f 2 --gff3-fasta-annotation=1 REL606_slidwindow400_CDS.faa > SOMETHING_slidwindow400_CDS.gff3' > basic_code_400.sh 
ls -d Ara* | grep -v "gff3" > SOMETHINGS 
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" basic_code_400.sh >> final_code_400.sh/g" > interim_code_400.sh 
bash interim_code_400.sh 
bash final_code_400.sh

#This returns a bunch of gff3 files. Only keep the windows that are common and map well to all lines and strains.
#Get rid of all sequences that map more than once in any of the evolved genomes - they are duplicates
#Also get rid of sequences that show an extended or truncated map in the evolved genomes, these are either spurious maps or includes insertion sequences

echo 'grep "mRNA" SOMETHING_slidwindow400_CDS.gff3 | grep "path2" | cut -f 2 -d '=' | cut -f 1 -d '.' | sed "s/^/=/g" | sed "s/$/;/g" > SOMETHING_400_doublemaps 
grep -v -f SOMETHING_400_doublemaps SOMETHING_slidwindow400_CDS.gff3 | grep "mRNA" > SOMETHING_400_step1 
awk -F '\t' '(($5-$4>410)||($5-$4<390))' SOMETHING_400_step1 > SOMETHING_400_longshortmaps 
grep -v -f SOMETHING_400_longshortmaps SOMETHING_400_step1 > SOMETHING_400_step2' > basic_code.sh
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" basic_code.sh >> final_code.sh/g" > interim_code.sh 
bash interim_code.sh 
bash final_code.sh 
rm *longshortmaps 
rm *doublemaps 

#check to see if the remaining "common sequences" map with a suspiciously low identity or query cover by scanning their ID and qCov values: 

cat *step2 | grep "mRNA" | grep "path1" | sed "s/=/;/g" | cut -f 12-15 -d ';' | sort -u | sort -t ';' -nk 2 | less 
#lowest coverage - 98 for 400
#lowest identity - 78.1
#the remanining maps look file

#Finally, get a set of windows that mapped to all strains, and weren't deleted in any of them. This represents the set of "common windows" we'll be using for analysis

cat Ara*400*step2 | cut -f 2,9 | cut -f 1 -d ';' | sed "s/.mrna1//g" | sed "s/ID=//g" | cut -f 2 | sort | uniq -c | sed "s/ //g" | grep "25window" | sed "s/25window/window/g" | sed "s/^/=/g" | sed "s/$/;/g" > common_400bp_windows.txt 

#The "25" number above in the previous line specifies - only select those that appear 25 times (for all 25 strains) in the record.

echo 'grep -f common_400bp_windows.txt SOMETHING_400_step2 > SOMETHING_400_commonwindows.gff3' > basic_code.sh 
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" basic_code.sh >> final_code.sh/g" > interim_code.sh 
bash interim_code.sh 
bash final_code.sh 

#Now to convert these gff3 files into gtfs 

echo 'grep "mRNA" SOMETHING_400_commonwindows.gff3 | sed "s/ID=/transcript_id=\"/g" | sed "s/.mrna1;Name=/\";gene_id \"/g" | sed "s/;Parent/\";Parent/g" | sed "s/\tmRNA\t/\tCDS\t/g" | awk -F '\t' '{OFS=FS}{print $1,$2=".",$3,$4,$5,$6,$7,$8,$9}' > SOMETHING_400_commonwindows.gtf' > basic_code.sh
sed "s/^/sed \"s\/SOMETHING\//g" SOMETHINGS | sed "s/$/\/g\" basic_code.sh >> final_code.sh/g" > interim_code.sh 
bash interim_code.sh 
bash final_code.sh 

#Generate the feature lists for two ancestors manually

sed "s/=/\"/g" common_400bp_windows.txt | sed "s/;/\"/g" > common_400bp_windows_2.txt 
grep -f common_400bp_windows_2.txt REL606_slidwindow400.gtf > REL606_400_commonwindows.gtf 
cp REL606_400_commonwindows.gtf REL607_400_commonwindows.gtf
