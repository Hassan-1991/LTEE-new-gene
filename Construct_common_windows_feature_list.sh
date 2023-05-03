#This code describes how to make the transcriptome comprised of 400-bp windows

#Step-1: Make an initial file with all 400-bp windows in the genome

awk '!/^>/ { printf "%s", $0; n="\n"} /^>/ {print n $0; n = "" } END {printf "%s", n }' REL606.fasta | grep -v "^>" | wc

#The value this returns is 4629813. We use this and the window size, 400, to write the next piece of code
#Construct a gtf file with all 400-bp windows across both strands of the genome
#Get the corresponding sequences

seq 1 400 4629813 > col1 #Start coordinate of each window
seq 1 400 4629813 | tail -n+2 > col2 #Stop coordinate of each window
paste col1 col2 | cat -n | sed "s/ //g" | awk -F '\t' '{OFS=FS}{print "REL606\t.\tCDS",$2,$3-1,".\t+\t0\tgene_id \"window_plus_",$1,"\";transcript_id \"window_plus_",$1,"\";"}' | head -n-1 | sed "s/_\t/_/g" | sed "s/\t\"/\"/g" > REL606_slidwindow400.gtf #Windows on plus strand
paste col1 col2 | cat -n | sed "s/ //g" | awk -F '\t' '{OFS=FS}{print "REL606\t.\tCDS",$2,$3-1,".\t-\t0\tgene_id \"window_minus_",$1,"\";transcript_id \"window_minus_",$1,"\";"}' | head -n-1 | sed "s/_\t/_/g" | sed "s/\t\"/\"/g" >> REL606_slidwindow400.gtf #Windows on minus strand
rm col* #remove intermediate files
gffread -E -w REL606_slidwindow400_CDS.faa -g REL606.fasta REL606_slidwindow400.gtf #get their sequences

#Step-2: From this file, get rid of all windows that map more than once in any genome, experienced insertions or deletions, diverged substantially, or got deleted from any evolved genome

#We first map the window sequences to each evolved genome with gmap

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6;
do gmap_build -D . -d "$i" "$i".faa #build gmap indices for each evolved genome
gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 REL606_slidwindow400_CDS.faa > "$i"_slidwindow400_CDS.gff3 #map to evolved genomes
grep "mRNA" "$i"_slidwindow400_CDS.gff3 | grep "path2" | cut -f 2 -d '=' | cut -f 1 -d '.' | sed "s/^/=/g" | sed "s/$/;/g" > "$i"_400_doublemaps #identify all windows that map more than once in any evolved genome
grep -v -f "$i"_400_doublemaps "$i"_slidwindow400_CDS.gff3 | grep "mRNA" > "$i"_400_step1 #subtract "doublemaps"
awk -F '\t' '(($5-$4>410)||($5-$4<390))' "$i"_400_step1 > "$i"_400_longshortmaps #identify the windows that map with a length +/-10bp
grep -v -f "$i"_400_longshortmaps "$i"_400_step1 > "$i"_400_step2; done #Subtract "longshortmaps"
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
do grep -f common_400bp_windows.txt "$i"_400_step2grep "mRNA" "$i"_400_commonwindows.gff3 | sed "s/ID=/transcript_id=\"/g" | sed "s/.mrna1;Name=/\";gene_id \"/g" | sed "s/;Parent/\";Parent/g" | sed "s/\tmRNA\t/\tCDS\t/g" | awk -F '\t' '{OFS=FS}{print $1,$2=".",$3,$4,$5,$6,$7,$8,$9}' > "$i"_400_commonwindows.gtf; done

#Generate the feature lists for two ancestors manually

sed "s/=/\"/g" common_400bp_windows.txt | sed "s/;/\"/g" > common_400bp_windows_2.txt 
grep -f common_400bp_windows_2.txt REL606_slidwindow400.gtf > REL606_400_commonwindows.gtf 
cp REL606_400_commonwindows.gtf REL607_400_commonwindows.gtf
