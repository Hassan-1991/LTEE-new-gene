#Step-1: Predict both intrinsic and rho-dependent terminators in ancestral LTEE genome.
#For Rho-dependent, use the RhoTermPredict.py script found here - https://github.com/MarcoDiSalvo90/RhoTermPredict/blob/master/RhoTermPredict_algorithm.py

python3 RhoTermPredict_algorithm.py

#Input genome name when prompted.
#For intrinsic, use ARNold web server:

http://rssf.i2bc.paris-saclay.fr/toolbox/arnold/

#Step-2: Convert these raw files into gtfs:

#ARNold:

tail -n+2 Arnold_rawoutput_REL606.txt | 
sed "s/ Rnamotif //g" | 
sed "s/ Both //g" | 
sed "s/ Erpin //g" | 
rev | 
cut -f 2- -d " " | 
rev | 
awk '{OFS=""}{print "REL606\t.\tCDS\t",$1,"\t",$1+length($3),"\t",$2}' | 
cat -n | 
sed "s/ REL606/\tREL606/g" | 
awk '{OFS=""}{print $2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",".\t",$7,"\t.\tgene_id \"Intrinsic_",$1,"\";transcript_id \"Intrinsic_",$1}' | 
sed "s/$/\";/g" > ARNold_REL606.gtf

#RhoTermPredict:
#This csv file has a weird ^M character after "plus"/"minus" that's apparently impossible to get rid of, even with sed shenanigans
#Just making the gtf without it, then adding a separate column with + and -

tail -n+2 predictions_coordinates_REL606.csv | 
sed "s/^/Rhodep_/g" | 
sed "s/\tplus/\t+/g" | 
sed "s/\tminus/\t-/g" | 
awk -F '\t' '{OFS=FS}{print "REL606\t.\tCDS",$2,$3,".\t.",$1}' | 
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$7,"gene_id \"",$8,"\";transcript_id \"",$8,"\";"}' | 
sed "s/\tRho/Rho/g" | 
sed "s/\t\";/\";/g" > interim

tail -n+2 predictions_coordinates_REL606.csv | 
grep "minus" | 
sed "s/^/-\t/g" | 
cut -f 1 >> interim_2

paste interim interim_2 | 
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6,$9,$7,$8}' | 
sed "s/_T/_/g" > RhoTerm_REL606.gtf

rm interim*

#Cat 'em together and sort by start coordinates

cat ARNold_REL606.gtf RhoTerm_REL606.gtf | sort -t '        ' -nk 4 > REL606_transterm_motifs.gtf

#GETFLANKS

bedtools flank -b 100 -s -i REL606_transterm_motifs.gtf -g REL606_genomefile.txt | awk '{ printf "%d.%d\t%s\n",(!a[$10]++? ++c:c),a[$10],$0 }' | cut -f 2- -d '.' | cut -f 1 -d ';' | awk -F '\t' '{OFS=FS}{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$1}' | sed "s/\"\t/_/g" | sed "s/$/\";/g" > REL606_transterm_motifs_100bpflanks.gtf

#GET SEQS

gffread

#Gmap to evolved, 50Ks first

#Next step: Get all cases where the 100-bp flanks map well and once in each genome.

#1-2:

awk -F '\t' '($3=="mRNA")' Ara-6_50000gen_11389_transtermmotifs_100bpflanks.gff3 | grep "mrna2" | cut -f 1 -d ';' | cut -f 2- -d '=' | sed "s/\.mrna2//g" | sed "s/^/=/g" | sed "s/$/\.mrna/g" > Ara-6_doublemaps
grep -v -f Ara-6_doublemaps Ara-6_50000gen_11389_transtermmotifs_100bpflanks.gff3 | awk -F '\t' '($3=="mRNA")' > Ara-6_50000gen_11389_transtermmotifs_100bpflanks_singlemaps.gff3

#2-3:

sed "s/=/;/g" Ara-6_50000gen_11389_transtermmotifs_100bpflanks_singlemaps.gff3 | awk -F ';' '($13<80||$15<80)' | cut -f 2 -d ';' | sed "s/\.mrna1//g" | sort -u | sed "s/^/=/g" | sed "s/$/\.mrna/g" > Ara-6_badmaps
grep -v -f Ara-6_badmaps Ara-6_50000gen_11389_transtermmotifs_100bpflanks_singlemaps.gff3 | awk -F '\t' '($3=="mRNA")' > Ara-6_50000gen_11389_transtermmotifs_100bpflanks_singlemaps_goodmaps.gff3

#3-4:

cut -f 1 -d ';' Ara-6_50000gen_11389_transtermmotifs_100bpflanks_singlemaps_goodmaps.gff3 | cut -f 9 | cut -f 2- -d '=' | cut -f 1 -d '.' | rev | cut -f 2- -d "_" | rev | sort | uniq -c | grep " 2 " | rev | cut -f 1 -d ' ' | rev | awk '{OFS=""}{print $1,"_1 "$1,"_2"}' | sed "s/ /\n/g" | sed "s/^/=/g" | sed "s/$/\.mrna/g" > Ara-6_singlemaps_goodmaps_bothmaps
grep -f Ara-6_singlemaps_goodmaps_bothmaps Ara-6_50000gen_11389_transtermmotifs_100bpflanks_singlemaps_goodmaps.gff3 > Ara-6_50000gen_11389_transtermmotifs_100bpflanks_singlemaps_goodmaps_bothmaps.gff3

#4-5:

cut -f 1,2 -d "_" Ara-6_singlemaps_goodmaps_bothmaps | sort -u | sed "s/$/_/g" > Ara-6_singlemaps_goodmaps_bothmaps.txt
sed "s/$/1/g" Ara-6_singlemaps_goodmaps_bothmaps.txt > Ara-6_singlemaps_goodmaps_bothmaps_upstream.txt
sed "s/$/2/g" Ara-6_singlemaps_goodmaps_bothmaps.txt > Ara-6_singlemaps_goodmaps_bothmaps_downstream.txt

grep -f Ara-6_singlemaps_goodmaps_bothmaps_upstream.txt Ara-6_50000gen_11389_transtermmotifs_100bpflanks_singlemaps_goodmaps_bothmaps.gff3 | cut -f 1 -d ';' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$6,$7,$8,$9}' | sed "s/_1\.mrna1/_interval_Ara-6;/g" > test1
grep -f Ara-6_singlemaps_goodmaps_bothmaps_downstream.txt Ara-6_50000gen_11389_transtermmotifs_100bpflanks_singlemaps_goodmaps_bothmaps.gff3 | cut -f 1 -d ';' | cut -f 4 > test2
paste test1 test2 | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4+1,$9-1,$5,$6,$7,$8}' > Ara-6_intervals.gtf

#What are intervals of "bad lengths"? The ideal approach would be to impose different ranges for arnold and rhoterm; but for my first volley I'm imposing a general cutoff of 30 to 100.

awk -F '\t' '{OFS=FS}{print $0,$5-$4}' Ara-6_intervals.gtf | awk -F '\t' '($10>30&&$10<100)' | cut -f 1-9 > goodintervals_interim

gffread -E -w goodintervals_interim.faa -g ../../../../all_genomes/Ara-6_50000gen_11389.fasta goodintervals_interim

#Another step (make gtfs out of applied-gd files):

grep -v "#" something*applied.gd | awk -F "applied_start=" '{print $2}' | cut -f 1 > something_col3
grep -v "#" something*applied.gd | awk -F "applied_end=" '{print $2}' | cut -f 1 > something_col4
grep -v "#" something*applied.gd | awk -F '\t' '{OFS="_"}{print $2,$1}' > something_col9
paste something_col3 something_col4 something_col9 | sed "s/^/REL606\t.\tCDS\t/g" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,".\t+\t.\tgene_id \"",$6,"\";transcript_id \"",$6,"\";"}' | sed "s/\t\";/\";/g" | sed "s/id \"\t/id \"/g" > something_applied.gtf
rm something_col*

#Big problem: It doesn't seem like ARNOLD can detect termination signals based on small pieces. It needs the whole genome

#Let's pause this here. Neither tool can do this.
#We need to start this transcription-first. Would be good to do it mutation first but...
#General way of doing it would be to look at the downstream regions of all cases with mutations.
#But I only care about the 3' ends of genes though.

