#Extract ORFs from REL606 genome using ncbi tool ORFfinder

time ORFfinder -c t -ml 30 -outfmt 1 -out REL606_ORFs_CDS.faa -in ../all_genomes/REL606.fasta

usearch -sortbylength REL606_allORFs.faa -fastaout REL606_allORFs_sorted.faa
usearch -cluster_smallmem REL606_allORFs_sorted.faa -id 1 -centroids REL606_allORFs_sorted_nr.faa

grep "^>" REL606_allORFs_sorted_nr.faa | cut -f 1 -d ' ' | tr -d ">" | sed "s/:/_/g" | awk -F '_' '{OFS=""}{print "REL606\t.\tCDS\t",$4,"\t",$5,"\t.\t+\t.\tgene_id \"",$2,"\";transcript_id \"",$2,"\";"}' > REL606_allORFs_firststep.gtf
awk -F '\t' '($5-$4>0)' REL606_allORFs_firststep.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4+1,$5+1,$6,$7,$8,$9}' > REL606_allORFs_secondstep.gtf
awk -F '\t' '($5-$4<0)' REL606_allORFs_firststep.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5+1,$4+1,$6,$7="-",$8,$9}' >> REL606_allORFs_secondstep.gtf

bedtools intersect -a REL606_allORFs_secondstep.gtf -b REL606.gtf -v | sed "s/ORF/intergenic_ORF/g" > REL606_allORFs_thirdstep.gtf
bedtools intersect -a REL606_allORFs_secondstep.gtf -b REL606_insremoved.gtf -S -wo | awk -F '\t' '($19>30)' | cut -f 1-9 | sort -u | sed "s/ORF/antisense_ORF/g" >> REL606_allORFs_thirdstep.gtf

gffread -E -w REL606_allORFs_thirdstep.faa -g REL606.fasta REL606_allORFs_thirdstep.gtf

#Gmap code skeleton (indices should be made as in Construct_common_windows_feature_list.sh. The index directory should be in the same directory where the code is run)

../../../../tools_scripts/gmap-2021-12-17/bin/gmap -D . -d Ara+1_50000gen_11392 -f 2 --gff3-fasta-annotation=1 REL606_allORFs_thirdstep.faa > Ara+1_50000gen_11392_gmapped_ORFs.gff3

#get the IDs of all ORFs that map twice in any genomes

cat *gff3 | grep "mRNA" | grep "path2" | cut -f 2 -d '=' | cut -f 1 -d '.' | sed "s/^/=/g" | sed "s/$/./g" | sort -u > doublemaps

#get IDs of ORFs that are shorter or longer than the ORF size range (92 bp - 3000 bp)

cat *gff3 | grep "mRNA" | grep "path1" | awk -F '\t' '{OFS=FS}{print $0,$5-$4}' | awk -F '\t' '($10>3000||$10<80)' | cut -f 2 -d '=' | cut -f 1 -d '.' | sed "s/^/=/g" | sed "s/$/./g" | sort -u > longshortmaps

#get IDs of ORFs that have poor pident (<90)/qcov (<90) scores

cat *gff3 | grep "mRNA" | grep "path1" | sed "s/=/;/g" | awk -F ';' '($13<90||$15<90)' | cut -f 2 -d ';' | sed "s/^/=/g" | cut -f 1 -d '.' | sed "s/$/./g" | sort -u > badmaps

cat longshortmaps badmaps doublemaps | sort -u > all_badmaps

#only get those that all 11 populations have in common

cat Ara*gff3 | grep -v -f all_badmaps | grep "mRNA" | cut -f 2- -d '=' | cut -f 1 -d '.' | sort | uniq -c | grep -v " 11 " | rev | cut -f 1 -d ' '  | rev | sed "s/^/=/g" | sed "s/$/./g" > missingmaps

#gff3togtf

for i in Ara+1 Ara+2 Ara+3 Ara+4 Ara+5 Ara-1 Ara-2 Ara-3 Ara-4 Ara-5 Ara-6; do awk -F '\t' '($3=="mRNA")' "$i"*gff3 | grep -v -f removeORFs | cut -f 1 -d ';' | sed "s/ID=//g" | awk -F '\t' '{OFS=""}{print $1,"\t",$2,"\t",$3=="CDS","\t",$4,"\t",$5,"\t",$6,"\t",$7,"\t",$8,"\t","transcript_id \"",$9,"\";gene_id \"",$9,"\";"}' | sed "s/.mrna1//g" > "$i"_commonORFs.gtf; done

#ancestors

cat Ara*common*gtf | cut -f 2 -d "\"" | s-u | sed "s/^/\"/g" | sed "s/$/\"/g" > all_common_ORFs.txt
grep -f all_common_ORFs.txt REL606_allORFs_thirdstep.gtf > REL606_commonORFs.gtf

#NEXT STEP: htseq



